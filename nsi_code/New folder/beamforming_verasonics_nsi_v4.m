

% 
% function [img_recons_final, backprojections] = dericks_DAS_beamforming(Nz, Nx, N_sensor, ix_cord, iz_cord, coord_sensor,...
%     dt, raw_data, speed)


%   ^   |
%   |   |
%   |   |
%   |   |
%   Nz  |
%   |   |
%   |   |
%   |   |
%   |   |
%   V   --------------------------------
%
%       <-------------Nx--------------->
%
%   =============INPUTS======================
%
%  Nz => # of points running horizontally [point space]
%
%  Nx => # of points running vertically [point space]
%
%  N_sensor => # of sensors
%
%  ix_cord => ['Nz' x 'Nx'] array, meshgrid of x coordinates [mm]
% 
%  iz_cord => ['Nz' x 'Nx'], array, meshgrid of z coordinates [mm]
%
%  coord_sensor => [2 x N_sensor] array, where [1 x ...] stores 'x'
%  and [2 x ...] stores 'z' [mm]
%
% =============OUTPUTS========================
%
% img_recons => ['Nz' x 'Nx'] array that stores the DAS image
%
% backprojections => ['Nz' x 'Nx' x 'N_sensor'] array that stores the
% individual backprojections of each 'N_sensor'
%
% 
%
% Execution of code:
%
%

% intm_sensor_data = raw_data';
% Uncomment below to view individual scanlines!
% function [img_recons_final] = tx_synth_DAS_beamforming_wavelengths_harmonic(test, Nz, Nx, N_sensor, ix_cord, iz_cord, coord_sensor,...
%     Fs, raw_data, f_c, size_data, element_per_line, speed)

function [img_recons_final] = beamforming_verasonics_nsi_v4(Nz, Nx, N_sensor, ix_cord, iz_cord, coord_sensor,...
    Fs, raw_data, f_c, size_data, element_per_line, speed, nsi_mode, pertubation, focus_index_z, start)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creation of Weights for NSI %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

curr_weight = ones(1,32);

switch (nsi_mode)
    case 0
        curr_weight(17:32) = -1;
    case 1
        curr_weight(17:32) = -1;
        curr_weight = curr_weight + pertubation;
    case 2
        curr_weight(1:16) = -1;
        curr_weight = curr_weight + pertubation;
    otherwise
        curr_weight = ones(1,32);
end

img_recons = zeros(Nz, Nx); % Create image
img_recons_final = zeros(Nz, Nx); % Create image pt. 2
for curr_x = 232:282 % {-2.5 [mm], 2.5 [mm]} , each scanline being 0.1 [mm] wide
    for curr_z = 1:Nz % Iterate through depth [mm]

        count = 1;
        active_elements = (curr_x-207):(curr_x-207)+31; % Declare elements for scanline
        
        for curr_rx = active_elements
            %%%%%%%%%%%%%%%%%%%%
            %% w/ no Focusing %%
            %%%%%%%%%%%%%%%%%%%%
            % tx_distance = iz_cord(curr_z,curr_x);
            % % rx_distance = sqrt((coord_sensor(1,curr_rx)-ix_cord(curr_z,curr_x)).^2 + (coord_sensor(2,curr_rx)-iz_cord(curr_z,curr_x)).^2);
            % rx_distance = sqrt((coord_sensor(curr_rx)-ix_cord(curr_z,curr_x)).^2 + (iz_cord(curr_z,curr_x)).^2);
            % 
            % path_length = tx_distance + rx_distance;
            % tau = Fs*(path_length/speed);


            %%%%%%%%%%%%%%%%%
            %% w/ Focusing %%
            %%%%%%%%%%%%%%%%%
            
            tx_distance = sqrt((coord_sensor(curr_rx)-ix_cord(curr_z,curr_x)).^2 + (focus_index_z).^2); % Calculate actual distance
            tx_distance_ref = focus_index_z; % Subtract off reference distance (alignment)
            delay = (tx_distance - tx_distance_ref)/speed; % Calculate continuous time sample
            tau = Fs*(delay+(2*iz_cord(curr_z,1)/speed)); % Calculate the decimal index

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Backpropagation / Interpolation %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            index = floor(tau); % Interpolate correctly based
            tau_r = tau-index;
            if((300 <= index) && (index+1 < size_data))
                v = (1-tau_r)*raw_data(curr_rx, index) + tau_r*raw_data(curr_rx, index+1);
                img_recons(curr_z,curr_x) = img_recons(curr_z,curr_x) + (1)*v*(curr_weight(count)); % Store beamformed data into index
            end

            count = count + 1;
        end
    end
    % Envelope Detect the Beamformed Column
    img_recons_final(:,curr_x) = abs(hilbert(img_recons(:,curr_x)));
end
