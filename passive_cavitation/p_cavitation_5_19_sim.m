clear all;
close all;
    

%% Simulation Task:

% We want to simulate the narrowband behavior of some point-source emitted
% at some cavitation time, t_c. 

% From there, we will try to beamform with TEA/NSI to see how this system
% behaves.
        
% =========================================================================
% SIMULATION
% =========================================================================

p0 = 150e6;                    % source pressure [Pa]
c0 = 1540;                     % sound speed [m/s]
rho0 = 998;                    % density [kg/m^3]
alpha_0 = 0.25;                % absorption coefficient [dB/(MHz^2 cm)]
sigma = 2;                     % shock parameter
source_freq = 5e6;             % frequency [Hz]
points_per_wavelength = 100;   % number of grid points per wavelength at f0
wavelength_separation = 15;    % separation between the source and detector
pml_size = 80;                 % PML size
pml_alpha = 1.5;               % PML absorption coefficient [Np/grid point]
CFL = 0.25;                    % CFL number



waveToMm = 0.3048/1.0510;

dx = 0.1*1e-3; % [m] % May need to change because of k-Wave run time...
% dz_wl = 0.5; % grid point sacing in z direction [wavelengths]
dz = 0.1*1e-3; % [m]


Fs = 1.923076923076923e+07; % From Receive struct from L9-4, decimSampleRate i believe!

Nx = 512; % Number pixels laterally (left-to-right)
Nz = 416; % Number pixels axially (depth-wise)

x_vec = (-Nx/2:((Nx/2)-1))*dx; % [m]
z_vec = (0:(Nz-1))*dz; % Creates in wavelengths

[ix_cord, iz_cord] = meshgrid(x_vec, z_vec); % Creation of meshgrid

kgrid = kWaveGrid(Nz, dz, Nx, dx);

% Fs = 100e6;


% define the properties of the propagation medium
x_px = wavelength_separation * points_per_wavelength;
x = x_px * dx;
mach_num = p0 / (rho0 * c0.^2);
k = 2 * pi * source_freq / c0;
BonA = 2 * (sigma / (mach_num * k * x) - 1);
medium.BonA = BonA;
medium.density = rho0*ones(Nz,Nx);
medium.sound_speed = c0*ones(Nz,Nx);


% ** Define transducer at the top of simulation space
sensor.mask = zeros(Nz, Nx);
sensor.mask(1, :) = 1;

% ** Define scatter at a certain position (single scatterer)

source.p_mask = zeros(Nz, Nx);
source.p_mask(Nx/2, Nz/2) = 1; % In this case, in the middle of our 

% create the time array
kgrid.makeTime(medium.sound_speed);

kgrid.t_array = (0:1500).*(1/Fs);

kgrid.t_array = (0:1500).*kgrid.dt;

input_args = {'PMLInside', false, 'PMLSize', 100,'PMLAlpha', pml_alpha 'Smooth', true, 'PlotPML', true};

%%
% =========================================================================
% SIMULATION, i.e., literally just backscattering LMFAO
% =========================================================================

% Here, we will model x_t as a Broadband signal (in this case, modeled as a
% very narrow gaussian source)

% We will backpropagate this signal onto our elements, creating a
% time-series stacked of the dimensions [128 x 1501], w/ Fs = 1.9231e+07 to
% match that of an L9-4 Transducer

% Laterally Locations of Transducer Elements:

sensor_location = (1e-3)* [ ...
-19.3548 -19.0500 -18.7452 -18.4404 -18.1356 -17.8308 -17.5260 -17.2212 ...
-16.9164 -16.6116 -16.3068 -16.0020 -15.6972 -15.3924 -15.0876 -14.7828 ...
-14.4780 -14.1732 -13.8684 -13.5636 -13.2588 -12.9540 -12.6492 -12.3444 ...
-12.0396 -11.7348 -11.4300 -11.1252 -10.8204 -10.5156 -10.2108 -9.9060 ...
-9.6012 -9.2964 -8.9916 -8.6868 -8.3820 -8.0772 -7.7724 -7.4676 ...
-7.1628 -6.8580 -6.5532 -6.2484 -5.9436 -5.6388 -5.3340 -5.0292 ...
-4.7244 -4.4196 -4.1148 -3.8100 -3.5052 -3.2004 -2.8956 -2.5908 ...
-2.2860 -1.9812 -1.6764 -1.3716 -1.0668 -0.7620 -0.4572 -0.1524 ...
0.1524 0.4572 0.7620 1.0668 1.3716 1.6764 1.9812 2.2860 ...
2.5908 2.8956 3.2004 3.5052 3.8100 4.1148 4.4196 4.7244 ...
5.0292 5.3340 5.6388 5.9436 6.2484 6.5532 6.8580 7.1628 ...
7.4676 7.7724 8.0772 8.3820 8.6868 8.9916 9.2964 9.6012 ...
9.9060 10.2108 10.5156 10.8204 11.1252 11.4300 11.7348 12.0396 ...
12.3444 12.6492 12.9540 13.2588 13.5636 13.8684 14.1732 14.4780 ...
14.7828 15.0876 15.3924 15.6972 16.0020 16.3068 16.6116 16.9164 ...
17.2212 17.5260 17.8308 18.1356 18.4404 18.7452 19.0500 19.3548 ];

plane_wave_data = zeros(128,1501);

z_s = 20e-3; % some random axial depth of scatterer
x_s = 0; % some random lateral depth of scatterer

for curr_sensor = 1:128 % Iterate through all sensors

    distance = sqrt((sensor_location(curr_sensor) - 0).^2 + z_s.^2);

    t_c0 = distance/c0;

    plane_wave_data(curr_sensor,:) = 5000*(1- (((t-t_c0).^2)./(sigma.^2))) .* exp(-(((t-t_c0).^2)./(2*sigma.^2)));    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beamforming Portion via TEA %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ========================= %
% Define the physical space %
% ========================= %

waveToMm = 0.3048/1.0510;
c0 = 1540; % Change depending on medium... water = 1540, phantom = 1490...?

dx = 0.1*1e-3; % [mm]
dz = 0.5; % grid point sacing in z direction [wavelengths]

Fs = 1.923076923076923e+07; % From Receive struct from L9-4, decimSampleRate i believe!

Nx = 512; % Number pixels laterally (left-to-right)
Nz = 416; % Number pixels axially (depth-wise)

x_vec = (-Nx/2:((Nx/2)-1))*dx; % [m]
z_vec = (0:(Nz-1))*dz; % Creates in wavelengths

z_vec = z_vec*waveToMm*1e-3; % Convert to [m]

[ix_cord, iz_cord] = meshgrid(x_vec, z_vec); % Creation of meshgrid

spatial_struct.Nx = Nx;
spatial_struct.Nz = Nz;
spatial_struct.Fs = 4*Fs;
spatial_struct.ix_cord = ix_cord;
spatial_struct.iz_cord = iz_cord;
spatial_struct.c0 = c0;


% ======================= %
%   Upsampling / RF Data  %
% ======================= %

U = 4;
x  = 1:size(plane_wave_data,2);
y  = plane_wave_data;
xq = 1:1/U:size(plane_wave_data,2);
yq = spline(x,y,xq);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call Beamforming Logic / RF Data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% pert_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
% pert_list = [0.05,0.07,0.09, 0.2,0.4,0.6];

% pert_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];


% function signature:
% [rect_image, nsi_image] = passive_TEA_v1(raw_data, sensor_location, pert) 

[rect_image, nsi_image] = passive_TEA_v1(yq,sensor_location,spatial_struct,0.5);


%%
% Simulated RF data is stored in the [128 x 1501] matrix, plane_wave_data!

% 
% 
% % ============================================
% % Stupid GIF stuff
% % ============================================
% 
% [gifData, gifMap] = imread('spinning_rei.gif', 'Frames', 'all');
% numFrames = size(gifData,4);
% 
% % x, y, width, height VVV
% 
% hFig = figure('MenuBar','none','ToolBar','none','NumberTitle','off', ...
%               'Name','Processing','Position',[600 400 500 800]);
% 
% % ---------------- GIF area ----------------
% hAxGif = axes('Parent', hFig, 'Position',[0.15 0.55 0.7 0.45]);
% hImg = imshow(gifData(:,:,:,1), gifMap, 'Parent', hAxGif);
% 
% hAxMedium = axes('Parent', hFig, 'Position',[0.35 0.15 0.25 0.5]);
% hImg_hm = imagesc(medium.density(100:170,:),'Parent', hAxMedium);
% 
% set(hAxMedium, 'XColor', 'None');   % show x-axis ticks & labels
% set(hAxMedium, 'YColor'); 
% 
% % ---------------- Progress bar area ----------------
% hAxBar = axes('Parent', hFig, 'Units','normalized', ...
%               'Position',[0.1 0.15 0.8 0.05]);
% 
% set(hAxBar, 'Visible','on');  % Hide tick  s, axes, everything
% 
% % Create bar frame (background)
% patch('Parent', hAxBar, ...
%       'XData',[0 1 1 0], 'YData',[0 0 1 1], ...
%       'FaceColor',[0.85 0.85 0.85], ...
%       'EdgeColor','none');
% 
% % Create bar fill
% hBar = patch('Parent', hAxBar, ...
%              'XData',[0 0 0 0], 'YData',[0 0 1 1], ...
%              'FaceColor',[0.2 0.6 1], ...
%              'EdgeColor','none');
% 
% % % ---------------- Status text ----------------
% hText = uicontrol('Style','text','Parent',hFig, ...
%                   'Units','normalized', ...
%                   'Position',[0.1 0.05 0.8 0.1], ...
%                   'FontSize',11, ...
%                   'HorizontalAlignment','center', ...
%                   'String','Initializing...');
% 
% 
% 
% 
% tic;
% fprintf('Starting 10 mm Target!\n');
% 
% frame = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Run simulation environment   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kspaceFirstOrder2D <= not CUDA based, w/ visualization

% kspaceFirstOrder2DG <= for CUDA based, w/ no visualization

plane_wave_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});

% 
% num_active_sensors = 32;
% na = Nx - num_active_sensors + 1;
% % focus_point_x = 0;
% % focus_point_z = depth_z(130);
% % elementPosMm = lat_x;
% curr_scan = 1;
% f0 = source_freq;
% 
% x_t = 10.*toneBurst(Fs, 5e6, 1); 
% 
% source.p = x_t;
% 
% 
% depth_list = 1:1; % 30 [mm] focus
% lateral_list = 1:1;
% 
% 
% 
% 
% full_data_pt1 = zeros(size(depth_list,2), 128, size(kgrid.t_array,2));

% for k = 5:5
%     for j = 1:size(depth_list,2)
% 
% 
%             % end
%             % if(k > 3 && k <= 6)
%                 % full_data_pt2(rj,i,:,:) = plane_wave_data;
%             % end
% 
%             % synth_plane_wave_data_20_mm_dof(k,j,i,:,:) = plane_wave_data;
% 
%             % Code for fun lol
%             % 
%             % plot(squeeze(full_data_pt1(3,64,:)))
%             % 
%         frame = frame + 1;
% 
%         if frame > numFrames
%             frame = 1;
%         end
%         progress = i/na;
%         set(hImg, 'CData', gifData(:,:,:,frame));
% 
%         y_min = max(depth_list(j)-10, 1);
%         y_max = min(depth_list(j)+10, Nz);
% 
%         set(hImg_hm, 'CData', medium.density(y_min:y_max,:));
% 
%         set(hAxMedium, 'YLim', [1 (y_max - y_min + 1)]);
% 
%         yticks = get(hAxMedium, 'YTick');
% 
%         set(hAxMedium, 'YTickLabel', (y_min + yticks - 1)*dz*1e3);
% 
%         ylabel(hAxMedium,'Depth (mm)');
% 
%         set(hBar, 'XData', [0 progress progress 0]);
% 
%         set(hText, 'String', sprintf('Lateral Sim: %d / %d Depth Sim: %d / %d', k,size(lateral_list,2),j,size(depth_list,2)));
%         drawnow;
% 
% 
%     end
% 
% 
% end
% elapsed = toc;
% fprintf('Finished Targets!\n Finished in %d\n', elapsed);

%% Store data:
folder = fullfile('simulation_data', 'intensity_map/','data_1_14/', 'post_sim_data');

filename = fullfile(folder, 'plane_wave_data.mat');

save(filename, 'full_data_pt1','-v7.3');


    % bw_data = zeros(size(synth_plane_wave_data_20_mm_dof,1), ...
    %     size(synth_plane_wave_data_20_mm_dof,3), ...
    %     size(synth_plane_wave_data_20_mm_dof,4), ...
    %     size(synth_plane_wave_data_20_mm_dof,5));
    % 
    % bw_data = (synth_plane_wave_data_20_mm_dof);
    % 
    % 
    % 
    % 
    % % dual_line_data = synth_plane_wave_data_20_mm_dof;
    % % 
    % % dual_line_data_left_pt1 = synth_plane_wave_data_20_mm_dof(1,1:5,:,:);
    % % dual_line_data_left_pt2 = synth_plane_wave_data_20_mm_dof(1,6:10,:,:);
    % % dual_line_data_left_pt3 = synth_plane_wave_data_20_mm_dof(1,11:15,:,:);
    % % dual_line_data_left_pt4 = synth_plane_wave_data_20_mm_dof(1,16:20,:,:);
    % % 
    % % dual_line_data_right_pt1 = synth_plane_wave_data_20_mm_dof(2,1:5,:,:);
    % % dual_line_data_right_pt2 = synth_plane_wave_data_20_mm_dof(2,6:10,:,:);
    % % dual_line_data_right_pt3 = synth_plane_wave_data_20_mm_dof(2,11:15,:,:);
    % % dual_line_data_right_pt4 = synth_plane_wave_data_20_mm_dof(2,16:20,:,:);
    % 
    % folder = fullfile('simulation_data', 'generic_data/');
    % % % % 
    % filename = fullfile(folder, 'timer.mat');
    % % % % 
    % save(filename, 'Nx', '-v7.3');
    % save(filename, 'dual_line_data_left_pt1','dual_line_data_left_pt2','dual_line_data_left_pt3','dual_line_data_left_pt4');

    % 
    % i = 1;
    % plot(squeeze(synth_plane_wave_data_20_mm_dof(3,i,64,:)))
    % i = i + 1;

    % save('simulation data\10_mm_to_24_mm_step_2_mm_focus_30_mm', 'help',dx);
    % 

    % 
    % filename = fullfile(folder, 'raw_data_pt2.mat');
    % 
    % raw_data_2 = synth_plane_wave_data_20_mm_dof(10:18,:,:,:);
    % 
    % save(filename, 'raw_data_2');

%% Working Simulation Code
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%% TEST_SIMUALTIONS %%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % 
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % %   DEFINE TIME DELAYS   %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % num_active_sensors = 128;
    % focus_point_x = 0;
    % focus_point_z = depth_z(130);
    % elementPosMm = lat_x;
    % curr_scan = 1;
    % f0 = source_freq;
    % 
    % x_t = 10.*toneBurst(Fs, 5e6, 2); 
    % 
    % delayed_wvfms = createDelayedWvfms_kWave(num_active_sensors, focus_point_x, focus_point_z, elementPosMm, c0, curr_scan, f0,x_t, kgrid.dt);
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % %   ASSIGN TIME DELAYS   %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % source.p = delayed_wvfms;
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % %  CREATE WIRE TARGETS  %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % medium.density = rho0*ones(Nz,Nx);
    % medium.sound_speed = c0*ones(Nz,Nx);
    % medium.density(66, ((Nx/2))) = rho_steel; % 10 mm
    % medium.density(130, ((Nx/2))) = rho_steel; % 20 mm
    % medium.density(195, ((Nx/2))) = rho_steel; % 30 mm
    % 
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%  COMMIT SIMUATION  %%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % source.p_mask(1,:) = 0;
    % sensor.mask(1,:) = 0;
    % sensor.mask(1,:) = 1;
    % 
    % source.p_mask(1, :) = 1;
    % 
    % plane_wave_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
    % 
    % 
    % 
    % figure
    % plot(plane_wave_data(64,:))
