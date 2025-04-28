function [img_recons_final, backprojections] = dericks_DAS_beamforming(Nz, Nx, N_sensor, ix_cord, iz_cord, coord_sensor,...
    dt, raw_data)  


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
intm_sensor_data = raw_data';
%initialize reconstruction space and backprojections
img_recons = zeros(Nz, Nx);
backprojections = zeros(Nz, Nx, N_sensor);

% start iterating through image space

for i = 1:Nz % each row, in this case determined by "Nz"
    for j = 1:Nx % each column, in this case determined by "Nx"
        for k = 1:N_sensor % iterate through each sensor
            % determine the distance between 
            % (j,i) <= use meshgrid
            % and
            % (j_k_sensor, i_k_sensor) <= given by coord_sensor

            x_component_dist = ix_cord(i,j) - coord_sensor(1, k);

            z_component_dist = iz_cord(i,j) - coord_sensor(2,k);

            distance = sqrt(x_component_dist^2 + z_component_dist^2);

            index = round((distance/ 1500) * (1/dt));

            if(index > 0)
                img_recons(i,j) = img_recons(i,j) + intm_sensor_data(index,k);
            
                backprojections(i,j,k) = backprojections(i,j,k) + img_recons(i,j);
            end
            
        end
    end
end


img_recons_final = (img_recons / N_sensor);
