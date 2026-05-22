
N_sensor = 32;

pert = 0;

x_vec = (-Nx/2:((Nx/2)-1))*dx;
z_vec = (0:(Nz-1))*dz;

lat_x = linspace(-Nx*dx/2, Nx*dx/2, Nx)*1e3;
depth_z = linspace(0, Nz*dz, Nz)*1e3;
    

raw_data = zeros(1,size(plane_wave_data,1),size(plane_wave_data,2));

raw_data(1,:,:) = plane_wave_data;

[ix_cord, iz_cord] = meshgrid(x_vec, z_vec);

coord_sensor = zeros(2,128);
coord_sensor(1,:) = ix_cord(1,:);

pert = 0.05;



[beamformed_images_nsi, beamformed_images_rect] = dof_analysis_v1_7(size(full_data_pt1,1), Nz, Nx, N_sensor, ix_cord, iz_cord, coord_sensor,...
    1/kgrid.dt, full_data_pt1, 5e6, size(full_data_pt1,3), 32, c0,pert);


number_samples = size(full_data_pt1,1);
full_images_nsi = beamformed_images_nsi;
% for j = 1:size(pert_list,2)

    total_image = zeros(Nz,Nx);
    curr_sample = zeros(1,Nz);
    curr_line = zeros(number_samples, Nz);

    for curr_x = 1:Nx
        for i = 1:number_samples
            curr_sample = ((total_images_nsi(j,i,:,curr_x)));
            curr_line(i,:)= transpose(squeeze(curr_sample));
        end
        compressed_line = max(curr_line, [], 1);   % 1 × numDepth
        total_image(:,curr_x) = compressed_line;
    end
    full_images_nsi(1,:,:) = total_image;
% end