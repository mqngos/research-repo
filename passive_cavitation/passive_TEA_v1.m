function  [rect_image, nsi_image] = passive_TEA_v1(raw_data, sensor_location, spatial_struct, pert) 

% ================= %
% Unpackages struct %
% ================= %

Nz = spatial_struct.Nz;
Nx = spatial_struct.Nx;
Fs = spatial_struct.Fs;
ix_cord = spatial_struct.ix_cord;
iz_cord = spatial_struct.iz_cord;
c0 = spatial_struct.c0;

organized_raw_data = raw_data;

intensity_image = zeros(Nz, Nx);
intensity_image_zm = zeros(Nz, Nx);
intensity_image_dc1 = zeros(Nz, Nx);
intensity_image_dc2 = zeros(Nz, Nx);
% all_zm_images = zeros(size(pert_list,2),Nz,Nx);
% all_dc1_images = zeros(size(pert_list,2),Nz,Nx);
% all_dc2_images = zeros(size(pert_list,2),Nz,Nx);
% all_nsi_images = zeros(size(pert_list,2),Nz,Nx);
% organized_raw_data = yq;
% curr_plot_grouped = zeros(466, size(organized_raw_data,1),size(organized_raw_data,2)); % Used to debug
curr_plot = zeros(128,size(organized_raw_data,2));
summed_sensors_rect = zeros(32,size(organized_raw_data(:,50:end-200),2));
summed_sensors_zm = zeros(32,size(organized_raw_data(:,50:end-200),2));
summed_sensors_dc1 = zeros(32,size(organized_raw_data(:,50:end-200),2));
summed_sensors_dc2 = zeros(32,size(organized_raw_data(:,50:end-200),2));


    start_window_index = 1;
    end_window_index = size(raw_data,2);

    dc_offset = pert;


for curr_x = 232:282 % {-2.5 [mm], 2.5 [mm]} , each scanline being 0.1 [mm] wide
% for curr_x = 256:257
    for curr_z = 1:Nz    
        % active_elements = 49:80; % Changed based on sub-aperture size, i.e. for 32-elements <=> 49:80
        active_elements = (curr_x-207):(curr_x-207)+31; % Declare elements for scanline

        if (size(active_elements,2) > 0)
            half_elements = floor(size(active_elements,2)/2);
            weights_zm = ones(1,size(active_elements,2));
            weights_zm(17:end) = -1;
            weights_dc1 = ones(1,size(active_elements,2));
            weights_dc1(17:end) = -1;
            weights_dc1 = weights_dc1 + dc_offset;
            weights_dc2 = -1*ones(1,size(active_elements,2));
            weights_dc2(17:end) = 1;
            weights_dc2 = weights_dc2 + dc_offset;
        end

        distances = zeros(1,size(active_elements,2));
        for i = 1:size(active_elements,2)

            dist = sqrt(iz_cord(curr_z,1).^2+((sensor_location(active_elements(i))-ix_cord(1,curr_x)).^2));

            distances(i) = dist;
        end

        d_r = min(distances);

        distances = distances - d_r;

        delays = (Fs.*distances)./(c0);
        
        delays_index = floor(delays);

        delays_r = delays-delays_index;

        curr_intensity = 0;
        curr_intensity_zm = 0;
        curr_intensity_dc1 = 0;
        curr_intensity_dc2 = 0;

        if (size(active_elements,2) > 0)
            for i = 1:size(active_elements,2)

                curr_sensor = active_elements(i);
                curr_time_series = zeros(1,size(organized_raw_data,2));

                for j = 1:size(organized_raw_data,2)

                    if(j+delays_index(i) > 0 && j+delays_index(i)+1<=size(organized_raw_data,2))
                        curr_plot(curr_sensor,j) = (1-delays_r(i)).*organized_raw_data(curr_sensor,j+delays_index(i)) + (delays_r(i)).*organized_raw_data(curr_sensor,1+j+delays_index(i));
                    end

                end

            end

            % start_window_index = 50; % Minimum needed, excludes burst at the start
            % end_window_index = size(curr_plot,2)-200;

            % start_window_index = 200; % Used for 10 [mm] approx
            % end_window_index = 300;

            % start_window_index = 460; % Used for 20 [mm] approx
            % end_window_index = 560;

            % start_window_index = 760; % Used for 30 [mm] approx
            % end_window_index = 860;

            % 

            

            
            for i = 1:size(active_elements,2)
                % summed_sensors_zm(i,:) = weights_zm(i).*curr_plot(active_elements(i),50:end-200);
                % summed_sensors_dc1(i,:) = weights_dc1(i).*curr_plot(active_elements(i),50:end-200);
                % summed_sensors_dc2(i,:) = weights_dc2(i).*curr_plot(active_elements(i),50:end-200);
                summed_sensors_rect(i,start_window_index:end_window_index) = curr_plot(active_elements(i),start_window_index:end_window_index);
                summed_sensors_zm(i,start_window_index:end_window_index) = weights_zm(i).*curr_plot(active_elements(i),start_window_index:end_window_index);
                summed_sensors_dc1(i,start_window_index:end_window_index) = weights_dc1(i).*curr_plot(active_elements(i),start_window_index:end_window_index);
                summed_sensors_dc2(i,start_window_index:end_window_index) = weights_dc2(i).*curr_plot(active_elements(i),start_window_index:end_window_index);
            end

            
            b_rect = abs(hilbert(sum(summed_sensors_rect,1)));
            b_dc1 = abs(hilbert(sum(summed_sensors_dc1,1)));
            b_dc2 = abs(hilbert(sum(summed_sensors_dc2,1)));
            b_zm = abs(hilbert(sum(summed_sensors_zm,1)));

            b_nsi = ((( b_dc1 + b_dc2 ) ./2 ) - b_zm);

            curr_inteinsity_nsi = sum((b_nsi).^2);
            
            curr_intensity = sum(b_rect.^2);
            curr_intensity_zm = sum((b_zm).^2);
            curr_intensity_dc1 = sum((b_dc1).^2);
            curr_intensity_dc2 = sum((b_dc2).^2);

        end
        % 
        intensity_image(curr_z, curr_x) = curr_intensity;
        intensity_image_zm(curr_z, curr_x) = curr_inteinsity_nsi;
        % intensity_image_dc1(curr_z, curr_x) = curr_intensity_dc1;
        % intensity_image_dc2(curr_z, curr_x) = curr_intensity_dc2;
    end
end

nsi_image = intensity_image_zm;
rect_image = intensity_image;