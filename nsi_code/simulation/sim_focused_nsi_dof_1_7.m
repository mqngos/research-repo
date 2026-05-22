    clear all;
    close all;
    
        
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
        
        % create the computational grid
        % PML_size = 20;          % size of the PML in grid points
        dx_wavelengths = 0.6306;            % grid point spacing in the x direction [wavelengths]
        dx = 0.3048e-3; %[mm]
        % dz_wavelengths = 0.5;            % grid point spacing in the y direction [wavelengths]
        % dz = dz_wavelengths*(c0/source_freq); %[mm]
        dz = 0.1e-3;
        Nx = 128;  % number of sensors
        Nz = 700;  % number of grid points in the y (column) direction
    
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
    
        % medium.density(128, Nx/2) = 100*rho0;
        % medium.sound_speed(128, Ny/2) = c0*2;
    
        
        % create initial pressure distribution using makeDisc
        
        
        % smooth the initial pressure distribution and restore the magnitude
        
        % define a binary line sensor (along Nx = 1)
        % 32 sensors long
    
        sensor.mask = zeros(Nz, Nx);
    
        sensor.mask(1, :) = 1;
    
        % source_pos = Nx/2;
        source.p_mask = zeros(Nz, Nx);
        source.p_mask(1, :) = 1;
        
        
        % create the time array
        kgrid.makeTime(medium.sound_speed);
    
        % kgrid.t_array = (0:1500).*(1/Fs);
        % Fs = 500e6;
        % kgrid.dt = 1/Fs;
        kgrid.t_array = (0:5500).*kgrid.dt;
    % 
        % kgrid.dt = 2e-8;
        
        T = 1 / source_freq;    % Period [s]
        
        % input_args = {'PMLInside', false, 'PMLSize', [80, 80], 'PMLAlpha', 1.5, ...
        %           'PlotPML', false, 'Smooth', false, 'PlotSim', true};
        input_args = {'PMLInside', false, 'PMLSize', 100,'PMLAlpha', pml_alpha 'Smooth', true, 'PlotPML', true};
    
        % %% run the simulation/DAQ #1 -- Zero-Mean Input
    
        % points_per_wavelength = 250;
        % 
        % % t = (0:points_per_wavelength-1).*();
        % 
        % % betteroni = smooth(10*sin(2 * pi * source_freq * t));
        % 
        % 
        % first_minimum_index_x = find(abs(betteroni) < 0.01);
        % 
        % cut_off = first_minimum_index_x(3);
        % 
        % 
        % source.p = zeros(size(kgrid.t_array));
        % 
        % source.p(1:cut_off) = betteroni(1:cut_off);
    
        Fs = 100e6;
    
        source.p = 10.*toneBurst(Fs, 5e6, 2);
    
    
        % %% Place Wire Targets
        % Properties of Steel
        rho_steel = 7700;
        c_steel = 5050;
    
        aperture_size = 32;
    
    
        
        
    
        source.p_mask = zeros(Nz, Nx);
        sensor.mask = zeros(Nz,Nx);
        sensor.mask(1, :) = 1;
        % source.p_mask(1,:) = 1;
    
    
        medium.density = rho0*ones(Nz,Nx);
        medium.sound_speed = c_steel; % 10 mm
        medium.sound_speed(66, ((Nx/2))) = c_steel; % 10 mm
        medium.density(66, ((Nx/2))) = rho_steel;
        % 
    
        medium.density = rho0*ones(Nz,Nx);
        medium.sound_speed = c0*ones(Nz,Nx);
        % medium.sound_speed(66, ((Nx/2))) = c_steel; % 10 mm
        medium.density(131, ((Nx/2))) = rho_steel;
    
        [ix_cord, iz_cord] =meshgrid(kgrid.y_vec, kgrid.x_vec);
    
        % lat_x = ix_cord(1,:)*1e3;
        % depth_z = iz_cord(:,1)*1e3-min(iz_cord(:,1))*1e3;
    
        iz_cord = kgrid.x(:,1) - min(kgrid.x(:,1));
    
        % %% Testing Environment for Focused Sensors
    
        % %% Actual Simualtion Stuff
    
        % ============================================
        % Stupid GIF stuff
        % ============================================
    
        [gifData, gifMap] = imread('spinning_rei.gif', 'Frames', 'all');
        numFrames = size(gifData,4);
    
        % x, y, width, height VVV
    
        hFig = figure('MenuBar','none','ToolBar','none','NumberTitle','off', ...
                      'Name','Processing','Position',[600 400 500 800]);
    
        % ---------------- GIF area ----------------
        hAxGif = axes('Parent', hFig, 'Position',[0.15 0.55 0.7 0.45]);
        hImg = imshow(gifData(:,:,:,1), gifMap, 'Parent', hAxGif);
    
        hAxMedium = axes('Parent', hFig, 'Position',[0.35 0.15 0.25 0.5]);
        hImg_hm = imagesc(medium.density(100:170,:),'Parent', hAxMedium);
    
        set(hAxMedium, 'XColor', 'None');   % show x-axis ticks & labels
        set(hAxMedium, 'YColor'); 
    
        % ---------------- Progress bar area ----------------
        hAxBar = axes('Parent', hFig, 'Units','normalized', ...
                      'Position',[0.1 0.15 0.8 0.05]);
    
        set(hAxBar, 'Visible','on');  % Hide tick  s, axes, everything
    
        % Create bar frame (background)
        patch('Parent', hAxBar, ...
              'XData',[0 1 1 0], 'YData',[0 0 1 1], ...
              'FaceColor',[0.85 0.85 0.85], ...
              'EdgeColor','none');
    
        % Create bar fill
        hBar = patch('Parent', hAxBar, ...
                     'XData',[0 0 0 0], 'YData',[0 0 1 1], ...
                     'FaceColor',[0.2 0.6 1], ...
                     'EdgeColor','none');
    
        % % ---------------- Status text ----------------
        hText = uicontrol('Style','text','Parent',hFig, ...
                          'Units','normalized', ...
                          'Position',[0.1 0.05 0.8 0.1], ...
                          'FontSize',11, ...
                          'HorizontalAlignment','center', ...
                          'String','Initializing...');
    
    
    
    
        tic;
        fprintf('Starting 10 mm Target!\n');
    
        %
        % To move 2 [mm], its around 13 pixels
        frame = 1;
    
        % synth_plane_wave_data_20_mm_dof = zeros(size(lateral_list,2),size(depth_list,2), 97, 128, size(kgrid.t_array,2));
    
    
        % full_data_pt3 = zeros(3,size(depth_list,2), 97, 128, size(kgrid.t_array,2));
        % 
        % full_data_pt4 = zeros(4,size(depth_list,2), 97, 128, size(kgrid.t_array,2));
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %   DEFINE TIME DELAYS   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%
    
        num_active_sensors = 32;
        na = Nx - num_active_sensors + 1;
        % focus_point_x = 0;
        % focus_point_z = depth_z(130);
        % elementPosMm = lat_x;
        curr_scan = 1;
        f0 = source_freq;
        
        x_t = 10.*toneBurst(Fs, 5e6, 1); 
    
        source.p = x_t;
      
        % focus_point_x = 0;
        % delayed_wvfms = createDelayedWvfms_kWave(num_active_sensors, focus_point_x, focus_point_z, elementPosMm, c0, curr_scan, f0,x_t, kgrid.dt);
    
    
    
    % 20 [mm] Foocus Simulaton
    % depth_list = [66,73,79,86,92,99,105,112, ...
    %                       118,125,131,137,144,150,157, ...
    %                       163,170,176,183,189,196];
    % lateral_list = 61:69;
    
    % 20 [mm] Foocus Simulaton (Upsampled)
    % depth_list = [101,102,105,108,111, 115, 118, 121, 125, 128, 131 ...
    %               134, 137, 141, 144, 147, 150, 154, 157, 160, 163 ...
    %               167, 170, 173, 176, 180, 183, 186, 189, 193 ,196];

    % depth_list = 101:5:302;
    % depth_list = [10,20,30];
    % depth_list = 201:501; % 20 [mm] focus
    depth_list = 1:1; % 30 [mm] focus
    % depth_list = 501;
    % depth_list = [101,201,301];
    % depth_list = 
    lateral_list = 1:1;
    
    % 10 [mm] Focus Simulation
    % depth_list = [34, 37, 40, 43, 47, 50, 53, 56, 60,63, 66,69, 72, 76,79,82,86,89,92,95,98]; 
    % lateral_list = 61:69;
    full_data_pt1 = zeros(size(depth_list,2), 128, size(kgrid.t_array,2));
    
        
    % new_depth_list = [131];
    
    for k = 5:5
    % for k = 5:5
        for j = 1:size(depth_list,2)
        % for j = size(depth_list,2):size(depth_list,2)
        % for j = 201:201
        % for j = 1:size()
            
            medium.density = rho0*ones(Nz,Nx);
            medium.sound_speed = c0*ones(Nz,Nx);
    
            % medium.density(depth_list(j), lateral_list(k)) = rho_steel;
            medium.density(201, 65) = rho_steel;
            % medium.density(102, lateral_list(k)) = rho_steel;
    
            source.p_mask(1,:) = 0;
            sensor.mask(1,:) = 0;
            sensor.mask(1,:) = 1;
            source.p_mask(1, 49:49+num_active_sensors-1) = 1; % Assign Delayed Waveforms
    
            focus_point_z = 30; % mm
            half_sensors = num_active_sensors/2;
            % focus_point_x = -dx*1e3/2 ;
            focus_point_x = 0;
    
                
            num_active_sensors = 32;
                % 
                % full_time_delays = zeros(na, size(lat_x, 2)); % Size of [Num. Acquistions x Num. Elements]
                % 
                % % d_angle = 1; % In Degrees
                % % start_angle = -16; % In Degrees
                % % end_angle = 16; % In Degrees
                % % 
                % % angles = start_angle:d_angle:end_angle; % in degrees
                % % 
                % % angles = angles * (pi/180); % Degrees to Radians
                % 
                % % tau_list = (sin(angles) * dx)/c0;
                % 
                % tx_distance = sqrt((focus_point_x*1e-3 - ix_cord(1,49:49+num_active_sensors-1)).^2 + (focus_point_z*1e-3).^2);
                % 
                % time_delays = tx_distance / c0;
                % 
                % time_delays = time_delays - min(time_delays);
    
                % ref_distance = sqrt((ix_cord(1,64) - focus_point_x*1e-3).^2 + (focus_point_z*1e-3).^2);
                % 
                % tx_distances = sqrt((ix_cord(1,49:49+num_active_sensors-1) - focus_point_x*1e-3).^2 + (focus_point_z*1e-3).^2);
                % 
                % tone_burst_offset = (tx_distances - min(tx_distances) )/ (c0 * Fs);
        
            x_f = focus_point_x*1e-3;
    
            z_f = focus_point_z*1e-3;
    
    
    
    
            R = sqrt((ix_cord(1,49:49+num_active_sensors-1) - x_f).^2 + z_f.^2);
            tau_TX = (max(R) - R) / (c0 * kgrid.dt);
    
    
                % tone_burst_offset = max(tone_burst_offset) - tone_burst_offset ;
    
            x_t = 10.*toneBurst(Fs, 5e6, 1);
                
                
    
            % x_t = 10.*toneBurst(1/(kgrid.dt), 5e6, 1, 'SignalOffset', tau_TX);
    
    
            
    
            source.p = x_t;
    
                % tone_burst_offset = dx * (1:num_active_sensors) * sin(start_angle * pi/180) / (c0 * kgrid.dt);
                % 
                % tone_burst_offset = tone_burst_offset - min(tone_burst_offset);
    
                % [ test] = createDelayedWvfms_kWave(num_active_sensors, focus_point_x, focus_point_z, ix_cord(1,:), c0, 48, f0,x_t, kgrid.dt);
                
                % source.p = delayed_wvfms;
                % 
            plane_wave_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
            
            % if(k > 0 && k <= 3)
            full_data_pt1(j,:,:) = plane_wave_data;
                % end
                % if(k > 3 && k <= 6)
                    % full_data_pt2(rj,i,:,:) = plane_wave_data;
                % end
    
                % synth_plane_wave_data_20_mm_dof(k,j,i,:,:) = plane_wave_data;
    
                % Code for fun lol
                % 
                % plot(squeeze(full_data_pt1(3,64,:)))
                % 
            frame = frame + 1;
    
            if frame > numFrames
                frame = 1;
            end
            progress = i/na;
            set(hImg, 'CData', gifData(:,:,:,frame));
    
            y_min = max(depth_list(j)-10, 1);
            y_max = min(depth_list(j)+10, Nz);
    
            set(hImg_hm, 'CData', medium.density(y_min:y_max,:));
    
            set(hAxMedium, 'YLim', [1 (y_max - y_min + 1)]);
    
            yticks = get(hAxMedium, 'YTick');
    
            set(hAxMedium, 'YTickLabel', (y_min + yticks - 1)*dz*1e3);
    
            ylabel(hAxMedium,'Depth (mm)');
    
            set(hBar, 'XData', [0 progress progress 0]);
    
            set(hText, 'String', sprintf('Lateral Sim: %d / %d Depth Sim: %d / %d', k,size(lateral_list,2),j,size(depth_list,2)));
            drawnow;
    
    
        end
    
    
    end
    elapsed = toc;
    fprintf('Finished Targets!\n Finished in %d\n', elapsed);

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
