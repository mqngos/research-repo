    clearvars;
    
    % =========================================================================
    % SIMULATION
    % =========================================================================
    
    p0 = 10e6;                     % source pressure [Pa]
    c0 = 1500;                     % sound speed [m/s]
    rho0 = 1000;                   % density [kg/m^3]
    alpha_0 = 0.25;                % absorption coefficient [dB/(MHz^2 cm)]
    sigma = 2;                     % shock parameter
    source_freq = 1e6;             % frequency [Hz]
    points_per_wavelength = 100;   % number of grid points per wavelength at f0
    wavelength_separation = 15;    % separation between the source and detector
    pml_size = 80;                 % PML size
    pml_alpha = 1.5;               % PML absorption coefficient [Np/grid point]
    CFL = 0.25;                    % CFL number
    
    % create the computational grid
    % PML_size = 20;          % size of the PML in grid points
    Nx = 256;  % number of grid points in the x (row) direction
    Ny = 256;  % number of grid points in the y (column) direction
    dx = 0.05e-3;            % grid point spacing in the x direction [m]
    dy = 0.05e-3;            % grid point spacing in the y direction [m]
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    
    medium.sound_speed = c0;
    medium.density = rho0;
    medium.alpha_power = 2;
    medium.alpha_coeff = alpha_0;
    
    
    % define the properties of the propagation medium
    x_px = wavelength_separation * points_per_wavelength;
    x = x_px * dx;
    medium.sound_speed = 1500;           % [m/s]
    mach_num = p0 / (rho0 * c0.^2);
    k = 2 * pi * source_freq / c0;
    BonA = 2 * (sigma / (mach_num * k * x) - 1);
    medium.BonA = BonA;
    
    % create initial pressure distribution using makeDisc
    
    
    % smooth the initial pressure distribution and restore the magnitude
    
    % define a binary line sensor
    sensor.mask = zeros(Nx, Ny);
    sensor.mask(1, :) = 1;
    
    source_pos = Nx/2;
    source.p_mask = zeros(Nx, Ny);
    source.p_mask(source_pos, Ny/2) = 1;
    
    
    % create the time array
    kgrid.makeTime(medium.sound_speed);
    
    T = 1 / source_freq;    % Period [s]
    
    % % Create the source term: one period of sine, then zero

    source.p = zeros(size(kgrid.t_array));
    one_period_indices = kgrid.t_array <= T;    % Logical mask for one period
    
    source_time_signal = zeros(size(kgrid.t_array));
    source_time_signal(one_period_indices) = p0 * cos(2 * pi * source_freq * kgrid.t_array(one_period_indices));
    
    num_source_points = sum(source.p_mask(:)); % should be 1 if single point
    source.p = repmat(source_time_signal, num_source_points, 1);
    
    % % Disc Simulation -- Source
    % source.p0 = 3 * makeDisc(Nx, Ny, Nx/2, Ny/2, 4) + 3 * makeDisc(Nx, Ny, Nx/4, Ny/2, 4) + 3 * makeDisc(Nx, Ny, 3*Nx/4, Ny/2, 4);
    
    
    
    
    % set the input arguements: force the PML to be outside the computational
    % grid; switch off p0 smoothing within kspaceFirstOrder2D
    input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false, 'PlotSim', true};
    
 
    % run the simulation -- inversion 
    % CHANGE:
    % Here, instead of having a sine func (instead of inversion), we will
    % modulate a cosine function instead with:

    % Acos(wt)
    % Bcos(3wt)

    % source frequency => Acos(wt) * Bcos(3wt) = ...
    % => \frac{AB}{2} \cdot cos(4wt) * cos(2wt)
    % This will only give us 4th and 2nd harmonics

    omega = 2 * pi * source_freq;

    help = cos(omega * kgrid.t_array(one_period_indices));

    help_1 = cos(3*omega * kgrid.t_array(one_period_indices));

    help_2 = help .* help_1;
    
    source.p(one_period_indices) = p0 * help_2;
    
    sensor_data_modulated_2nd_4th = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

    % run the simulation -- inversion 
    % CHANGE:
    % Here, instead of having a sine func (instead of inversion), we will
    % modulate a cosine function instead with:

    % Acos(wt)
    % Bcos(2wt)

    % source frequency => Acos(wt) * Bcos(2wt) = ...
    % => \frac{AB}{2} \cdot cos(3wt) * cos(wt)
    % This will only give us 3rd and 1st harmonics


    help = cos(omega * kgrid.t_array(one_period_indices));

    help_1 = cos(2*omega * kgrid.t_array(one_period_indices));

    help_2 = help .* help_1;
    
    source.p(one_period_indices) = p0 * help_2;
    
    sensor_data_modulated_1st_3rd = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});


dz = dy;

x_vec = (0:(Nx-1))*dx;

z_vec = (0:(Ny-1))*dz;

[ix_cord, iz_cord] = meshgrid(x_vec, z_vec);

% find sensor coordinates based off sensor.mask
coords_sensor_x = zeros(1,size(sensor_data_modulated_1st_3rd,1));
coords_sensor_z = zeros(1,size(sensor_data_modulated_1st_3rd,1));
coords_sensor = [coords_sensor_x;coords_sensor_z];
count = 1;

dt = kgrid.dt;
for i = 1:Nx
    for j = 1:Ny
        if(sensor.mask(i,j) == 1)
            coords_sensor(1, count) = ix_cord(i,j);
            coords_sensor(2, count) = iz_cord(i,j);   
            count = count + 1;
        end
    end
end

% convert coords_sensor from point space to [mm]

% coords_sensor = coords_sensor * dx;


[img_recons_final_modulated_1st_3rd, backprojections_1st_3rd] = dericks_DAS_beamforming(Ny, Nx, size(sensor_data_modulated_1st_3rd,1), ix_cord, iz_cord, coords_sensor,...
    dt, sensor_data_modulated_1st_3rd);

[img_recons_final_modulated_2nd_4th, backprojections_2nd_4th] = dericks_DAS_beamforming(Ny, Nx, size(sensor_data_modulated_2nd_4th,1), ix_cord, iz_cord, coords_sensor,...
    dt, sensor_data_modulated_2nd_4th);

% Create for loop to iterate through all sensors:


sensor_data_modulated_filtered_2nd_gone = zeros(size(sensor_data_modulated_2nd_4th));

for i = 1:size(sensor_data_modulated_2nd_4th, 1)
    % For sensor, i, take the fft
    L = size(sensor_data_modulated_2nd_4th(1,:), 2);
    fft_sensor_data = fft(sensor_data_modulated_2nd_4th(i,:) ,L);
    fft_sensor_data_post_filter = fft_sensor_data;
    % zero out specified indexes:

    % => for 4th harmonics, (1:40) and N signal to 40, assuming a zeropad
    % of 0

    fft_sensor_data_post_filter(1:40) = 0;
    fft_sensor_data_post_filter(L-40+2:end) = 0;
    % take inverse fft
    ifft_sensor_data_filtered = ifft(fft_sensor_data_post_filter);
    sensor_data_modulated_filtered_2nd_gone(i,:) = ifft_sensor_data_filtered;
end

figure
plot(1:1207, abs(fft_sensor_data/L));
title("Original Spectra of Modulated Signal (2nd and 4th Domination)");

figure
plot(1:1207, abs(fft_sensor_data_post_filter/L));
title("Original Spectra of Modulated Signal Second Harmonic Filtered Out");

[img_recons_final_modualtion_2nd_gone, backprojections_modulation_2nd] = dericks_DAS_beamforming(Ny, Nx, size(sensor_data_modulated_filtered_2nd_gone,1), ix_cord, iz_cord, coords_sensor,...
    dt, sensor_data_modulated_filtered_2nd_gone);

sensor_data_modulated_filtered_4th_gone = zeros(size(sensor_data_modulated_2nd_4th));

for i = 1:size(sensor_data_modulated_2nd_4th, 1)
    % For sensor, i, take the fft
    L = size(sensor_data_modulated_2nd_4th(1,:), 2);
    fft_sensor_data = fft(sensor_data_modulated_2nd_4th(i,:) ,L);
    fft_sensor_data_post_filter = fft_sensor_data;
    % zero out specified indexes:

    % => for 4th harmonics, (1:40) and N signal to 40, assuming a zeropad
    % of 0

    fft_sensor_data_post_filter(40:61) = 0;
    first = L-40+2;
    second = L-61+2;
    fft_sensor_data_post_filter(second:first) = 0;
    % take inverse fft
    ifft_sensor_data_filtered = ifft(fft_sensor_data_post_filter);
    sensor_data_modulated_filtered_4th_gone(i,:) = ifft_sensor_data_filtered;
end

figure
plot(1:1207, abs(fft_sensor_data_post_filter/L));
title("Original Spectra of Modulated Signal Fourth Harmonic Filtered Out");

[img_recons_final_modualtion_4th_gone, backprojections_modulation_4th] = dericks_DAS_beamforming(Ny, Nx, size(sensor_data_modulated_filtered_4th_gone,1), ix_cord, iz_cord, coords_sensor,...
    dt, sensor_data_modulated_filtered_4th_gone);


sensor_data_modulated_filtered_1st_gone = zeros(size(sensor_data_modulated_1st_3rd));

for i = 1:size(sensor_data_modulated_1st_3rd, 1)
    % For sensor, i, take the fft
    L = size(sensor_data_modulated_1st_3rd(1,:), 2);
    fft_sensor_data = fft(sensor_data_modulated_1st_3rd(i,:) ,L);
    fft_sensor_data_post_filter = fft_sensor_data;
    % zero out specified indexes:

    % => for 4th harmonics, (1:40) and N signal to 40, assuming a zeropad
    % of 0

    fft_sensor_data_post_filter(3:24) = 0;
    first = L-3+2;
    second = L-24+2;
    fft_sensor_data_post_filter(second:first) = 0;
    % take inverse fft
    ifft_sensor_data_filtered = ifft(fft_sensor_data_post_filter);
    sensor_data_modulated_filtered_1st_gone(i,:) = ifft_sensor_data_filtered;
end

figure
plot(1:1207, abs(fft_sensor_data/L));
title("Original Spectra of Modulated Signal (1st and 3rd Domination)");

figure
plot(1:1207, abs(fft_sensor_data_post_filter/L));
title("Original Spectra of Modulated Signal First Harmonic Filtered Out");

[img_recons_final_modualtion_1st_gone, backprojections_modulation_1st_gone] = dericks_DAS_beamforming(Ny, Nx, size(sensor_data_modulated_filtered_1st_gone,1), ix_cord, iz_cord, coords_sensor,...
    dt, sensor_data_modulated_filtered_1st_gone);





% store in 1st row

% repeat, store in filtered modulation


x_axis_mm = (-(Nx-1)*(dx/2):dx:(Nx-1)*(dx/2));
z_axis_mm = (0:dy:(Ny-1)*dy);

figure
imagesc(x_axis_mm, z_axis_mm, abs(hilbert(img_recons_final_modualtion_1st_gone)));
title("Backprojection Sum w/ 3rd Harmonic Dominant");
xlabel('x [mm]');
ylabel('z [mm]');

figure
imagesc(x_axis_mm, z_axis_mm, abs(hilbert(img_recons_final_modualtion_2nd_gone)));
title("Backprojection Sum w/ 4th Harmonic Dominant");
xlabel('x [mm]');
ylabel('z [mm]');

figure
imagesc(x_axis_mm, z_axis_mm, abs(hilbert(img_recons_final_modualtion_4th_gone)));
title("Backprojection Sum w/ 2nd Harmonic Dominant");
xlabel('x [mm]');
ylabel('z [mm]');


figure
imagesc(x_axis_mm, z_axis_mm, abs(hilbert(img_recons_final_modulated_2nd_4th)));
title("Backprojection Sum w/ both 2nd and 4th Dominant");
xlabel('x [mm]');
ylabel('z [mm]');

figure
imagesc(x_axis_mm, z_axis_mm, abs(hilbert(img_recons_final_modulated_1st_3rd)));
title("Backprojection Sum w/ both 1st and 3rd Dominant");
xlabel('x [mm]');
ylabel('z [mm]');