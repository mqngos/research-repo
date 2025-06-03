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
    source_time_signal(one_period_indices) = p0 * sin(2 * pi * source_freq * kgrid.t_array(one_period_indices));
    
    num_source_points = sum(source.p_mask(:)); % should be 1 if single point
    source.p = repmat(source_time_signal, num_source_points, 1);
    
    % % Disc Simulation -- Source
    % source.p0 = 3 * makeDisc(Nx, Ny, Nx/2, Ny/2, 4) + 3 * makeDisc(Nx, Ny, Nx/4, Ny/2, 4) + 3 * makeDisc(Nx, Ny, 3*Nx/4, Ny/2, 4);
    
    
    
    
    % set the input arguements: force the PML to be outside the computational
    % grid; switch off p0 smoothing within kspaceFirstOrder2D
    input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false, 'PlotSim', true};
    
    % run the simulation -- linear
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    
    
    
    % run the simulation -- inversion 
    
    source.p(one_period_indices) = p0 * sin(2 * pi * source_freq * kgrid.t_array(one_period_indices) - pi);
    
    sensor_data_inverted = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});


dz = dy;

x_vec = (0:(Nx-1))*dx;

z_vec = (0:(Ny-1))*dz;

[ix_cord, iz_cord] = meshgrid(x_vec, z_vec);

% find sensor coordinates based off sensor.mask
coords_sensor_x = zeros(1,size(sensor_data,1));
coords_sensor_z = zeros(1,size(sensor_data,1));
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


summed_data = sensor_data + sensor_data_inverted;


% === Define sine wave properties ===

% % Harmonics of Original Wave
f0 = source_freq;          % Frequency of sine wave [Hz]
Fs = 50000000 * (f0);         % Sampling frequency [Hz]
T = 1 / Fs;        % Sampling period [s]
L = size(sensor_data,2);          % Number of samples
t = (0:L-1)*T;     % Time array

L_pad = L * 100;           % zero-pad by factor of 10

Y = fft(sensor_data(1,:), L_pad);        % zero-padded FFT
P2 = abs(Y / L);          % normalize by original length
P1 = P2(1:L_pad/2+1);
P1(2:end-1) = 2 * P1(2:end-1);

f = Fs * (0:(L_pad/2)) / L_pad;

figure;
plot(f, P1)
xlabel('Frequency (MHz)')
ylabel('|P1(f)|')
title('Harmonics of Original Wave')
grid on

% % Harmonics of Original Wave
f0 = source_freq;          % Frequency of sine wave [Hz]
Fs = 50000000 * (f0);         % Sampling frequency [Hz]
T = 1 / Fs;        % Sampling period [s]
L = size(sensor_data,2);          % Number of samples
t = (0:L-1)*T;     % Time array

L_pad = L * 100;           % zero-pad by factor of 10

Y = fft(summed_data(1,:), L_pad);        % zero-padded FFT
P2 = abs(Y / L);          % normalize by original length
P1 = P2(1:L_pad/2+1);
P1(2:end-1) = 2 * P1(2:end-1);

f = Fs * (0:(L_pad/2)) / L_pad;

figure;
plot(f, P1)
xlabel('Frequency (MHz)')
ylabel('|P1(f)|')
title('Harmonics of Summed Wave')
grid on

% From above, we have a first-harmonic dominant signal and a
% second-harmonic dominant signal, so let us modulate these two signals to
% get a third harmonic and first harmonic signal!

%Create a for-loop iteration, start from each sensor, k:

modulated_sensor_data = zeros(size(sensor_data));
debug = zeros(size(source.p));

for k = 1:256
    % Iterate each data point in the sensors
    for i = 1:1207
        % Grab first data point from original wave:
        original_data_point = sensor_data(k,i);
        inverted_data_point = summed_data(k,i);
        modulated_sensor_data(k,i) = original_data_point * inverted_data_point;
        
    end
end


% % Harmonics of Modulated Wave
f0 = source_freq;          % Frequency of sine wave [Hz]
Fs = 50000000 * (f0);         % Sampling frequency [Hz]
T = 1 / Fs;        % Sampling period [s]
L = size(sensor_data,2);          % Number of samples
t = (0:L-1)*T;     % Time array

L_pad = L * 100;           % zero-pad by factor of 10

Y = fft(modulated_sensor_data(1,:), L_pad);        % zero-padded FFT
P2 = abs(Y / L);          % normalize by original length
P1 = P2(1:L_pad/2+1);
P1(2:end-1) = 2 * P1(2:end-1);

f = Fs * (0:(L_pad/2)) / L_pad;

figure;
plot(f, P1)
xlabel('Frequency (MHz)')
ylabel('|P1(f)|')
title('Harmonics of Modulated Wave')
grid on

% Apply a highpass filter after applying modulation, grab the third
% harmonic:

sensor_data_modulated_filtered_2nd_gone = zeros(size(modulated_sensor_data));
%for i = 1:size(sensor_data, 1)
for i = 1:1
    
    % For sensor, i, take the fft
    L = size(modulated_sensor_data(1,:), 2);
    fft_sensor_data = fft(modulated_sensor_data(i,:) ,L);
    fft_sensor_data_post_filter = fft_sensor_data;
    % zero out specified indexes:

    % => for 4th harmonics, (1:40) and N signal to 40, assuming a zeropad
    % of 0

    fft_sensor_data_post_filter(1:26) = 0;
    fft_sensor_data_post_filter(L-26+2:end) = 0;
    % take inverse fft
    ifft_sensor_data_filtered = ifft(fft_sensor_data_post_filter);
    sensor_data_modulated_filtered_2nd_gone(i,:) = ifft_sensor_data_filtered;
end

% Visualize Post-Processing
f0 = source_freq;          % Frequency of sine wave [Hz]
Fs = 50000000 * (f0);         % Sampling frequency [Hz]
T = 1 / Fs;        % Sampling period [s]
L = size(sensor_data,2);          % Number of samples
t = (0:L-1)*T;     % Time array

L_pad = L * 100;           % zero-pad by factor of 10

Y = fft(sensor_data_modulated_filtered_2nd_gone(1,:), L_pad);        % zero-padded FFT
P2 = abs(Y / L);          % normalize by original length
P1 = P2(1:L_pad/2+1);
P1(2:end-1) = 2 * P1(2:end-1);

f = Fs * (0:(L_pad/2)) / L_pad;

figure;
plot(f, P1)
xlabel('Frequency (MHz)')
ylabel('|P1(f)|')
title('Harmonics of Modulated Wave (2nd Harmonic Filtered)')
grid on



[img_recons_final, backprojections] = dericks_DAS_beamforming(Ny, Nx, size(sensor_data,1), ix_cord, iz_cord, coords_sensor,...
    dt, sensor_data);

[img_recons_final_summation, backprojections_summation] = dericks_DAS_beamforming(Ny, Nx, size(sensor_data,1), ix_cord, iz_cord, coords_sensor,...
    dt, summed_data);

[img_recons_final_modulation, backprojections_modulation] = dericks_DAS_beamforming(Ny, Nx, size(sensor_data_modulated_filtered_2nd_gone,1), ix_cord, iz_cord, coords_sensor,...
    dt, sensor_data_modulated_filtered_2nd_gone);



x_axis_mm = (-(Nx-1)*(dx/2):dx:(Nx-1)*(dx/2));
z_axis_mm = (0:dy:(Ny-1)*dy);

figure
imagesc(x_axis_mm, z_axis_mm, abs(hilbert(img_recons_final_modulation)));
title("Original Modulation Sum");
xlabel('x [mm]');
ylabel('z [mm]');


figure
imagesc(x_axis_mm, z_axis_mm, abs(hilbert(img_recons_final_summation)));
title("Pulse Inversion Inversion Backproj. Sum");
xlabel('x [mm]');
ylabel('z [mm]');

% % Hello

% plot the wave as a function of the x-axis (slice x-wise)
figure
plot(1:256, img_recons_final_modulation(128,:));

%Plot the the wave as a function of the z-axis (slice z-wise)
figure
plot(1:256, img_recons_final_modulation(:,128));