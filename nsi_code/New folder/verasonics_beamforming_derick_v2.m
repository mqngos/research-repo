% master_data.Trans.spacing = 1.0510
function [nsi_image, rect_image] = verasonics_beamforming_derick_v2(data, focus_index_z, pert)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize L9-4 Transducer Thingies... %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trans_connections = [
1 3 67 8 9 11 75 16 17 19 83 24 25 27 91 32 ...
33 35 99 40 41 43 107 48 49 51 115 56 57 59 123 64 ...
2 4 68 72 10 12 76 80 18 20 84 88 26 28 92 96 ...
34 36 100 104 42 44 108 112 50 52 116 120 58 60 124 128 ...
66 5 6 70 74 13 14 78 82 21 22 86 90 29 30 94 ...
98 37 38 102 106 45 46 110 114 53 54 118 122 61 62 126 ...
65 69 7 71 73 77 15 79 81 85 23 87 89 93 31 95 ...
97 101 39 103 105 109 47 111 113 117 55 119 121 125 63 127 ]';

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Image Medium  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Organizing our RF data  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

organized_raw_data = zeros(128, 1792);


raw_data = data;

for j = 1:128
    organized_raw_data(j,:) = raw_data(trans_connections(j),:);
end

% organized_raw_data = data;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Upsampling / RF Data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

U = 4;
x  = 1:size(organized_raw_data,2);
y  = organized_raw_data;
xq = 1:1/U:size(organized_raw_data,2);
yq = spline(x,y,xq);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call Beamforming Logic / RF Data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[bfm_image] = beamforming_verasonics_nsi_v4(Nz, Nx, 128, ix_cord, iz_cord, sensor_location, U*Fs, ...
    yq, 5, size(yq,2), 32, c0, 3,0, focus_index_z,-1 );
[zm] = beamforming_verasonics_nsi_v4(Nz, Nx, 128, ix_cord, iz_cord, sensor_location, U*Fs, ...
    yq, 5, size(yq,2), 32, c0, 0,pert, focus_index_z,-1 );
[dc1] = beamforming_verasonics_nsi_v4(Nz, Nx, 128, ix_cord, iz_cord, sensor_location, U*Fs, ...
    yq, 5, size(yq,2), 32, c0, 1,pert, focus_index_z,-1 );
[dc2] = beamforming_verasonics_nsi_v4(Nz, Nx, 128, ix_cord, iz_cord, sensor_location, U*Fs, ...
    yq, 5, size(yq,2), 32, c0, 2,pert, focus_index_z,-1 );

zm = abs(hilbert(zm));
dc1 = abs(hilbert(dc1));
dc2 = abs(hilbert(dc2));

nsi_image = ((dc1+dc2)./2) - zm;

rect_image = abs(hilbert(bfm_image));

% [a,b] = (max(rect_image(:,:)));
% 
% [c,d] = max(a);
% % close([1,2]);
% figure;
% imagesc(rect_image)
% 
% figure
% plot(256-3:256+3,a)
% ylim([0,c])
% fprintf("\nindex of maxmimum lateral: %d\n", d+(256-3-1));


[a,b] = (max(nsi_image(:,:)));

[c,d] = max(a);
% close([1,2]);
% figure
% plot(a)

% figure;
% imagesc(nsi_image)

% figure
% plot(a)
% ylim([0,c])
% fprintf("\nindex of maxmimum lateral: %d\n", d);
% 
% fprintf('done!\n')
% figure
% imagesc(nsi_image(:,232:282))
% title('nsi bruh')


% 0.0001