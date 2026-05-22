%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOADING YOUR DATA FROM MACHINE. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


master_data = load('rf data\master_data_1_21_2026.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD IN RF DATA FROM STRUCT. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all;
% sample_pool = 246:700; % 3_23 data
% sample_pool = 230:700; % 4_5 data


sample_pool = 1:207;

full_orgdata = zeros(size(sample_pool,2),128,1792);
given_rect_images = zeros(size(sample_pool,2),416,512);
given_nsi_images = zeros(size(sample_pool,2),416,512);
for n = sample_pool
    filename = fullfile('data_4_26_v2\',sprintf('%d.mat', n));
    testing_data = load(filename);
    raw_data = testing_data.yeah';
    full_orgdata(n-(sample_pool(1)-1),:,:) = raw_data;
    given_rect_images(n-(sample_pool(1)-1),:,:) = testing_data.rect_img;
    given_nsi_images(n-(sample_pool(1)-1),:,:) = testing_data.nsi_img;
    fprintf('Loaded %d / %d into org. raw_data \n',n, sample_pool(end));
end

% full_orgdata_pt1 = full_orgdata;



sample_pool = 29:100;

full_orgdata = zeros(size(sample_pool,2),128,1792);
given_rect_images = zeros(size(sample_pool,2),416,512);
given_nsi_images = zeros(size(sample_pool,2),416,512);
for n = sample_pool
    filename = fullfile('data_4_18_v3\',sprintf('%d.mat', n));
    testing_data = load(filename);
    raw_data = testing_data.yeah';
    full_orgdata(n-(sample_pool(1)-1),:,:) = raw_data;
    given_rect_images(n-(sample_pool(1)-1),:,:) = testing_data.rect_img;
    given_nsi_images(n-(sample_pool(1)-1),:,:) = testing_data.nsi_img;
    fprintf('Loaded %d / %d into org. raw_data \n',n, sample_pool(end));
end

full_orgdata_pt2 = full_orgdata;

sample_pool = 21:473;

full_orgdata = zeros(size(sample_pool,2),128,1792);
given_rect_images = zeros(size(sample_pool,2),416,512);
given_nsi_images = zeros(size(sample_pool,2),416,512);
for n = sample_pool
    filename = fullfile('data_4_17\',sprintf('%d.mat', n));
    testing_data = load(filename);
    raw_data = testing_data.yeah';
    full_orgdata(n-(sample_pool(1)-1),:,:) = raw_data;
    given_rect_images(n-(sample_pool(1)-1),:,:) = testing_data.rect_img;
    given_nsi_images(n-(sample_pool(1)-1),:,:) = testing_data.nsi_img;
    fprintf('Loaded %d / %d into org. raw_data \n',n, sample_pool(end));
end

full_orgdata_pt3 = full_orgdata;

full_orgdata = cat(1, full_orgdata_pt1, full_orgdata_pt2, full_orgdata_pt3);


full_orgdata = cat(1, full_orgdata_pt1, full_orgdata_pt2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     RUN BEAMFORMING PIPELINE.          %
%  uses verasonics_beamforming_derick_v2  %
%  uses beamforming_verasonics_nsi_v4     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pert = 0.05;
% x

% pert_list = 0.1:0.1:1; % DC offsets [0.1,0.2, ... 1]
% pert_list_1 = [0.04,0.06,0.08];
% pert_list = [pert_list_1, pert_list_2];
% pert_list = [0.1,0.5];
% pert_list = [0:0.2:1];
pert_list = [0.05,0.07,0.09, 0.2,0.4,0.6];
% pert_list = pert_list(2:end);
% pert_list = [0.1,0.3,0.5,0.7];
temp_rect_images = zeros(size(sample_pool,2),416,512);
temp_nsi_images = zeros(size(sample_pool,2),416,512);
final_nsi_images = zeros(size(pert_list,2),416,512);
final_rect_image = zeros(416,512);
focus_index_z = 20e-3;
pert = 0.5;
sample_pool = 1:size(full_orgdata,1);
for k = 1:size(pert_list,2)

    pert = pert_list(k); 
    for n = 1:size(sample_pool,2)
        [curr_nsi_image, curr_rect_image] = verasonics_beamforming_derick_v2(squeeze(full_orgdata(n,:,:)), focus_index_z, pert);
        % [curr_nsi_image, curr_rect_image] = verasonics_beamforming_derick_v2(plane_wave_data, focus_index_z, pert);
        
        if(k == 1)
            temp_rect_images(n,:,:) = curr_rect_image;
        end
        temp_nsi_images(n,:,:) = curr_nsi_image;
        fprintf('%d/%d left, Curr Pert: %d/%d \n', n, size(sample_pool,2), k, size(pert_list,2));
    end
    
    % At this point, temp_rect_images/temp_nsi_images [455 x 416 x 512] 
    % Average the images together, (forgot to divide by N, but its just a scaling factor LOL)

    % Rect Averaging
    if(k == 1)
        avg_img = zeros(416,512);
        for i = 1:size(full_orgdata,1)
            avg_img = avg_img + squeeze(temp_rect_images(i,:,:));
        end
        final_rect_image = avg_img;
    end


    % NSI Averaging
    avg_img = zeros(416,512);
    for i = 1:size(full_orgdata,1)
        avg_img = avg_img + squeeze(temp_nsi_images(i,:,:));
    end
    
    final_nsi_images(k,:,:) = avg_img;
    
end

%%%%%%%%%%%%%%%%%
%% MISC. STUFF %%
%%%%%%%%%%%%%%%%%

final_nsi_images_30mm_real = final_nsi_images_30mm;

final_nsi_images_30mm = ans;

test = zeros(1,416,512);
test(1,:,:) = final_rect_image;
[final_nsi_images; test];

figure
plot(final_rect_image(56:end,255))

figure
plot(final_nsi_images(56:end,255))

figure
i= 138;
i = i + 1;
imagesc(x_vec, z_vec(40:end), squeeze(temp_nsi_images(i,40:end,:)))
colormap gray
% xlim([256-30,256+30])
% figure
% imagesc(final_rect_image(50:end,:))


figure
imagesc(squeeze(final_nsi_images(4,56:end,:)))



figure
i = 0;
i = i + 1;
hm = squeeze(final_nsi_images(i,56:end,:));
hold on;
% figure
plot(z_vec(56:end), squeeze(hm(:,255))/max(squeeze(hm(:,255))))
% 
% figure
% imagesc(to_plot)
% xlim([254,258])
% to_plot_nsi_0_5 = avg_img;
% to_plot_rect = avg_img;
% to_plot_nsi_0_1 = avg_img;
to_plot_nsi_0_05 = avg_img;


figure
imagesc(final_rect_image)

figure
i = 0;
i = i + 1;
imagesc(squeeze(ans(i,:,:)))

test = zeros(1,416,512);

test(1,:,:) = final_rect_image;

[final_nsi_images; test]; 

final_nsi_images_30mm = ans;



resize(final_rect_image, [1,416,512]);

figure
plot(z_vec*1e3,to_plot_rect(:,259)/max(to_plot_rect(:,259)))
hold on;
plot(z_vec*1e3,to_plot_nsi_0_5(:,258)/max(to_plot_nsi_0_5(:,258)))
hold on;
plot(z_vec*1e3,to_plot_nsi_0_1(:,258)/max(to_plot_nsi_0_1(:,258)))
hold on;
plot(z_vec*1e3,to_plot_nsi_0_05(:,258)/max(to_plot_nsi_0_05(:,258)))
hold off;
legend('Rect','DC = 0.5','DC = 0.1','DC = 0.05')
