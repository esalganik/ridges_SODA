close all; clear; clc;
lat_soda(1) = 73.00; lon_soda(1) = -148.34; lat_soda(2) = 75.46; lon_soda(2) = -145.64; lat_soda(3) = 77.74; lon_soda(3) = -139.14; % Morring coordinates

% data import
M = readmatrix('201904_UMDRDA_v02_DistanceBetweenRidges.txt'); % Duncan and Farell (2022) http://doi.org/10.5281/ZENODO.7129192
lon_rda = M(:,1); lat_rda = M(:,2); s_rda = M(:,4); % Lon (1), Lat (2), Distance Between Ridges (m) (4)

project = '0.2mThresh_201904_10km_segment_TotalDC_modOIB.nc'; % Mchedlishvili (2023) https://doi.org/10.1594/PANGAEA.959728
s = ncread(project,'x_e'); lat = ncread(project,'lat'); lon = ncread(project,'lon');

t_0 = (datetime('0000-01-01 00:00:00')); % Brenner (2011) http://hdl.handle.net/1773/46919
time = readmatrix('C:\Users\evsalg001\Documents\MATLAB\datasets\Data_for_Brenner_et_al\iceGeometryWeekly.csv','Range','B2:B30'); t_sodaA = t_0 + days(time-1);
time = readmatrix('C:\Users\evsalg001\Documents\MATLAB\datasets\Data_for_Brenner_et_al\iceGeometryWeekly.csv','Range','B56:B98'); t_sodaB = t_0 + days(time-1);
time = readmatrix('C:\Users\evsalg001\Documents\MATLAB\datasets\Data_for_Brenner_et_al\iceGeometryWeekly.csv','Range','B107:B155'); t_sodaC = t_0 + days(time-1);
s_sodaA = readmatrix('C:\Users\evsalg001\Documents\MATLAB\datasets\Data_for_Brenner_et_al\iceGeometryWeekly.csv','Range','L2:L30');
s_sodaB = readmatrix('C:\Users\evsalg001\Documents\MATLAB\datasets\Data_for_Brenner_et_al\iceGeometryWeekly.csv','Range','L56:L98');
s_sodaC = readmatrix('C:\Users\evsalg001\Documents\MATLAB\datasets\Data_for_Brenner_et_al\iceGeometryWeekly.csv','Range','L107:L155');
s_sodaA_april = mean(s_sodaA(27:29)); s_sodaB_april = mean(s_sodaA(25:27)); s_sodaC_april = mean(s_sodaA(26:29)); s_soda_april = mean([s_sodaA_april s_sodaB_april s_sodaC_april]);
s_soda = mean([s_sodaA_april s_sodaB_april s_sodaC_april]); % note: window is 400 km per month based on ice drift velocity

for k = 1:15
    for i = 1:3
        step = k; % averaging over "1+step" 25 km grid cells
        [~,idx] = min(((lat(:)-lat_soda(i)).^2+(lon(:)-lon_soda(i)).^2));
        [idx_2D(1),idx_2D(2)] = ind2sub(size(lat),idx); % GPS index of a mooring
        
        lat_min = min(lat(idx_2D(1)-step:idx_2D(1)+step,idx_2D(2)-step:idx_2D(2)+step),[],'all');
        lat_max = max(lat(idx_2D(1)-step:idx_2D(1)+step,idx_2D(2)-step:idx_2D(2)+step),[],'all');
        lon_min = min(lon(idx_2D(1)-step:idx_2D(1)+step,idx_2D(2)-step:idx_2D(2)+step),[],'all');
        lon_max = max(lon(idx_2D(1)-step:idx_2D(1)+step,idx_2D(2)-step:idx_2D(2)+step),[],'all');
        
        s_rda_filt = s_rda; lat_rda_filt = lat_rda; lon_rda_filt = lon_rda;
        s_rda_filt(lat_rda < lat_min | lat_rda > lat_max | lon_rda < lon_min | lon_rda > lon_max) = NaN;
        lat_rda_filt(lat_rda < lat_min | lat_rda > lat_max | lon_rda < lon_min | lon_rda > lon_max) = NaN;
        lon_rda_filt(lat_rda < lat_min | lat_rda > lat_max | lon_rda < lon_min | lon_rda > lon_max) = NaN;
        lon_rda_filt(isnan(lat_rda_filt)) = []; s_rda_filt(isnan(lat_rda_filt)) = []; lat_rda_filt(isnan(lat_rda_filt)) = [];
        s_rda_filt_all(i) = mean(s_rda_filt,'all');
        s_A_alex(i) = nanmean(s(idx_2D(1)-step:idx_2D(1)+step,idx_2D(2)-step:idx_2D(2)+step),'all'); % Alex, 3x3 grid
    end
    s_rda_soda(k) = mean(s_rda_filt_all);
    s_alex_soda(k) = mean(s_A_alex);
    wgs84 = wgs84Ellipsoid("km");
    scale(k) = stdist(double([lat_min lat_max]),double([lon_min lon_max]),wgs84);
end

figure
% scale = 75:50:50*k+25;
plot(scale,s_rda_soda,'rx'); hold on
plot(scale,s_alex_soda,'bx');
plot(scale,s_soda*ones(length(scale),2),'k--');
hYLabel = ylabel('Ridge spacing (m)'); hXLabel = xlabel('Horizontal scale (km)'); set([hXLabel hYLabel gca],'FontSize',9,'FontWeight','normal'); ylim([0 600]); 
leg = legend('UMD-RDA','ATL07','SODA moorings','box','off'); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];