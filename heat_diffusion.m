clc
clear 
close all

%% to calculate the beamspot in W/cm^2

beamRadius = 0.3;  % beamspot radius in cm
beamArea = pi*(beamRadius^2)
powerNeeded_1W_cm2 = 1/beamArea  % this is to convert for 1W/cm^2 as is often quioted in the papers

beamArea * powerNeeded

%% to read in data from .nc file (11th June 2020)

% matrix = rand(201,201);
% matrix(matrix(:, 115) > 0.5) = 0;
% 
% minus_NP = data; 
% 
% figure, 
% imagesc(minus_NP(:, :, 32)); 

%% matlab question ask 

% matrix = rand(201,201);
% matrix(matrix(:, 115) > 0.5) = 0; 

%% calcualtion for mask multiplication and final NP OD concentration 

% GNR_OD = 1.07
% tunour_Ua = 2.3 
% tumour_plusAU = 1.07+2.3 
% muliplied = tumour_plusAU * 3
% final_OD_GNRs = 10.1 - 2.3


%% to read data from original run in paper
% Image the absorb cube after loading in manually
% load('/Users/charliejeynes/Projects/git_NP_PTT/just_tumour.mat'); 
% minus_NP = absorbCube;
% figure, 
% % minus_NP = log10(minus_NP); 
% imagesc(rot90(minus_NP(:, :, 100))); %check what it looks like
% %caxis([0,8.5])% create and label the colorbar
% caxis([1e6,5e7])% create and label the colorbar
% cmap = jet();
% colormap(cmap);
% cb=colorbar
% cb.Label.String = 'absorbance density (W/m^3) (log10)';
% 
% absorbCube = absorbCube; % divide by 2 top take account change in spot size
% a = absorbCube(97, 139, 32)
% figure, 
% imagesc(rot90(absorbCube(:, :, 101))); %check what it looks like


%% Hack the absorb cube to give plus and minus
% % REMEMBER TO LOAD THESE IN SEQUENTIALLY
% load('just_tumour.mat'); 
% minus_NP = absorbCube; 
% clear('absorbCube'); 
% 
% %----this can be corrected with the mask------
% % load('3_tumour_plus_Au.mat'); 
% % plus_NP = absorbCube ; % THIS IS FOR THE MASK
% % % plus_NP = absorbCube .* 2.5;   % THIS IS WITHOUT THE MASK it works
% % % plus_NP_2d = plus_NP(: ,:, 100); % THIS IS FOR THE MASK
% % clear('absorbCube'); 
% 
% load('22.5_14.6_e7.mat'); 
% plus_NP = absorbCube ; % THIS IS FOR THE MASK
% % plus_NP = absorbCube .* 2.5;   % THIS IS WITHOUT THE MASK it works
% % plus_NP_2d = plus_NP(: ,:, 100); % THIS IS FOR THE MASK
% clear('absorbCube'); 

%% have a look at the numbers 
% 
% imtool(plus_NP(: ,: ,100))

%%  make the mask over the tumour to correct for the GNR OD

% yes this works to index the tumour!!
% plus_NP_2d = plus_NP(: ,:, 100);
% index = (plus_NP_2d(:, 1:115) > 5e6);
% plus_NP_2d(index) = plus_NP_2d(index) * 4; % THIS IS WHERE THE OD FOR THE NPS IS TRIPLED SO END IS OD = 8
% figure, imagesc(index)
% colorbar
% plus_NP(: ,:, 100) = plus_NP_2d(:, :, 1); 
% % plus_NP(plus_NP(:, 1:115) > 1e5) = 0;
% figure, imagesc(plus_NP(: ,:, 100))
% colorbar

% this works to  put mask to zero
% plus_NP_2d = plus_NP(: ,:, 100);
% plus_NP_2d(plus_NP_2d(:, 1:115) > 5e6) = true;
% % figure, imagesc(plus_NP_2d)
% % colorbar
% plus_NP(: ,:, 100) = plus_NP_2d(:, :, 1); 
% % plus_NP(plus_NP(:, 1:115) > 1e5) = 0;
% figure, imagesc(plus_NP(: ,:, 100))
% colorbar


% %% fig 2a no perfusion 
% 
% % define medium properties
% medium.density              = 1079;     % [kg/m^3]
% medium.thermal_conductivity = 0.52;     % [W/(m.K)]
% medium.specific_heat        = 3540;     % [J/(kg.K)]
% 
% % do minusNP
% minus_NP_noP  = minus_NP(:, :, 100) ./ (medium.density .* medium.specific_heat);
% minus_NP_noP = minus_NP_noP ./ (2); % account for our laser which is 2 times more intense,
% rot_minus_NP_noP = rot90(minus_NP_noP, 1); 
% 
% % do plusNP
% plus_NP_noP  = plus_NP(:, :, 100) ./ (medium.density .* medium.specific_heat);
% plus_NP_noP = plus_NP_noP ./ (2); % account for our laser which is 2 times more intense,
% rot_plus_NP_noP = rot90(plus_NP_noP, 1); 

%% plot all starting no perfusion in subfigure
% start = 0; 
% finish = 3; 
% figure
% hold on
% subplot(1,3,1)
% imagesc(rot_minus_NP_noP)
% caxis([start,finish])% create and label the colorbar
% cb = colorbar;  
% cb.Label.String = 'Heat change, dT, (^oC)';
% title('Control'); 
% 
% subplot(1,3,2)
% imagesc(rot_plus_NP_noP)
% caxis([start,finish])% create and label the colorbar
% cb = colorbar;  
% cb.Label.String = 'Heat change, dT, (^oC)';
% title('With NanoRods (OD=4)'); 
% 
% subplot(1,3,3)
% x = 1:201;
% plot(x, rot_minus_NP_noP(:, 100), x, rot_plus_NP_noP(:, 100)); 
% xlabel('mm')
% ylabel('dT heat change (^oC)')
%  xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% % % xticklabels({'20','18','16','14','12','10','8','6', '4', '2', '0'});
% xticks(0:20:200)
% % set(gca, 'XDir','reverse')
% % set(gca, 'xlim', 200)
% xlim([1 200])
% hold off
%%

clc
clear 
close all
%%

% to read in data from .nc file (11th June 2020)
file = '/Users/charliejeynes/Projects/dia/dia/output/mcrt/absorption_dens.nc'
% source = '/Users/charliejeynes/Projects/dia/dia/output/mcrt/hits.nc'
ncid = netcdf.open(file); 
data = netcdf.getVar(ncid,0); 
netcdf.close(ncid);

%check if data cube is giving the same sum absorption

% oneOone = data; 
 

test = sum(sum(sum(data)))

%%

reslst = [1,11,51,101,301,501]; 
absorption = [a,b,c,d,de,e]

figure, 
plot(reslst, absorption, '-x')

figure, 
semilogy(reslst, absorption, '-x')
ylabel('sum absorption of the datacube')
xlabel('resolution (x^3)')
title('sum absorption of the datacube vs resolution of the grid - 10^4 photons simulated')
%%

res = size(data);  % get the resolution (size) of the data cube

% specify the medium parameters
% create the computational grid
Nx = res(1);           % number of grid points in the x (row) direction
Ny = res(2);           % number of grid points in the y (column) direction
% Nz = 201; 

% calculate grid point spacing 
dx =  (20 * 1e-3) / res(1) % this is the height (or width / depth) of the 'world' in mm * m / the resolution aka grid points
dy =  (20 * 1e-3) / res(2)

dx = dx;        % grid point spacing in the x direction [m]
dy = dy;        % grid point spacing in the y direction [m]
% dz =  1e-4;

kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define medium properties

airStart = 40; 

medium.density = ones(Nx, Ny); 
medium.density(:, :)  = 1079;     % of tissue [kg/m^3]
% medium.density(:, airStart:end)  = 1.255; % of air [kg/m^3]
% medium.density(airStart:end, :)  = 1.255; % of air [kg/m^3]
figure, imagesc(medium.density)

medium.thermal_conductivity = ones(Nx, Ny); 
medium.thermal_conductivity(:, :)  = 0.52;     % tissue [W/(m.K)]
% medium.thermal_conductivity(:, airStart:end)  = 26.02;  % of air [W/(m.K)]
% medium.thermal_conductivity(airStart:end,:)  = 26.02;  % of air [W/(m.K)]

medium.specific_heat = ones(Nx, Ny);
medium.specific_heat(:, :)   = 3540;     % of tissue [J/(kg.K)]
% medium.specific_heat(:, airStart:end)  = 0718;        % of air [J/(kg.K)]
% medium.specific_heat(airStart:end, :)  = 718;        % of air [J/(kg.K)]

% define medium properties related to perfusion
medium.blood_density                = ones(Nx, Ny);     % [kg/m^3]
medium.blood_density(:, :)          = 1060;     % [kg/m^3]
% medium.blood_density(:, airStart:end)    = 0;     % [kg/m^3]
medium.blood_density(airStart:end, :)    = 0;     % of air [kg/m^3]

medium.blood_specific_heat          = ones(Nx, Ny);     % [J/(kg.K)]
medium.blood_specific_heat(:, :)    = 3617;     % [J/(kg.K)]
% medium.blood_specific_heat(:, airStart:end)    = 0;  % [J/(kg.K)]
% medium.blood_specific_heat(airStart:end, :)    = 0;  % of air [J/(kg.K)]

medium.blood_perfusion_rate         = ones(Nx, Ny);    % [1/s]
medium.blood_perfusion_rate(:, :)   = 0.01;     % [1/s]
% medium.blood_perfusion_rate(:, airStart:end) = 0;     % [1/s]
% medium.blood_perfusion_rate(airStart:end, :) = 0;     % of air [1/s]

medium.blood_ambient_temperature    = ones(Nx, Ny);       % [degC]
medium.blood_ambient_temperature(:, :)    = 37;       % [degC]
% medium.blood_ambient_temperature(:, airStart:end)   = 22;       % [degC]
% medium.blood_ambient_temperature(airStart:end, :)   = 22;      % of air [degC]

minus_NP = data; % divide by 2 top take account change in spot size

figure,
subplot(1, 2, 1)
log_minus_NP = log10(minus_NP); 
% imagesc(rot90(minus_NP(:, :, 100))); %check what it looks like
imagesc(log_minus_NP(:, :, (round(res(3)/2)))); %check what it looks like
caxis([0,7])% create and label the colorbar
% caxis([1e6,5e7])% create and label the colorbar
cmap = jet();
colormap(cmap);
cb=colorbar
cb.Label.String = 'absorbance density (W/m^3) (log10)';

% a = minus_NP(57, 33, 32)
subplot(1,2,2) 
% imagesc(rot90(minus_NP(:, :, 32),2));
imagesc(minus_NP(:, :, (round(res(3)/2))));%check what it looks like
colorbar

% get kdiff for the control 1 SECOND

T0 = 37 .* ones(Nx, Ny); 
% T0(:, 120:201) = 37 ; 
figure, imagesc(T0)
source.T0 = T0; 

% testSource = zeros(65,65); 
% testSource(30:39, 20:30) = 5e6; 
% source.Q = testSource;

source.Q = (minus_NP(:, : ,(round(res(3)/2)))); 
% set input args
input_args = {'PlotScale', [37, 50]};

% set up the sensor 
sensor.mask = zeros(Nx, Ny);
sensor.mask(22, :) = 1;

% create kWaveDiffusion object
kdiff_control_1 = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});

% % take time steps (temperature can be accessed as kdiff.T)
Nt = 600; 
dt = 1;
kdiff_control_1.takeTimeStep(Nt, dt);
% % 
% % plot the current temperature field
figure;
kdiff_control_1.plotTemp;
colorbar;

%% check if res is giving different outputs

% oneOone = data; 
sixtyfive = data; 
% test = oneOone(50, :, 50); 

figure, 
plot(1:101, oneOone(:, 50, 50), 1:65, sixtyfive(:, 50, 50))
legend('100x100x100 res', '65x65x65 res')




%% get kdiff for the control 600 SECONDs (10 minutes)
T0 = 37 .* ones(201, 201); 
T0(:, 120:201) = 37 ; 
figure, imagesc(T0)
source.T0 = T0; 
source.Q = (minus_NP(:, : ,100)./2); % account for our laser which is 2 times more intense,
% set input args
input_args = {'PlotScale', [37, 55]};

% set up the sensor 
sensor.mask = zeros(Nx, Ny);
sensor.mask(100, :) = 1;

% create kWaveDiffusion object
kdiff_control_600 = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});

% % take time steps (temperature can be accessed as kdiff.T)
Nt = 600; 
dt = 1;
kdiff_control_600.takeTimeStep(Nt, dt);
% % 
% % plot the current temperature field

figure;
kdiff_control_600.plotTemp;
colorbar;

%% get kdiff for the NPs 1 SECOND
T0 = 37 .* ones(201, 201); 
T0(:, 120:201) = 37 ; 
figure, imagesc(T0)
source.T0 = T0; 
source.Q = (plus_NP(:, : ,100)./2); % account for our laser which is 2 times more intense,
% set input args
input_args = {'PlotScale', [37, 55]};

% set up the sensor 
sensor.mask = zeros(Nx, Ny);
sensor.mask(100, :) = 1;

% create kWaveDiffusion object
kdiff_NP_1 = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});

% % take time steps (temperature can be accessed as kdiff.T)
Nt = 1; 
dt = 1;
kdiff_NP_1.takeTimeStep(Nt, dt);
% % 
% % plot the current temperature field
figure;
kdiff_NP_1.plotTemp;
colorbar;

%% get kdiff for the NPs 600 SECONDs
T0 = 37 .* ones(201, 201); 
T0(:, 120:201) = 37 ; 
figure, imagesc(T0)
source.T0 = T0; 
source.Q = (plus_NP(:, : ,100)./2); % account for our laser which is 2 times more intense,
% set input args
input_args = {'PlotScale', [37, 55]};

% set up the sensor 
sensor.mask = zeros(Nx, Ny);
sensor.mask(100, :) = 1;

% create kWaveDiffusion object
kdiff_NP_600 = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});

% % take time steps (temperature can be accessed as kdiff.T)
Nt = 600; 
dt = 1;
kdiff_NP_600.takeTimeStep(Nt, dt);
% % 
% % plot the current temperature field

figure;
kdiff_NP_600.plotTemp;
colorbar;

%%
imagesc(kdiff_NP_600.T)
colorbar

%% THIS IS THE START OF PLOTTING ALL THE DATA
%-------- plot in the middle of the tumour every second------

% time_series = kdiff_control_600.sensor_data(110, :); 
% 
% figure
% plot(1:600, time_series)

real_data_direct_control_x = [0 20 40 60 80 100 120 240 400 500 600]; 
real_data_direct_control_y = [0 2.5 4 4.7 5 5.6 5.8 5.9 6.5 6.5 6.5];

real_data_direct_NP_x = [0 20 40 60 80 100 120 240 400 500 600]; 
real_data_direct_NP_y = [0 11 15 18 19.5 21 21.5 22.5 22.5 22.5 22.5]; 
 
control_time_series_600 = kdiff_control_600.sensor_data(110, :) - 37; 
NP_time_series_600 = kdiff_NP_600.sensor_data(110, :) - 37; % the 110 here should relate to the column 

fh = figure; 
hold on

lineProps.width = 1;
lineProps.edgestyle = ':';
lineProps.style = '-';
% lineProps.col = 'b';
% figure; title('Custom line properties')
err =  2.5*ones(size(real_data_direct_control_y)); %[4 4 4 4 4 2 2 2 2 2 2]
mseb(real_data_direct_control_x, real_data_direct_control_y, err, lineProps); % control data
lineProps.width = 1;
lineProps.edgestyle = '--';
lineProps.style = '--';
mseb(real_data_direct_NP_x, real_data_direct_NP_y,err,lineProps);% GNR injected data

% errorfill(real_data_direct_control_x, real_data_direct_control_y, 0.25, 'b-+')
% errorfill(real_data_direct_NP_x, real_data_direct_NP_y, 0.25, 'b-+')

plot(1:600, control_time_series_600, '-.r','LineWidth', 5)
plot(1:600, NP_time_series_600, '--.r','LineWidth', 5)

%     real_data_direct_control_x, real_data_direct_control_y, 'x', real_data_direct_NP_x, real_data_direct_NP_y, 'x')
 
% errorbar(real_data_direct_control_x,real_data_direct_control_y,err)
% errorbar(real_data_direct_NP_x,real_data_direct_NP_y,err)


ylabel('Heat rise, dT (^oC)'); 
legend('Control data','NanoRod data','Control simulation' , 'NanoRod simulation' )
xlabel('Time (seconds)')
% title('Heat rise with constant irradiation vs time')
ylim([0 25]); 
supersizeme(fh, 2.5);
%% plot out a line profile of the sensor mask 600 seconds

% control_profile_600 = kdiff_control_600.sensor_data(:, 600) - kdiff_control_600.sensor_data(:, 1); 
% NP_profile_600 = kdiff_NP_600.sensor_data(:, 600) - kdiff_NP_600.sensor_data(:, 1); % the 110 here should relate to the column 
% 
% figure
% hold on
% plot(1:201, control_profile_600, 1:201, NP_profile_600)
% xlabel('mm')
% ylabel('dT heat change (^oC)')
%  xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% % % xticklabels({'20','18','16','14','12','10','8','6', '4', '2', '0'});
% xticks(0:20:200)
% % set(gca, 'XDir','reverse')
% % set(gca, 'xlim', 200)
% xlim([1 200])
% title('Line profile');
% hold off

% %% plot out a line profile of the sensor mask 600 seconds
% 
% control_profile_600 = kdiff_control_600.sensor_data(:, 600) - kdiff_control_600.sensor_data(:, 1); 
% NP_profile_600 = kdiff_NP_600.sensor_data(:, 600) - kdiff_NP_600.sensor_data(:, 1); % the 110 here should relate to the column 
% 
% figure
% hold on
% plot(1:201, control_profile_600, 1:201, NP_profile_600)
% xlabel('mm')
% ylabel('dT heat change (^oC)')
%  xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% % % xticklabels({'20','18','16','14','12','10','8','6', '4', '2', '0'});
% xticks(0:20:200)
% % set(gca, 'XDir','reverse')
% % set(gca, 'xlim', 200)
% xlim([1 200])
% title('Line profile');
% hold off

%% Figure 2 make figure of W/m3
start = 1e5; 
finish = 2e7; 

fh = figure; 
hold on
subplot(1,3,1)
imagesc(rot90((minus_NP(:, :, 100)/2), 1))
yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
ylabel('mm')
yticks(0:20:200)
xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
xticks(0:20:200)
xlabel('mm')
caxis([start,finish])% create and label the colorbar
cb = colorbar;  
cb.Label.String = 'energy deposited  (W/m^3)';
title('Control');

subplot(1,3,2)
imagesc(rot90((plus_NP(:, :, 100)./2), 1))
yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
ylabel('mm')
yticks(0:20:200)
xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
xticks(0:20:200)
xlabel('mm')
caxis([start,finish])% create and label the colorbar
cb = colorbar;  
cb.Label.String = 'energy deposited  (W/m^3)';
title('With NanoRods'); 

subplot(1,3,3)
x = 1:201;
plot(x, rot90((minus_NP(100, :, 100)./2), 1), x, rot90((plus_NP(100, :, 100)./2), 1)); 
xlabel('mm')
ylabel('energy deposited  (W/m^3)')
 xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% % xticklabels({'20','18','16','14','12','10','8','6', '4', '2', '0'});
xticks(0:20:200)
% set(gca, 'XDir','reverse')
% set(gca, 'xlim', 200)
xlim([1 200])
title('Line profile');
hold off

supersizeme(fh, 2.5);

%% this is heat change after 1 second 
% this is the no perfusion bit

% r = 1
% c = 3
% start = 0; 
% finish = 4.5; 
% fh = figure; 
% hold on
% subplot(r,c,1)
% imagesc(rot_minus_NP_noP)
% yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% ylabel('mm')
% yticks(0:20:200)
% xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% xticks(0:20:200)
% xlabel('mm')
% caxis([start,finish])% create and label the colorbar
% cb = colorbar;  
% cb.Label.String = 'Heat change, dT, (^oC)';
% title('Control'); 
% 
% subplot(r,c,2)
% imagesc(rot_plus_NP_noP)
% yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% ylabel('mm')
% yticks(0:20:200)
% xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% xticks(0:20:200)
% xlabel('mm')
% caxis([start,finish])% create and label the colorbar
% cb = colorbar;  
% cb.Label.String = 'Heat change, dT, (^oC)';
% title('With NanoRods'); 
% 
% subplot(r,c,3)
% x = 1:201;
% plot(x, rot_minus_NP_noP(:, 100), x, rot_plus_NP_noP(:, 100)); 
% xlabel('mm')
% ylabel('heat change, dT, (^oC)')
%  xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% % % xticklabels({'20','18','16','14','12','10','8','6', '4', '2', '0'});
% xticks(0:20:200)
% % set(gca, 'XDir','reverse')
% % set(gca, 'xlim', 200)
% xlim([1 200])
% title('Line profile');
% 
% supersizeme(fh, 2.5);

%% figure 3. this is the perfusion bitplot out kdiff_control and kdiff_NP 

% after 1 second
r = 2
c = 3

start = 0; 
finish = 3; 

fh = figure; 
subplot(r,c,1)
imagesc(rot90((kdiff_control_1.T - 37),1));
yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
ylabel('mm')
yticks(0:20:200)
xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
xticks(0:20:200)
xlabel('mm')
caxis([start,finish])% create and label the colorbar
cb = colorbar;  
cb.Label.String = 'Heat change, dT, (^oC)';
title('Control'); 

subplot(r,c,2)
imagesc(rot90((kdiff_NP_1.T - 37),1))
yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
ylabel('mm')
yticks(0:20:200)
xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
xticks(0:20:200)
xlabel('mm')
caxis([start,finish])% create and label the colorbar
cb = colorbar;  
cb.Label.String = 'Heat change, dT, (^oC)';
title('With NanoRods'); 

subplot(r,c,3)
rotated_control_1 = rot90((kdiff_control_1.T - 37),1); 
rotated_NP_1 = rot90((kdiff_NP_1.T - 37),1); 
x = 1:201;
plot(x, rotated_control_1(:, 100), x, rotated_NP_1(:, 100)); 
xlabel('mm')
ylabel('dT heat change (^oC)')
 xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% % xticklabels({'20','18','16','14','12','10','8','6', '4', '2', '0'});
xticks(0:20:200)
% set(gca, 'XDir','reverse')
% set(gca, 'xlim', 200)
xlim([1 200])
title('Line profile');
hold off

% after ten minutes 

start = 0; 
finish = 25; 

subplot(r,c,4)
imagesc(rot90((kdiff_control_600.T - 37),1));
yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
ylabel('mm')
yticks(0:20:200)
xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
xticks(0:20:200)
xlabel('mm')
caxis([start,finish])% create and label the colorbar
cb = colorbar;  
cb.Label.String = 'Heat change, dT, (^oC)';
title('Control'); 

subplot(r,c,5)
imagesc(rot90((kdiff_NP_600.T - 37),1))
yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
ylabel('mm')
yticks(0:20:200)
xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
xticks(0:20:200)
xlabel('mm')
caxis([start,finish])% create and label the colorbar
cb = colorbar;  
cb.Label.String = 'Heat change, dT, (^oC)';
title('With NanoRods'); 

subplot(r,c,6)
rotated_control_600 = rot90((kdiff_control_600.T - 37),1); 
rotated_NP_600 = rot90((kdiff_NP_600.T - 37),1); 
x = 1:201;
plot(x, rotated_control_600(:, 100), x, rotated_NP_600(:, 100)); 
xlabel('mm')
ylabel('dT heat change (^oC)')
 xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% % xticklabels({'20','18','16','14','12','10','8','6', '4', '2', '0'});
xticks(0:20:200)
% set(gca, 'XDir','reverse')
% set(gca, 'xlim', 200)
xlim([1 200])
title('Line profile');
hold off

supersizeme(fh, 2.5);

%% check the unroted plot is correct

% figure, plot(x, (kdiff_control_600.T(100, :) - 37), x, (kdiff_NP_600.T(100, :) - 37))

%% plot the sensor data as a line profile 
% figure;
% plot(1:201, kdiff.sensor_data(:, 360), 1:201, kdiff.sensor_data(:, 1));
% 
% dT = kdiff_control_600.sensor_data(:, :) - kdiff_control_600.sensor_data(:, 1);
% 
% dT = kdiff.sensor_data; 
% figure, 
% plot(1:201, dT(:, 1), 1:201, dT(:, 20), 1:201, dT(:, 60), 1:201, dT(:, 360), 1:201, dT(:, 600)); 
% xlabel('mm')
% ylabel('dT heat change (^oC)')
% % xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% xticklabels({'20','18','16','14','12','10','8','6', '4', '2', '0'});
% xticks(0:20:200)
% set(gca, 'XDir','reverse')


%% Figure 3: Visulaise the CEM43 

r = 2
c = 3

% plot out the CEM43 after 10 minutes

start = 0; 
finish = 3000; 

rotated_cem43_control_600 = rot90((kdiff_control_600.cem43 / 60),1); 
rotated_cem43_NP_600 = rot90((kdiff_NP_600.cem43 / 60),1); 

fh = figure; 
subplot(r,c,1)
imagesc(rotated_cem43_control_600);
yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
ylabel('mm')
yticks(0:20:200)
xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
xticks(0:20:200)
xlabel('mm')
caxis([start,finish])% create and label the colorbar
cb = colorbar;  
cb.Label.String = 'CEM43 (mins)';
title('Control'); 

subplot(r,c,2)
imagesc(rotated_cem43_NP_600)
yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
ylabel('mm')
yticks(0:20:200)
xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
xticks(0:20:200)
xlabel('mm')
caxis([start,finish])% create and label the colorbar
cb = colorbar;  
cb.Label.String = 'CEM43 (mins)';
title('With NanoRods '); 

subplot(r,c,3) 
x = 1:201;
plot(x, rotated_cem43_control_600(:, 100), x, rotated_cem43_NP_600(:, 100)); 
xlabel('mm')
ylabel('CEM43 (mins)')
 xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% % xticklabels({'20','18','16','14','12','10','8','6', '4', '2', '0'});
xticks(0:20:200)
% set(gca, 'XDir','reverse')
% set(gca, 'xlim', 200)
xlim([1 200])
title('Line profile');
hold off


% --------------plot out estimated survival fraction after ten minutes 

% -------------survival fraction at 43 with time on cancer cells in vitro

survival_fraction = [100 65 50 25 15 7 4 3]'; % this is the cancer data at 43 degrees vs time
time_hours = [0 60 90 120 150 210 240 300]'; 
f = fit(time_hours,survival_fraction,'exp1');  % this is the exponetial function , 'StartPoint',[100,3]
f(300)

% p = polyfix(time_hours,survival_fraction,2,[0,10000],[100,3]); 
% figure, plot(time_hours,survival_fraction,'.',time_hours,polyval(p,time_hours)); 
% p
% 
% x=310
% f = p(1)*x^2 + p(2)*x + p(3)
% f(x)
% -------example on polyfix----------

% x = linspace(0,2,100)';y = sin(pi*x)+ 0.1*randn(100,1);
% p = polyfix(x,y,3,[0,2],[0,0]);figure, plot(x,y,'.',x,polyval(p,x));
% -------------------------------------


%this was to check that the fitting function was OK
% figure, semilogy(time_hours, survival_fraction, 'x') 
% ylim([1 100])
% 
% figure, plot(time_hours, survival_fraction, 'x')
% ylim([1 100])

% plot(f,time_hours,survival_fraction)
% ans = f(150)


% ---------------- convert our CEM43 to survival fraction 
% -------this bit rotates data and get the SF using exponential function
rotated_SF_control_600 = f(rot90((kdiff_control_600.cem43 / 60),1)); % / 60 from seconds to minutes
rotated_SF_control_600 =  100 - reshape(rotated_SF_control_600, 201, 201); 
rotated_SF_NP_600 = f(rot90((kdiff_NP_600.cem43 / 60),1)); % / 60 from seconds to minutes
rotated_SF_NP_600 = 100 -reshape(rotated_SF_NP_600, 201, 201); 
% ---------------------------------------------------------------------
% --------this bit plots out the % cell kill
start = 0; 
finish = 100; 

subplot(r,c,4)
imagesc(rotated_SF_control_600);
yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
ylabel('mm')
yticks(0:20:200)
xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
xticks(0:20:200)
xlabel('mm')
caxis([start,finish])% create and label the colorbar
cb = colorbar;  
cb.Label.String = 'estimated % dead cells';
title('Control'); 

subplot(r,c,5)
imagesc(rotated_SF_NP_600)
yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
ylabel('mm')
yticks(0:20:200)
xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
xticks(0:20:200)
xlabel('mm')
caxis([start,finish])% create and label the colorbar
cb = colorbar;  
cb.Label.String = 'estimated % dead cells';
title('With NanoRods '); 

subplot(r,c,6)
x = 1:201;
plot(x, rotated_SF_control_600(:, 100), x, rotated_SF_NP_600(:, 100)); 
xlabel('mm')
ylabel('estimated % dead cells')
xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% % xticklabels({'20','18','16','14','12','10','8','6', '4', '2', '0'});
xticks(0:20:200)
% set(gca, 'XDir','reverse')
% set(gca, 'xlim', 200)
xlim([1 200])
ylim([0 110])
title('Line profile');
hold off

% start = 0; 
% finish = 40;
% figure
% hold on
% subplot(1,2,1)
% imagesc(kdiff_control_600.cem43)
% caxis([start,finish])% create and label the colorbar
% cb = colorbar;
% subplot(1,2,2)
% imagesc(kdiff_NP_600.cem43)
% 
% figure
% plot(1:201, kdiff_control_600.cem43(100, :)); 


supersizeme(fh, 2.5); 



%% -----------put a mask around the % dead cells to for percentage normal cell death -------- 

% --------draw a red line around the tumour out line
 
% % plus_NP_2d = plus_NP(: ,:, 100);
% % imagesc(plus_NP_2d)
% % yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% % ylabel('mm')
% % yticks(0:20:200)
% % xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% % xticks(0:20:200)
% % xlabel('mm')
% % caxis([start,finish])% create and label the colorbar
% % cb = colorbar;  
% % cb.Label.String = 'estimated % dead cells';
% % title('With NanoRods '); 
% 
% % --------make the mask TRUE
% rot_plus_NP_2d = rot90(plus_NP(: ,:, 100),1);
% mask = false(size(rot_plus_NP_2d)); % rot it to fit with the image
% index = (rot_plus_NP_2d(:, 1:115) > 6e6);
% mask(index) = true;
% figure, imagesc(mask)
% 
% % --------display the image
% % imagesc(rotated_SF_NP_600)
% yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% ylabel('mm')
% yticks(0:20:200)
% xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% xticks(0:20:200)
% xlabel('mm')
% caxis([start,finish])% create and label the colorbar
% cb = colorbar;  
% cb.Label.String = 'estimated % dead cells';
% title('With NanoRods ');
% 
% hold on
% visboundaries(mask,'Color','r');
% hold off

%%
%------get the number of voxels of normal tissue thaty have some death in them
% --------make the mask FALSE
% --------superimpose that on the data 
% rot_plus_NP_2d = rot90(plus_NP(: ,:, 100),1);
% mask = true(size(rot_plus_NP_2d)); % rot it to fit with the image
% index = (rot_plus_NP_2d(:, 1:115) > 6e6);
% mask(index) = false;
% figure, imagesc(mask)

%%
%------get the number of voxels of normal tissue thaty have some death in them
% --------make the mask FALSE
% --------superimpose that on the data 
start = 0
finish = 100

fh = figure 
subplot(1,2,1)
rot_plus_NP_2d = rot90(plus_NP(: ,:, 100),1);
mask = true(size(rot_plus_NP_2d)); % rot it to fit with the image
index = (rot_plus_NP_2d(:, 1:115) > 6e6);
rotated_SF_NP_600(index) = false;

imagesc(rotated_SF_NP_600)
yticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
ylabel('mm')
yticks(0:20:200)
xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
xticks(0:20:200)
xlabel('mm')
caxis([start,finish])% create and label the colorbar
cb = colorbar;  
cb.Label.String = 'estimated % dead cells';
title('With NanoRods '); 

subplot(1,2,2)
x = 1:201;
plot(x, rotated_SF_NP_600(:, 100)); 
xlabel('mm')
ylabel('estimated % dead cells')
xticklabels({'0','2','4','6','8','10','12','14', '16', '18', '20'});
% % xticklabels({'20','18','16','14','12','10','8','6', '4', '2', '0'});
xticks(0:20:200)
% set(gca, 'XDir','reverse')
% set(gca, 'xlim', 200)
xlim([1 200])
ylim([0 110])
title('Line profile');


% start = 0; 
% finish = 40;
% figure
% hold on
% subplot(1,2,1)
% imagesc(kdiff_control_600.cem43)
% caxis([start,finish])% create and label the colorbar
% cb = colorbar;
% subplot(1,2,2)
% imagesc(kdiff_NP_600.cem43)
% 
% figure
% plot(1:201, kdiff_control_600.cem43(100, :)); 

supersizeme(fh, 2.5);
%%
% display the image but with the tumour == 0

SF_mask_zero = rotated_SF_NP_600; 
SF_mask_zero = SF_mask_zero(mask); 
figure, imagesc(SF_mask_zero)







%%

% %% 2d from the example NO air added
% % create the computational grid
% Nx = 201;           % number of grid points in the x (row) direction
% Ny = 201;           % number of grid points in the y (column) direction
% % Nz = 201; 
% dx = 1e-4;        % grid point spacing in the x direction [m]
% dy = 1e-4;        % grid point spacing in the y direction [m]
% % dz =  1e-4;
% kgrid = kWaveGrid(Nx, dx, Ny, dy);
% 
% % define medium properties
% medium.density = 1079;     % of tissue [kg/m^3]
% medium.thermal_conductivity = 0.52;     % tissue [W/(m.K)
% medium.specific_heat        = 3540;     % [J/(kg.K)]
% 
% % define medium properties related to perfusion
% medium.blood_density                = 1060;     % [kg/m^3]
% medium.blood_specific_heat          = 3617;     % [J/(kg.K)]
% medium.blood_perfusion_rate         = 0.01;     % [1/s]
% medium.blood_ambient_temperature    = 37;       % [degC]
% 
% % set Gaussian initial temperature distribution [degC]
% % width = 4 * dx;
% % source.T0 = 37 + 5 .* exp( -(kgrid.x ./ width).^2 - (kgrid.y ./width).^2 );
% T0 = 37 .* ones(201, 201); 
% T0(:, 120:201) = 20 ; 
% figure, imagesc(T0)
% source.T0 = T0; 
% source.Q = (absorbCube(:, : ,100)); 
% % set input args
% input_args = {'PlotScale', [37, 60]};
% 
% % set up the sensor 
% sensor.mask = zeros(Nx, Ny);
% sensor.mask(100, :) = 1;
% 
% % create kWaveDiffusion object
% kdiff = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});
% 
% % % take time steps (temperature can be accessed as kdiff.T)
% Nt = 300; 
% dt = 0.5;
% kdiff.takeTimeStep(Nt, dt);
% % % 
% % % plot the current temperature field
% figure;
% kdiff.plotTemp;
% colorbar;
% 
% %% calculate using bioExact
% 
% % define medium properties
% medium.density = 1079;     % of tissue [kg/m^3]
% medium.thermal_conductivity = 0.52;     % tissue [W/(m.K)
% medium.specific_heat        = 3540;     % [J/(kg.K)]
% 
% % define medium properties related to perfusion
% medium.blood_density                = 1060;     % [kg/m^3]
% medium.blood_specific_heat          = 3617;     % [J/(kg.K)]
% medium.blood_perfusion_rate         = 0.01;     % [1/s]
% medium.blood_ambient_temperature    = 37;       % [degC]
% 
% % calculate perfusion coefficient from the medium parameters
% P = medium.blood_density .* medium.blood_perfusion_rate .* ...
%     medium.blood_specific_heat ./ (medium.density .* medium.specific_heat);
% 
% % calculate diffusivity from the medium parameters
% D = medium.thermal_conductivity / (medium.density * medium.specific_heat);
% 
% % Initial temperature 
% T0 = 37 .* ones(201, 201); 
% T0(:, 120:201) = 20 ;
% 
% % calculate normalised heat source
% S = (absorbCube(:, :, 100)./8) ./ (medium.density .* medium.specific_heat);
% 
% S(100, 100)
% 
% % compute Green's function solution using bioheatExact
% T_exact = bioheatExact(T0, S, ...
%     [D, P, medium.blood_ambient_temperature], 1e-4, 360);
% 
% 
% start = 37; 
% finish = 40; 
% figure
% hold on
% % subplot(1,3,1)
% imagesc(T_exact)
% caxis([start,finish])% create and label the colorbar
% cb = colorbar;  
% cb.Label.String = 'final heat (^oC)';
% title('tumour (Ua=2.3, Us=21.2)'); 
% 
% 
% %% 3D that has an error
% 
% % create the computational grid
% Nx = 201;           % number of grid points in the x (row) direction
% Ny = 201;           % number of grid points in the y (column) direction
% Nz = 201; 
% dx = 1e-4;        % grid point spacing in the x direction [m]
% dy = 1e-4;        % grid point spacing in the y direction [m]
% dz =  1e-4;
% kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
% 
% % define medium properties
% medium.density              = 1079;     % [kg/m^3]
% medium.thermal_conductivity = 0.52;     % [W/(m.K)]
% medium.specific_heat        = 3540;     % [J/(kg.K)]
% 
% % set Gaussian initial temperature distribution [degC]
% % width = 4 * dx;
% % source.T0 = 37 + 5 .* exp( -(kgrid.x ./ width).^2 - (kgrid.y ./width).^2 );
% T0 = 37 .* ones(201, 201, 201); 
% T0(:, 120:201, :) = 20 ; 
% source.T0 = T0; 
% source.Q = absorbCube; 
% % set input args
% input_args = {'PlotScale', [37, 50]};
% 
% % create kWaveDiffusion object
% kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args{:});
% 
% % % take time steps (temperature can be accessed as kdiff.T)
% Nt = 300;
% dt = 0.5;
% kdiff.takeTimeStep(Nt, dt);
% % % 
% % % plot the current temperature field
% figure;
% kdiff.plotTemp;
% 
% 
% % % calculate diffusivity from medium parameters
% % D = medium.thermal_conductivity / (medium.density * medium.specific_heat);
% % 
% % % compute Green's function solution using bioheatExact
% % T_exact = bioheatExact(source.T0, 0, [D, 0, 0], kgrid.dx, Nt * dt);