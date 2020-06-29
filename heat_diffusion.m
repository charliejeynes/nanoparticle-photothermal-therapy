
% Charlie Jeynes - June 2020
% This script contains a number of functions which analyses output from 
% 'ARCTORUS' (see github.com/freddywordingham) 
% and uses this out for input into heat simulations using 'k-wave' (see
% k-wave.org

clc
clear 
close all
%%

% test = simulations{1,1}.plotTemp


%%

% specify your directory of .nc simulation files here
directory = '/Users/charliejeynes/Projects/dia/sim_data/grid_res_test/'; 
% directory = '/Users/charliejeynes/Projects/dia/sim_data/move_back_tumour/'; 


file_names = get_file_paths(directory);
number_of_files = length(file_names); 
%global number_of_files


% setup NOTE the grid resolution is set to 201 - change if this is not the
% case in your simulations. 
heat_profile = zeros(number_of_files,201); % this is a matrix of line profiles through each of the heat simulation files
cem43 = zeros(201, 201, number_of_files); 
lesions = zeros(201, 201, number_of_files); 
kwave_object_simulations = cell(1,number_of_files); 

for i = 1:number_of_files
    
    file = [directory file_names(i).name];
    
    datacube = read_nc(file);
    twoD_slice = get_twoD_slice(datacube);
    
    
    % this is the start of the heat diffusion simulation
    mm_of_world = 20; % mm_of_world is the dimension of the world in mm set by bounds in ARC
    [Nx, Ny, dx, dy] = set_grid_heatDiffusion(datacube, mm_of_world); 
    medium = get_heatDiffusion_constants(Nx, Ny); 
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    [input_args, source, sensor] = create_source_sensor(Nx, Ny, twoD_slice);
    % create kWaveDiffusion object
    heat_simulation = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:} ); % inputs for kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});
    % take time steps (temperature can be accessed as kdiff.T)
    Nt=300; % Nt=time in seconds
    dt=1;   % dt = number of timesteps
    heat_simulation.takeTimeStep(Nt, dt); % .takeTimeStep(Nt, dt) Nt=time in seconds, dt = number of timesteps
    
    plot_heatDiffusion(heat_simulation) % plot the last temperature field
    
    %create kwave simulation object for each file
    kwave_object_simulations{1, i} = heat_simulation; 
    
    % collect up the outputs for each file to compare results later 
%     heat_profile(i, :) = heat_simulation.sensor_data(:, Nt); % add heat profiles to a matrix
%     cem43(:, :, i) = heat_simulation.cem43(:, :); 
%     lesions(:, :, i) = heat_simulation.lesion_map(:, :); 
   
end    

% compare results of all the .nc files in the folder

% plot_heat_profiles(heat_profile)
% image_cem43(cem43, number_of_files)
% image_lesions(lesions, number_of_files)

function file_names = get_file_paths(directory)
    
    file_names = dir(directory);
    
    mask = ismember({file_names.name}, {'.', '..'});
    file_names(mask) = [];   %get rid of . and .. directories
    
    
end


function datacube = read_nc(file)
    ncid = netcdf.open(file); 
    datacube = netcdf.getVar(ncid,0); 
    netcdf.close(ncid);
end


function twoD_slice = get_twoD_slice(datacube)
    
    middle = round(size(datacube, 2)/2); 
    twoD_slice = datacube(:, :, middle); 
    
    figure,
    imagesc(twoD_slice);
    colorbar

end

function [Nx, Ny, dx, dy] = set_grid_heatDiffusion(datacube, mm_of_world)

    res = size(datacube);  % get the resolution (size) of the data cube

    % specify the medium parameters
    % create the computational grid
    Nx = res(1);           % number of grid points in the x (row) direction
    Ny = res(2);           % number of grid points in the y (column) direction
    % Nz = 201; 

    % calculate grid point spacing 
    dx =  (mm_of_world * 1e-3) / res(1) % this is the height (or width / depth) of the 'world' in mm * m / the resolution aka grid points
    dy =  (mm_of_world * 1e-3) / res(2)

    dx = dx;        % grid point spacing in the x direction [m]
    dy = dy;        % grid point spacing in the y direction [m]
    % dz =  1e-4;

end

function medium = get_heatDiffusion_constants(Nx, Ny)

    % define medium properties for the heat diffusion equation

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

end 

function [input_args, source, sensor] = create_source_sensor(Nx, Ny, twoD_slice)

    T0 = 37 .* ones(Nx, Ny); 
    % T0(:, 120:201) = 37 ; 
    %figure, imagesc(T0)
    source.T0 = T0; 

    % testSource = zeros(65,65); 
    % testSource(30:39, 20:30) = 5e6; 
    % source.Q = testSource;

    source.Q = twoD_slice; 
    % set input args
    input_args = {'PlotScale', [37, 50]};

    % set up the sensor 
    sensor.mask = zeros(Nx, Ny);
    sensor.mask(:, :) = 1;
    
end 

function plot_heatDiffusion(heat_simulation)

    figure;
    heat_simulation.plotTemp;
    colorbar;

    
end

function plot_heat_profiles(heat_profile)

    figure, 
    plot(heat_profile');
    title('heat diffusion profiles as the tumour moves back from skin surface by 2mm x 3 times')
    xlabel('resolution (mm*10)')
    ylabel('Temperature ^oC')
    set(gca,'FontSize',20)
    
end

function image_cem43(cem43, number_of_files)
    
%     global number_of_files
    figure, 
    title('Cumulative minutes at 43^o C')
    for i=1:number_of_files
        subplot(1, 3, i)
        imagesc(cem43(:, :, i))
        %caxis([0,7])% create and label the colorbar
        %caxis([0,2500])% create and label the colorbar
        cmap = jet();
        colormap(cmap);
        cb=colorbar; 
        cb.Label.String = 'cem43';
        colorbar
       
    end

end


function image_lesions(lesions, number_of_files)
    
%    global number_of_files
    figure, 
    for i=1:number_of_files
        subplot(1, 3, i)
        imagesc(lesions(:, :, i))
        caxis([0,1])% create and label the colorbar
        cmap = jet(2);
        colormap(cmap);
        cb=colorbar; 
        cb.Label.String = 'thermally ablated';
        colorbar
    end

end


function get_from_sim_object()
    
%    global number_of_files
    figure, 
    for i=1:3
       subplot(1,3,i)
       simulations{1,i}.plotTemp; 
       colorbar
    end
        
end 



