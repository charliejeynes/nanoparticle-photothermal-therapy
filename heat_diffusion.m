
% Charlie Jeynes - June 2020
% This script contains a number of functions which analyses output from 
% 'ARCTORUS' (see github.com/freddywordingham) 
% and uses this out for input into heat simulations using 'k-wave' (see
% k-wave.org

clc
clear 
close all
%% test 

% save('with_air_boundary_heat_sim.mat', 'heat_simulation')

% six_minutes = temperature_image{1,1};
% six_minutes = reshape(six_minutes(:, 10), 201, 201); 
% imagesc(six_minutes)
% save('six_minutes.mat', 'six_minutes')

%% input parameters 

% specify your directory of .nc simulation files here

%directory = '/Users/charliejeynes/Projects/dia/sim_data/misc/'; 
%directory = '/Users/charliejeynes/Projects/dia/sim_data/grid_res_test/'; 
%directory = '/Users/charliejeynes/Projects/dia/sim_data/move_back_tumour/'; 
%directory = '/Users/charliejeynes/Projects/dia/sim_data/move_back_tumour/power1W/'; 
%directory = '/Users/charliejeynes/Projects/dia/sim_data/compare_to_hirscht/';
%directory ='/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power0.28/'; 
%directory = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power_changing/'; 
%directory = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power_changing_2mm_back/'; 
%directory = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power_changing_4mm_back/';
%directory = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power_changing_optimised_for_max_tumour_depth/';
%directory = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_tumour4mm_power_changing/'; 
%directory = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_tumour6mm_power_changing/'; 
%directory = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_tumour2mm_power_changing_NO_GNRS/';
%directory = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_tumour2mm_power_changing_NO_GNRS_test/';
%directory = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_tumour4mm_power_changing_NO_GNRS/';
%directory = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_tumour6mm_power_changing_NO_GNRS/';
%directory = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_tumour6mm_power_changing_NO_GNRS_optimised/';
directory = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_tumour6mm_power_changing_with_GNRS_optimised/';

file_names = get_file_paths(directory);
number_of_files = length(file_names); 


% run the simulation
% run_simulation(number_of_files, directory, file_names)

% function run_simulation(number_of_files, directory, file_names)


heat_profile = zeros(number_of_files,201); % this is a matrix of line profiles through each of the heat simulation files
temperature_image = cell(1,number_of_files); 
cem43 = zeros(201, 201, number_of_files); 
lesions = zeros(201, 201, number_of_files); 
kwave_object_simulations = cell(1,number_of_files); 

    for i = 1:number_of_files

        % setup NOTE the grid resolution is set to 201 - change if this is not the
        % case in your simulations.
        
        
        file = [directory file_names(i).name];

        %datacube = read_biomolecules_dotmat(); % this here is to run the
        %biomolecules paper 'just_tumour' data - have to hash the next line
        %if going to run instead
        datacube = read_nc(file);
        twoD_slice = get_twoD_slice(datacube);
        grid_resolution = get_grid_resolution(twoD_slice); 
        grid_resolution_lst(i) = get_grid_resolution(twoD_slice);


        % this is the start of the heat diffusion simulation
        mm_of_world = 20; % mm_of_world is the dimension of the world in mm set by bounds in ARC
        [Nx, Ny, dx, dy] = set_grid_heatDiffusion(datacube, mm_of_world); 
        medium = get_heatDiffusion_constants(Nx, Ny); 
        %unhash the next two lines if you want an heterogenous medium (air/tissue)
%         airStart = 140; 
%         medium = get_heatDiffusion_constants_hetero(Nx, Ny, airStart);  % Use this one if you want a heterogenous medium e.g. tissue/air
        kgrid = kWaveGrid(Nx, dx, Ny, dy);
        [input_args, source, sensor] = create_source_sensor(Nx, Ny, twoD_slice);
        % create kWaveDiffusion object
        heat_simulation = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:} ); % inputs for kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});
        % take time steps (temperature can be accessed as kdiff.T)
        Nt=1; % Nt=number of timesteps 
        dt= 300;   % dt = time in seconds
        heat_simulation.takeTimeStep(Nt, dt); % .takeTimeStep(Nt, dt) 

        plot_heatDiffusion(heat_simulation) % plot the last temperature field

        %create kwave simulation object for each file
        kwave_object_simulations{1, i} = heat_simulation; 



        % collect up the outputs for each file to compare results later 
        % heat_profile(i, :) = heat_simulation.sensor_data(:, Nt); % add heat profiles to a matrix
        cem43(:, :, i) = heat_simulation.cem43(:, :); 
        temperature_image{:, i} = heat_simulation.sensor_data;
        lesions(:, :, i) = heat_simulation.lesion_map(:, :); 

    end    
% end

% get data from the kwaveDiffusion object stored in a cell array
get_from_sim_object(kwave_object_simulations, number_of_files, Nt, grid_resolution_lst)

% compare results of all the .nc files in the folder
plot_heat_profiles(heat_profile)
image_cem43(cem43, number_of_files)
image_lesions(lesions, number_of_files)
plot_line_profile_from_image(cem43, number_of_files)

% get data out of the kwave object files for export to python 
temperature_image_list = extract_temperature_images(temperature_image, number_of_files, grid_resolution); 

% save cem43, legions and temperature final temperature for python to read
save_sims_as_dotmat(cem43, temperature_image_list)


function file_names = get_file_paths(directory)
    
    file_names = dir(directory);
    
    mask = ismember({file_names.name}, {'.', '..', '.DS_Store'});
    file_names(mask) = [];   %get rid of . and .. directories
    
    
end

function datacube = read_biomolecules_dotmat()
    
    datacube = load('just_tumour.mat'); 
    datacube = datacube.absorbCube; 
    datacube = datacube / 2;  % 1.67; 
    
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
%    medium.blood_density(airStart:end, :)    = 0;     % of air [kg/m^3]

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

function medium = get_heatDiffusion_constants_hetero(Nx, Ny, airStart)

    % define medium properties for the heat diffusion equation
    air_density = 12.55;            % 1.255; % of air [kg/m^3]
    air_specific_heat = 718;        % 718 of air [J/(kg.K)]
    air_thermal_conductivity =  26.02;    % 26.02 of air [J/(kg.K)]  
    
    medium.density = ones(Nx, Ny); 
    medium.density(:, :)  = 1079;     % of tissue [kg/m^3]
%     medium.density(:, airStart:end)  = 1.255; % of air [kg/m^3]
    medium.density(airStart:end, :)  = air_density ;  %1.255; % of air [kg/m^3]
    figure, imagesc(medium.density)
    title('medium')

    medium.thermal_conductivity = ones(Nx, Ny); 
    medium.thermal_conductivity(:, :)  = 0.52;     % tissue [W/(m.K)]
%     medium.thermal_conductivity(:, airStart:end)  = 26.02;  % of air [W/(m.K)]
    medium.thermal_conductivity(airStart:end,:) = air_thermal_conductivity;  %26.02;  % of air [W/(m.K)]

    medium.specific_heat = ones(Nx, Ny);
    medium.specific_heat(:, :)   = 3540;     % of tissue [J/(kg.K)]
%     medium.specific_heat(:, airStart:end)  = 0718;        % of air [J/(kg.K)]
    medium.specific_heat(airStart:end, :)  = air_specific_heat;        % of air [J/(kg.K)]

    % define medium properties related to perfusion
    medium.blood_density                = ones(Nx, Ny);     % [kg/m^3]
    medium.blood_density(:, :)          = 1060;     % [kg/m^3]
%     medium.blood_density(:, airStart:end)    = 0;     % [kg/m^3]
    medium.blood_density(airStart:end, :)    = air_density;     % of air [kg/m^3]

    medium.blood_specific_heat          = ones(Nx, Ny);     % [J/(kg.K)]
    medium.blood_specific_heat(:, :)    = 3617;     % [J/(kg.K)]
%     medium.blood_specific_heat(:, airStart:end)    = 0;  % [J/(kg.K)]
    medium.blood_specific_heat(airStart:end, :) = air_specific_heat;  % of air [J/(kg.K)]

    medium.blood_perfusion_rate         = ones(Nx, Ny);    % [1/s]
    medium.blood_perfusion_rate(:, :)   = 0.01;     % [1/s]
%     medium.blood_perfusion_rate(:, airStart:end) = 0;     % [1/s]
    medium.blood_perfusion_rate(airStart:end, :) = 0.01;     % of air [1/s]

    medium.blood_ambient_temperature    = ones(Nx, Ny);       % [degC]
    medium.blood_ambient_temperature(:, :)    = 37;       % [degC]
%     medium.blood_ambient_temperature(:, airStart:end)   = 22;       % [degC]
    medium.blood_ambient_temperature(airStart:end, :)   = 22;      % of air [degC]

end 

function [input_args, source, sensor] = create_source_sensor(Nx, Ny, twoD_slice)

    T0 = 37 .* ones(Nx, Ny); 
%     T0(140:201, :) = 22; 
    figure, imagesc(T0)
    title('source heat');
    colorbar;
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
        subplot(1, number_of_files, i)
        imagesc(cem43(:, :, i))
        %caxis([0,7])% create and label the colorbar
        %caxis([0,2500])% create and label the colorbar
        cmap = jet();
        colormap(cmap);
        cb=colorbar; 
        cb.Label.String = 'cem43';
       
    end

end


function image_lesions(lesions, number_of_files)
    
%    global number_of_files
    figure, 
    for i=1:number_of_files
        subplot(1, number_of_files, i)
        imagesc(lesions(:, :, i))
        caxis([0,1])% create and label the colorbar
        cmap = jet(2);
        colormap(cmap);
        cb=colorbar; 
        cb.Label.String = 'thermally ablated';
        
    end

end


function get_from_sim_object(kwave_object_simulations, number_of_files, Nt, grid_resolution_lst)
    
%    global number_of_files
    
    figure,
    for i=1:number_of_files
       subplot(2,number_of_files,i)
       kwave_object_simulations{1,i}.plotTemp; 
       %caxis([0,7])% create and label the colorbar
       %caxis([0,2500])% create and label the colorbar
       cmap = jet();
       colormap(cmap);
       cb=colorbar; 
%        cb.Label.String = 'Temperature ^oC';
       title(['resolution =', num2str(grid_resolution_lst(i))]); 
       
       
    end
%     currentFigure = gcf;
%     title(currentFigure.Children(end), 'Irradiation for', num2str(Nt) , 'seconds');    
end 


function grid_resolution = get_grid_resolution(twoD_slice)

    grid_resolution = size(twoD_slice, 1); 

end

function save_sims_as_dotmat(cem43, temperature_image_list)

   
%     struct_kwave_object = struct(cem43)
%     save('/Users/charliejeynes/Projects/dia/sim_data/sim_matlab_files/cem43.mat', 'cem43'); 
    save('/Users/charliejeynes/Projects/dia/sim_data/sim_matlab_files/cem43.mat', 'cem43');
    save('/Users/charliejeynes/Projects/dia/sim_data/sim_matlab_files/temperature_image_list.mat', 'temperature_image_list');

end 

function plot_line_profile_from_image(cem43, number_of_files)

    figure, 
    hold on
    for i=1:number_of_files
        semilogy(cem43(:,100, i))
    end
        
end 

function temperature_image_list = extract_temperature_images(temperature_image, number_of_files, grid_resolution)
    
    temperature_image_list = zeros(grid_resolution, grid_resolution, number_of_files); 
    for i=1:number_of_files
        image = temperature_image{1, i}; 
        image_reshaped = reshape(image, grid_resolution, grid_resolution); 
        temperature_image_list(:, :, i) = image_reshaped; 
        
    end
    
end
        
        
function get_cem_boundary(cem43, number_of_files)


    cem_boundaries = []; 
    for i=1:length(number_of_files)
        
        cem_to_mask = cem43(:, :, i); 
        cem_mask = cem_to_mask > 240; 
        figure, imshow(cem_mask)
        cem_boundary = bwboundaries(cem_mask); 
        cem_boundary_non_cell = cem_boundary{1,1}; 
        figure, imshow(cem_boundary{1,1}); 
        dummy = zeros(201, 201); 
        dummy(cem_boundary_non_cell) = 1 ; 
        figure, imshow(dummy); 
        
        imshow(cem_mask); hold on
        plot(cem_boundary_non_cell(:, 2), cem_boundary_non_cell(:, 1), 'lineWidth', 2); 
        
        
    end


end


% function six_minutes = get_temperature_data(temperature_image)
% 
%     six_minutes = temperature_image{1,10}
%     figure, 
%     hold on
%     for i=1:number_of_files
%         imagesc(cem43(:,100, i))
%     end
%         
% end
