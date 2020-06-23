clc
clear 
close all

%% test
figure, 
plot(heat_simulation.sensor_data(:, 1:10:400)); 

%%

directory = '/Users/charliejeynes/Projects/dia/sim_data/'; 

file_names = get_file_paths(directory);

for i = 5 %3:size(file_names,1)
    file = [directory file_names(i).name];
    datacube = read_nc(file);
    twoD_slice = get_twoD_slice(datacube); 
    imagesc(twoD_slice);
    [Nx, Ny, dx, dy] = set_grid_heatDiffusion(datacube, 20); % mm_of_world is the dimension of the world in mm set by bounds in ARC
    medium = get_heatDiffusion_constants(Nx, Ny); 
    
    %this is the start of the heat diffusion simulation
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    [input_args, source, sensor] = create_source_sensor(Nx, Ny, twoD_slice);
    % create kWaveDiffusion object
    heat_simulation = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:} ); % inputs for kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});
    % take time steps (temperature can be accessed as kdiff.T)
    heat_simulation.takeTimeStep(400, 1); % .takeTimeStep(Nt, dt) Nt=time in seconds, dt = number of timesteps
    % plot the last temperature field
    figure;
    heat_simulation.plotTemp;
    colorbar;
    
    figure, 
    plot(heat_simulation.sensor_data(:, 300)); 
    
end    

function file_names = get_file_paths(directory)
    
    file_names = dir(directory);
end


function datacube = read_nc(file)
    ncid = netcdf.open(file); 
    datacube = netcdf.getVar(ncid,0); 
    netcdf.close(ncid);
end


function twoD_slice = get_twoD_slice(datacube)
    
    middle = round(size(datacube, 2)/2); 
    twoD_slice = datacube(:, :, middle); 

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
    sensor.mask(:, 100) = 1;


end 


