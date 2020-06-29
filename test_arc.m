% Charlie Jeynes - June 2020
% This script contains a number of functions which analyses output from 
% 'ARCTORUS' (see github.com/freddywordingham) for sanity checking 


clc
clear 
close all


%%

% specify your directory of .nc simulation files here
directory = '/Users/charliejeynes/Projects/dia/sim_data/grid_res_test/'; 
% directory = '/Users/charliejeynes/Projects/dia/sim_data/move_back_tumour/'; 


file_names = get_file_paths(directory);
number_of_files = length(file_names); 


sum_absorption_list = []; 
grid_resolution_lst = [];
lines_profiles = cell(number_of_files, 1); 

for i = 1:number_of_files
    
    file = [directory file_names(i).name];
    
    datacube = read_nc(file);
    twoD_slice = get_twoD_slice(datacube);
    grid_resolution = get_grid_resolution(twoD_slice); 
    grid_resolution_lst(i) = get_grid_resolution(twoD_slice);
    sum_absorption_list(i) = sum_output_datacube(twoD_slice, grid_resolution); 
    line_profile = get_line_profile(twoD_slice, grid_resolution); 
    lines_profiles{i, 1} = line_profile; 
    
    
end

figure, 
semilogy(sum_absorption_list)

% figure, 
% hold on
% for i=1:number_of_files
% semilogy(1:length(lines_profiles{i, 1}), lines_profiles{i, :})
% legend()
% end



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

function grid_resolution = get_grid_resolution(twoD_slice)

    grid_resolution = size(twoD_slice, 1); 

end

function summ = sum_output_datacube(twoD_slice, grid_resolution)
    
    summ = sum(sum(twoD_slice)) / grid_resolution^2; 
    
end

function line_profile = get_line_profile(twoD_slice, grid_resolution)

    half = round((grid_resolution/2)); 
%     plus = round((grid_resolution/2)) + 5; 
%     minus = round((grid_resolution/2)) - 5); 
    line_profile = twoD_slice(:, half); 

end 