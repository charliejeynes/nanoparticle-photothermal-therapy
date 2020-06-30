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

bounds = 20 * 1e-3;  % this is the bounds of the world in m^3
sum_absorption_list = []; 
grid_resolution_lst = [];
lines_profiles = cell(number_of_files, 1); 
twoD_slice_lst = cell(number_of_files, 1); 

for i = 1:number_of_files
    
    file = [directory file_names(i).name];
    
    datacube = read_nc(file);
    twoD_slice = get_twoD_slice(datacube);
    twoD_slice_lst{i, 1} = twoD_slice; 
    grid_resolution = get_grid_resolution(twoD_slice); 
    grid_resolution_lst(i) = get_grid_resolution(twoD_slice);
    sum_absorption_list(i) = sum_output_datacube(twoD_slice, bounds, grid_resolution); 
    line_profile = get_line_profile(twoD_slice, grid_resolution); 
    lines_profiles{i, 1} = line_profile; 
    
    
end

figure, 
plot(grid_resolution_lst, sum_absorption_list, 'x')
xlabel('resolution X^2)')
ylabel('sum absorption * cell area [W/m^2)') 

figure, 
plot_twoD_slices(twoD_slice_lst, number_of_files, grid_resolution_lst)

% figure, 
% hold on
% for i=1:number_of_files
% semilogy(1:length(lines_profiles{i, 1}), lines_profiles{i, :})
% legend()
% end



function file_names = get_file_paths(directory)
    
    file_names = dir(directory);
    
    mask = ismember({file_names.name}, {'.', '..', '.DS_Store'});
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

function summ = sum_output_datacube(twoD_slice, bounds, grid_resolution)
    

    scaler = bounds^2 / grid_resolution^2; 
    matrix = twoD_slice * scaler; 
    summ = sum(sum(matrix)); 

%     summ = sum(sum(twoD_slice)) / grid_resolution^2; 
    
end

function line_profile = get_line_profile(twoD_slice, grid_resolution)

    half = round((grid_resolution/2)); 
%     plus = round((grid_resolution/2)) + 5; 
%     minus = round((grid_resolution/2)) - 5); 
    line_profile = twoD_slice(:, half); 

end 

function plot_twoD_slices(twoD_slice_lst, number_of_files, grid_resolution_lst)

        
        figure, 
        for i=1:number_of_files
            
            subplot(2, 3, i); 
            logg = log10(twoD_slice_lst{i, 1}); 
            imagesc(logg); 
            caxis([1,7])% create and label the colorbar
            cmap = jet();
            colormap(cmap);
            cb=colorbar; 
            cb.Label.String = 'Log10 absorption density';
            title(['resolution =', num2str(grid_resolution_lst(i))]);  
                       
        end
      
end