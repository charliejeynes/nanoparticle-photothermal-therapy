#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 10:24:19 2020

@author: charliejeynes

This gets data from ARCTORUS .nc output files and visualises and analyses the 
results 


"""

import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import numpy as np
import pandas as pd
import netCDF4
import sys
import os
import scipy.io as spio
import cv2 as cv2
from skimage.measure import label, regionprops
import glob




plt.close('all')

#filepath = '/Users/charliejeynes/Projects/dia/sim_data/move_back_tumour/power1W/'

#filepath = '/Users/charliejeynes/Projects/git_NP_PTT/'

#filepath = '/Users/charliejeynes/Projects/git_NP_PTT/six_minutes.mat'

#filepath = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power0.28/'

# arc_nc_filepath = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power0.28/'

#arc_nc_filepath = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power_changing/'

#arc_nc_filepath = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power_changing_2mm_back/'

# arc_nc_filepath = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power_changing_4mm_back/'

#arc_nc_filepath = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power_changing_optimised_for_max_tumour_depth/'

# arc_nc_filepath = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_tumour1mm_power_changing/'

#arc_nc_filepath = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_tumour4mm_power_changing/'; 

# arc_nc_filepath = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_tumour6mm_power_changing/'; 

arc_nc_filepath = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_tumour2mm_power_changing_NO_GNRS/'

matlab_files_path = '/Users/charliejeynes/Projects/dia/sim_data/sim_matlab_files/'

x = 12
max_depth_tumour_mm = [x, x, x]
max_tumour_depth = 120
power = [0.3, 0.5, 0.7]
# max_tumour_depth = 100
# power = [0.3, 0.5, 0.7]
# max_tumour_depth = 80
# power = [0.6, 0.8, 1]

# name_to_save_file = '0mm_back_5mins_heat'
# name_to_save_file = '2mm_back_5mins_heat'
# name_to_save_file = '4mm_back_5mins_heat'
name_to_save_file = 'spot6mm_1mmtumour_5mins_heat'


def list_nc_files(filepath):
    '''
    lists the file path for all .nc simulation files in the subdirectory
    
    '''

    files_list = []
    for file in os.listdir(filepath):
        if file.startswith('.'):
            file = []
        else:
            fullpath = os.path.join(filepath, file)
            files_list.append(fullpath)
            
    return sorted(files_list)



def get_full_file_path(filepath, data):
    
    ''' join the filepathto the .mat file name
    '''
    
    return os.path.join(filepath, data)



def read_nc(file):
    
    '''
    
    Parameters
    ----------
    files_list : LIST OF STRINGS
        reads in datacubes from .nc file output in a subdirectory  

    Returns
    -------
    cross_sections : LIST OF 2D ARRAYS
        returns a 2D cross section of each datacube, returning a list

    '''

    nc = netCDF4.Dataset(file)
    print(nc.variables.keys()) # get all variable names
    datacube = nc.variables['data']  # access variable
    print(datacube)
    datacube = datacube[:,:,:] 
    datacube = datacube.data   # this gets the vakues from the masked array
    nc.close()
    return  datacube


def read_in_heat_simulation(filepath):
    
    return spio.loadmat(filepath, squeeze_me=True)


def get_data_from_spio_dic(data, name_of_data):
    
    ''' this passes in the 'data' which is a 3d matrix from matlab, and the name 
    of the data as a string , so that it can look it up in the dict
    
    '''
    return data[name_of_data]

       
def convert_3Dmatrix_to_lst(heat_sims_matrix):
    
    '''
    The 3Ddimension in the matrix is data from different simulations. It was easier to store 
    and save in this for at in matlab. BUt now the 3D dimension is converted to a list. 
    Takes a 3D matlab array and makes the 3rd dimension a list in an n*m array. 
    This is then in the right format for the rest of the functions in this
    scripts

    '''
    
    lst_of_heat_sims = []
    for third_dimension in range(heat_sims_matrix.shape[-1]):
        
        #profile = np.mean(cem43[:, 90:110, i], axis=0)
        heat_sim = heat_sims_matrix[:, :, third_dimension]
        lst_of_heat_sims.append(heat_sim)
            
    return lst_of_heat_sims   
    

def get_grid_resolution():
    print('to do')

def convert_resolution_to_mm():
    print('to do')
    
def resolution_in_mm_vector(): 
    resolution_in_mm = np.arange(0,201, 1) / 10
    # print(resolution_in_mm, resolution_in_mm.size) 
    return resolution_in_mm




def get_2D_section(datacube, slice_pos_in_y):
    
    section = datacube[:, slice_pos_in_y, :]    
    
    return section 



def get_sections_as_lst(files_list):
    cross_sections = []
    for file in files_list:
        datacube = read_nc(file)
        section = get_2D_section(datacube, 100)
    
        cross_sections.append(section)
    return cross_sections
    
 

def log_data(cross_sections):
    
    '''
    
    Parameters
    ----------
    cross_sections : LIST OF 2D ARRAYS
        list of 2D numpy arrays (cross sections from a 3D datacube).

    Returns
    -------
    log_cross_sections : LIST OF 2D ARRAYS
        log10 transformation of each array in the list.
        
    '''
    log_cross_sections = []
    for section in cross_sections:
        print(section)
        logged = np.log10(section)
        log_cross_sections.append(logged)   
    return log_cross_sections    


def image_2D_sections(cross_sections, log_cross_sections):
    '''

    Parameters
    ----------
    cross_sections : list of 2D array
        
    log_cross_sections : list of 2D array
        
    Returns
    -------
    figure with the 2D array as an image.

    '''
    fig, ax = plt.subplots(3, 2, sharex='col', sharey='row')
    # axes are in a two-dimensional array, indexed by [row, col]
    for r in range(3):
        for c in range(2):
            if c%2 == 0:
                img = ax[r, c].imshow(cross_sections[r], cmap='jet', vmin=1e5, vmax=1e7)
                fig.colorbar(img, ax=ax[r,c])
            else:
                img = ax[r, c].imshow(log_cross_sections[r],cmap='jet')
                fig.colorbar(img, ax=ax[r,c])




def get_line_profile(cross_sections):
    
    line_profiles = []
    # try:
    #     for section in cross_sections:
    #         profile = np.mean(section[90:110, :], axis=1)
    #         line_profiles.append(profile)
    #         plt.plot(profile)
       
    # except: 
    #     for section in cross_sections:
    #         #profile = np.mean(cem43[:, 90:110, i], axis=0)
    #         profile = section[130, :]
    #         plt.plot(profile)
    #         line_profiles.append(profile)     
    for section in cross_sections:
        profile = np.mean(section[90:110, : ], axis=0)
        line_profiles.append(profile)
        #plt.plot(profile)
    
    return line_profiles
        
 
 
def plot_image_profile(images, profiles, logScale = True):
    
    '''
    Parameters
    ----------
    cross_sections : list of 2D array
        
    log_cross_sections : list of 2D array
        
    Returns
    -------
    figure with the 2D array as an image next to the line profile in y.

    '''
    
    fig, ax = plt.subplots(3, 2, figsize=(20,13))
    fig.tight_layout()
    # axes are in a two-dimensional array, indexed by [row, col]
    for r in range(3):
        for c in range(2):
            if c%2 == 0:
                img = ax[r, c].imshow(images[r], cmap='jet')#, vmin=1, vmax=7)
                fig.colorbar(img, ax=ax[r,c])
            elif logScale == True:
                ax[r, c].semilogy(profiles[r])
                ax[r, c].axis(ymin=10,ymax=1e7)
            # elif plot_cem_data == True:
            #     ax[r, c].semilogy(cem_data_line_profiles[r])
            #     ax[r, c].axis(ymin=10,ymax=1e7)
            else:
                ax[r, c].plot(profiles[r])
                ax[r, c].axis(ymin=37,ymax=70)
                
    
                
                
def plot_hirsch_data():
    
    depth_from_skin = [0,  0.5, 1,    2,  3,   4,   5,   6] # in mm 
    control_1min =    [5,  6, 5,    4,  3,  2.5,  2,   1] # heat rise in C
    control_6min =    [12, 13.5, 11.5, 9,  9,  8.0,  7.5, 6] # heat rise in C

    # biomolecules_data = [8.2, 6.4, 4.9, 3.7, 2.7, 1.8, 1.4]

    plt.figure(1)
    plt.plot(depth_from_skin, control_6min, '-b', label= 'Hirscht et al') 
    # plt.plot(depth_from_skin, biomolecules_data, '-r', label= 'biomolecules paper')
    plt.ylim(0, 15)
    plt.ylabel('Temperature rise (C)')
    plt.xlabel('Depth from skin surface (mm)')
    plt.legend(loc="upper right")
    

def covert_spotsize_to_power(spot_radius = 0.25, power_quoted = 4):
    
    ''' this takes the spot diameter (usually ~5mm) in mm  and power quoted 
    in a paper (W/cm^2) and converts it to the 
    the power in the spot
    '''

    #spot_radius = 0.25 # cm
    spot_area = 3.14 / 4 * spot_radius**2
    #power_quoted = 4 # W/cm^2
    
    power_in_spot = power_quoted / spot_area
    return power_in_spot

    
def rotate_2D_section(lst_of_heat_sims):
    
    rotated_heat_sims = []
    for each_simulation in lst_of_heat_sims:
        rotated_heat_sims.append(np.rot90(each_simulation))
    return rotated_heat_sims


def scale_beam_spot_in_arc():
    
    spot = 10 # in m
    scaleby = 0.5e-3 # to scale to mm
    scaled_spot = spot * scaleby 


def plot_arc_nc_files(filepath):
    files_list = list_nc_files(filepath)
    cross_sections = get_sections_as_lst(files_list)
    log_cross_sections = log_data(cross_sections)   
    image_2D_sections(cross_sections, log_cross_sections)
    line_profiles = get_line_profile(cross_sections) 
    plot_image_profile(log_cross_sections, line_profiles)
arc_figs = plot_arc_nc_files(arc_nc_filepath)

def get_arc_nc_data(filepath):
    files_list = list_nc_files(filepath)
    cross_sections = get_sections_as_lst(files_list)
    log_cross_sections = log_data(cross_sections)   
    line_profiles = get_line_profile(cross_sections) 
    return log_cross_sections, line_profiles, cross_sections, log_cross_sections
log_cross_sections, line_profiles, cross_sections, log_cross_sections = get_arc_nc_data(arc_nc_filepath)


def plot_kwave_sim_files(matlab_files_path, cem_or_temp_rise_dot_mat, cem_or_temp_rise, logScale = False):
    matlab_files_path = get_full_file_path(matlab_files_path, cem_or_temp_rise_dot_mat)
    spio_dict = read_in_heat_simulation(matlab_files_path)
    matlab_file = get_data_from_spio_dic(spio_dict, cem_or_temp_rise)
    matlab_data_lst = convert_3Dmatrix_to_lst(matlab_file)
    rot_data = rotate_2D_section(matlab_data_lst)
    #log_cem43 = log_data(rot_lst_cem43)
    #image_2D_sections(lst_cem43, log_cem43)
    data_line_profiles = get_line_profile(rot_data) 
    plot_image_profile(rot_data, data_line_profiles, logScale)
    
cem_or_temp_rise_dot_mat = 'temperature_image_list.mat'
cem_or_temp_rise = 'temperature_image_list'

plot_kwave_sim_files(matlab_files_path, 'temperature_image_list.mat', 'temperature_image_list', logScale = False)
plot_kwave_sim_files(matlab_files_path, 'cem43.mat', 'cem43', logScale = True)

def get_kwave_sim_data(matlab_files_path, cem_or_temp_rise_dot_mat, cem_or_temp_rise):
    matlab_files_path = get_full_file_path(matlab_files_path, cem_or_temp_rise_dot_mat)
    spio_dict = read_in_heat_simulation(matlab_files_path)
    matlab_file = get_data_from_spio_dic(spio_dict, cem_or_temp_rise)
    matlab_data_lst = convert_3Dmatrix_to_lst(matlab_file)
    rot_data = rotate_2D_section(matlab_data_lst)
    #log_cem43 = log_data(rot_lst_cem43)
    #image_2D_sections(lst_cem43, log_cem43)
    data_line_profiles = get_line_profile(rot_data) 
    log_data_line_profiles = log_data(data_line_profiles)
    return matlab_data_lst, rot_data, data_line_profiles, log_data_line_profiles

cem_data_lst, cem_rot_data, cem_data_line_profiles, log_data_line_profiles = get_kwave_sim_data(
                                                           matlab_files_path, 'cem43.mat', 'cem43')


def plot_biomoleculesPaper_data(filepath):
    six_minutes = read_in_heat_simulation(filepath)
    six_minutes = get_data_from_spio_dic(six_minutes, 'six_minutes')
    line_profile = six_minutes[100, :]
    line_profile_norm = line_profile - 37
    line_profile_rev = line_profile_norm[::-1]
    # # plt.figure(3)
    # plt.plot(np.arange(1, 202, 1), line_profile_rev)
    # # plt.figure(4)
    # plt.plot(np.arange(0, 20, 0.1), line_profile_rev[0:200])
    plt.figure(1)
    plt.plot(np.arange(0, 6, 0.1), line_profile_rev[80:140], '-r', label='biomolecules paper')
    plt.legend()
    # plt.plot((np.arange(0, 5.5, 0.05), line_profile_norm[:110])


def get_survival_fitting_function():
    '''
    plot out estimated survival fraction after ten minutes 

    survival fraction at 43 with time on cancer cells in vitro
    
    '''
    survival_fraction = [100, 65, 50, 25, 15, 7, 4, 3] # % this is the cancer data at 43 degrees vs time
    time_min = [0, 60, 90, 120, 150, 210, 240, 300] 
    # f = fit(time_min,survival_fraction,'exp1');  # this is the exponetial function , 'StartPoint',[100,3]
    # f(300)
 
    
    
def cem_apply_threshold(cem_rot_data, threshold_number):
    
    thresholded_cem_mask_lst = []
    for cem in range(len(cem_rot_data)):
        sim = cem_rot_data[cem] 
        mask = sim > 240
        # sim[mask] = 10
        # sim[~mask] = 20
        thresholded_cem_mask_lst.append(mask)
    return thresholded_cem_mask_lst
thresholded_cem_mask_lst = cem_apply_threshold(cem_rot_data, 240) 



def make_boundary_from_mask(thresholded_cem_mask_lst, where_air_starts):
    
    cem_xy_boundaries = []
    for mask in thresholded_cem_mask_lst:   
        image = np.ascontiguousarray(mask, dtype=np.uint8)
        contours, hierarchy = cv2.findContours(image,cv2.RETR_TREE,
                                               cv2.CHAIN_APPROX_SIMPLE)[-2:]
        xy = np.array(contours)
        boundary = np.squeeze(xy)
        if boundary.size > 0: # only this if there is a cem boundary, don't if empty
            truncated_boundary = np.where(boundary[:, 0] < where_air_starts)
            cem_xy_boundaries.append(boundary[truncated_boundary]) 
        else:
            cem_xy_boundaries.append(np.array([])) #this empty list is to keep track of the files
       
    return cem_xy_boundaries
cem_xy_boundaries = make_boundary_from_mask(thresholded_cem_mask_lst, 140)




def plot_absorbDensity_and_cem_profiles(arc_nc_filepath, cem_xy_boundaries, max_depth_tumour_mm, power):
    
    [absorbDensity_cross_sections, 
     absorbDensity_line_profiles,
     cross_sections, 
     log_cross_sections] =  get_arc_nc_data(arc_nc_filepath)
    
    cem_data_lst, cem_rot_data, 
    cem_data_line_profiles, 
    log_data_line_profiles = get_kwave_sim_data(matlab_files_path, 'cem43.mat', 'cem43')
    
    fig, ax = plt.subplots(3, 2, figsize=(20,13))
    #fig.tight_layout()
    x=resolution_in_mm_vector()
    # axes are in a two-dimensional array, indexed by [row, col]
    for r in range(3):
        for c in range(2):  
            if c%2 == 0: # this images the 2d cross sections of simulations, with cem43 lesion overlayed  
                img = ax[r, c].imshow(log_cross_sections[r], cmap='jet', vmin=0.1, vmax=7.5)
                ax[r, c].tick_params(axis='both', which='major', labelsize=20)
                ax[r, c].tick_params(axis='both', which='minor', labelsize=8)

                if cem_xy_boundaries[r].size: # this checks if the list is empty and if so, ignores it
                    ax[r,c].plot(cem_xy_boundaries[r][:, 0], cem_xy_boundaries[r][:, 1], linewidth=5, color='blue')
                cbar = fig.colorbar(img, ax=ax[r,c])
                cbar.ax.set_ylabel('absorption density (W/m$^3$)', fontsize=20, labelpad= 25)
                cbar.ax.tick_params(labelsize=20)
                # cbar.ax.set_yticklabels([1,3,5,7])
                
                
                # Where we want the ticks, in pixel locations
                ticks = np.arange(0, 201, 1)
                # What those pixel locations correspond to in data coordinates.
                # Also set the float format here
                ticklabels = ["{:6.0f}".format(i) for i in ticks/10]
                # ticklabels = ticks/10
                
                
                ax[r, c].set_xticks(ticks)
                ax[r, c].set_xticklabels(ticklabels)
                ax[r, c].set_yticks(ticks)
                ax[r, c].set_yticklabels(ticklabels)

                

                n = 40  # Keeps every nth label
                [l.set_visible(False) for (i,l) in enumerate(ax[r, c].xaxis.get_ticklabels()) if i % n != 0]
                [l.set_visible(False) for (i,l) in enumerate(ax[r, c].yaxis.get_ticklabels()) if i % n != 0]
                
                ax[2, 0].set_xlabel('(mm)', fontsize=25, labelpad = 0)
                ax[r, c].set_ylabel('(mm)', fontsize=25, labelpad = 0)
                
                ax[r, c].set_title('Power = %1.2f' %power[r], fontsize=25)
                
            else: 
                ax[r, c].semilogy(resolution_in_mm_vector(), absorbDensity_line_profiles[r], linewidth=5, color='black', label= 'absorption density (W/m$^3$)' )
                ax[r, c].semilogy(resolution_in_mm_vector()[:140], cem_data_line_profiles[r][:140], linewidth=5, color='blue', label='CEM43')
                ax[r, c].semilogy(resolution_in_mm_vector(), np.repeat(240, 201), '--', linewidth=5, color='purple', label='ablation threshold (CEM43>240)')
                ax[r, c].semilogy(np.repeat(max_depth_tumour_mm[r], 7), np.array([1,1e2,10e3, 10e4, 10e5, 10e6, 10e7]), '--', linewidth=5, color='red', label='max depth of tumour')
                ax[r, c].axis(ymin=10,ymax=3e7)
                ax[0, 1].legend(loc='upper left', fontsize=17)
                ax[r, c].tick_params(axis='both', which='major', labelsize=20)
                # Make a plot with major ticks that are multiples of 20 and minor ticks that
                # are multiples of 5.  Label major ticks with '%d' formatting but don't label
                # minor ticks.
                ax[r, c].xaxis.set_major_locator(MultipleLocator(2))
                ax[r, c].xaxis.set_major_formatter(FormatStrFormatter('%d'))
                # For the minor ticks, use no labels; default NullFormatter.
                ax[r, c].xaxis.set_minor_locator(MultipleLocator(0.2))
                ax[2, 1].set_xlabel('(mm)', fontsize=25)
    plt.tight_layout()

plot_absorbDensity_and_cem_profiles(arc_nc_filepath, cem_xy_boundaries, max_depth_tumour_mm, power)  
          



def cem_at_max_tumour_depth(max_tumour_depth, power):
    
    '''this measures the cem43 value at the very back of the tumour, and then fits an exponential 
    to the points, so that the optmised power can be extapolated, which is predicted to ablated the 
    tumour, but with minimal tissue damage'''
    
    cem_data_lst, cem_rot_data, cem_data_line_profiles, log_data_line_profiles = get_kwave_sim_data(
                                                           matlab_files_path, 'cem43.mat', 'cem43')
    
     # reolution of the grid
    cem_at_max_tumour_depth = []
    
    for cem in cem_data_line_profiles:
        cem_at_max_tumour_depth.append(cem[max_tumour_depth])
    
    fig, ax = plt.subplots(1,1)    
    ax.semilogy(power, cem_at_max_tumour_depth, 'x')
    ax.set_xlabel('Power (W)')
    ax.set_ylabel('CEM43 value at max depth of tumour')
    
    x = np.array(power)
    y = np.array(cem_at_max_tumour_depth)

    p = np.polyfit(x, np.log(y), 1)
    ax.semilogy(x, np.exp(p[0] * x + p[1]), 'g--')
    ax.legend(['data','exponential fit'])
    plt.show(fig)
    
    
    optimised_power = (np.log(240) - p[1] ) / p[0]
    return optimised_power    
optimised_power = cem_at_max_tumour_depth(max_tumour_depth, power) #
print('The optimised power is', optimised_power)      
    
def combine_figs_as_pdf():
    filepath = '/Users/charliejeynes/Projects/git_NP_PTT/saved_figures/'
    print('to do')

def optimised_power_and_time_to_ablate_tumour_with_depth():
    
    '''
    max tumour depth is in unit of the resolution of the grid
    
    power is a list of the powers (floats in Watts) used for each simulation

    '''
    
    max_depth_of_tumour = [0, 2, 4] #this is in mm
    ten_mins_optimised_power_to_ablate_tumour = [0.19, 0.35, 0.77] # these numbers are the results of running this script on the simulations changing the tumour depth
    five_mins_optimised_power_to_ablate_tumour = [0.21, 0.39, 0.865]
    two_mins_optimised_power_to_ablate_tumour = [0.24, 0.499, 1.15]
    
    
    
    fig, ax = plt.subplots(1, 1)
    ax.semilogy(max_depth_of_tumour, ten_mins_optimised_power_to_ablate_tumour, '-x', color='blue', label='ten minutes irradiance')
    ax.semilogy(max_depth_of_tumour, five_mins_optimised_power_to_ablate_tumour, '-x', color='red', label='five minutes irradiance')
    ax.semilogy(max_depth_of_tumour, two_mins_optimised_power_to_ablate_tumour, '-x', color='green', label='two minutes irradiance')
    ax.set_ylabel('optimised_power_to_ablate_tumour (W)')
    ax.set_xlabel('max_depth_of_tumour (mm)')
    ax.legend()
optimised_power_and_time_to_ablate_tumour_with_depth()  
    
def ratio_ablated_normal_tissue_to_tumour(thresholded_cem_mask_lst, radius_of_tumour=10):
    
    '''
    **only works for the optimised simualtions**
    this takes the thresholded cem mask list for the 3 optimised laser powers to ablate to the back of the tumour
    and plots against the depth of the tumour (top most surface)
    The aim is to show how much normal tissue is ablated with respect to depth the destroy the tumour. 
    
    '''
    
    #tumour_2d_area = 3.14 * (radius_of_tumour ** 2)
    
    num_pixels_in_tumour = 3.14 * (radius_of_tumour ** 2)
    
    num_pixel_in_cem_ablation = []
    for cem in thresholded_cem_mask_lst:
        num_pixel_in_cem_ablation.append(np.sum(cem))
        
    ratio = []
    for ablation in num_pixel_in_cem_ablation:
        ratio.append(ablation / num_pixels_in_tumour)
    
    max_depth_of_tumour = [0, 2, 4] 
    
    fig, ax = plt.subplots(1, 1)
    ax.semilogy(max_depth_of_tumour, ratio, 'x')
    ax.set_ylabel('areal ratio of ablated normal tissue \n to tumour tissue')
    ax.set_xlabel('max_depth_of_tumour (mm)')
ratio_ablated_normal_tissue_to_tumour(thresholded_cem_mask_lst, radius_of_tumour=10)
    
def save_ablation_masks_per_depth(name_to_save_file):
    np.save(name_to_save_file, thresholded_cem_mask_lst)
save_ablation_masks_per_depth(name_to_save_file)    

# def open_saved_ablation_masks_into_dict():
#     ablated_masks_dict = {}
#     for np_name in glob.glob('*.np[yz]'):
#         ablated_masks_dict[np_name] = np.load(np_name)

# def unpack_ablation_mask_dict():
#     lst = []
#     for name in ablated_masks_dict:
#         for i in range(3):
#             lst.append(sum(sum(ablated_masks_dict[name][1, : ,:])))

def get_npz_filenames_lst():
    lst_npz_names = []
    for np_name in glob.glob('*.np[yz]'):
        lst_npz_names.append(np_name)
    return lst_npz_names   
        
def get_ablated_areas_from_npz_files():
    tumour_depth_names_lst = []
    num_pixel_in_ablation_lst = []
    for npz in get_npz_filenames_lst():
        test_lst = np.load(npz)        
        for cem in test_lst:
            num_pixel_in_ablation_lst.append(np.sum(cem))   
            tumour_depth_names_lst.append(npz)
    return tumour_depth_names_lst, num_pixel_in_ablation_lst

def combine_cem_thresholds_into_dataframe():
    
    power = [0.6, 0.8, 1, 0.3, 0.5, 0.7, 0.1, 0.2, 0.3 ]
    d = {'power (W)': power, 'tumour depth': get_ablated_areas_from_npz_files()[0], 
                            'cem43>240 area': get_ablated_areas_from_npz_files()[1]}
    
    return pd.DataFrame(data=d)
    
         
def plot_area_of_ablation_vs_power():
    
    import seaborn as sns
    sns.set()
    
    
    df = combine_cem_thresholds_into_dataframe()
    
    df['tumour depth'] = df['tumour depth'].str.replace('_back_5mins_heat.npy', '')
    
    ax = sns.lineplot(x="power (W)", y="cem43>240 area", data=df, hue="tumour depth",
                        style='tumour depth', markers=True, dashes=False)
    plt.show()
# plot_area_of_ablation_vs_power()   
    
    
def plot_optimised_power_to_ablate_tumour_vs_diameter():
    
    tumour_diameter = [2, 4, 6]
    optimised_power = [0.21, 0.37, 0.78]
    
    fig, ax = plt.subplots()
    ax.plot(tumour_diameter, optimised_power, '-x')
    ax.set_xlabel('tumour diameter (mm)')
    ax.set_ylabel('optimised power (W) for tumour treatment')
plot_optimised_power_to_ablate_tumour_vs_diameter()
    
    
    
    
    
    
  
   