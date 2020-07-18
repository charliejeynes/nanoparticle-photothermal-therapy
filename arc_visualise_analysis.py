#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 10:24:19 2020

@author: charliejeynes

This gets data from ARCTORUS .nc output files and visualises and analyses the 
results 


"""

import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
import netCDF4
import sys
import os
import scipy.io as spio
import cv2 as cv2
from skimage.measure import label, regionprops


plt.close('all')

#filepath = '/Users/charliejeynes/Projects/dia/sim_data/move_back_tumour/power1W/'

#filepath = '/Users/charliejeynes/Projects/git_NP_PTT/'

#filepath = '/Users/charliejeynes/Projects/git_NP_PTT/six_minutes.mat'

#filepath = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power0.28/'

# arc_nc_filepath = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power0.28/'

arc_nc_filepath = '/Users/charliejeynes/Projects/dia/sim_data/spot6mm_power_changing/'

matlab_files_path = '/Users/charliejeynes/Projects/dia/sim_data/sim_matlab_files/'
a = os.listdir(arc_nc_filepath) 


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
    datacube = datacube[:,:, :] 
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
resolution_in_mm_vector()

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
    
plot_hirsch_data()

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
power_in_spot = covert_spotsize_to_power(spot_radius = 0.15, power_quoted = 1)
    
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
plot_arc_nc_files(arc_nc_filepath)

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




def plot_absorbDensity_and_cem_profiles(arc_nc_filepath, cem_xy_boundaries):
    
    [absorbDensity_cross_sections, 
     absorbDensity_line_profiles,
     cross_sections, 
     log_cross_sections] =  get_arc_nc_data(arc_nc_filepath)
    
    cem_data_lst, cem_rot_data, 
    cem_data_line_profiles, 
    log_data_line_profiles = get_kwave_sim_data(matlab_files_path, 'cem43.mat', 'cem43')
    
    fig, ax = plt.subplots(3, 2, figsize=(20,13))
    fig.tight_layout()
    x=resolution_in_mm_vector()
    # axes are in a two-dimensional array, indexed by [row, col]
    for r in range(3):
        for c in range(2):  
            if c%2 == 0:
                img = ax[r, c].imshow(log_cross_sections[r], cmap='jet', vmin=1, vmax=7)
                ax[r, c].tick_params(axis='both', which='major', labelsize=20)
                ax[r, c].tick_params(axis='both', which='minor', labelsize=8)
                if cem_xy_boundaries[r].size: # this checks if the list is empty and if so, ignores it
                    ax[r,c].plot(cem_xy_boundaries[r][:, 0], cem_xy_boundaries[r][:, 1], linewidth=5, color='blue')
                cbar = fig.colorbar(img, ax=ax[r,c])
                cbar.ax.set_ylabel('absorption density (W/m$^3$)', fontsize=20, labelpad= 25)
                cbar.ax.tick_params(labelsize=20)

            else: 
                ax[r, c].semilogy(resolution_in_mm_vector(), absorbDensity_line_profiles[r], linewidth=5, color='black')
                ax[r, c].semilogy(resolution_in_mm_vector(), cem_data_line_profiles[r], linewidth=5, color='blue')
                ax[r, c].axis(ymin=10,ymax=1e7)
                ax[r, c].tick_params(axis='both', which='major', labelsize=20)
                ax[r, c].set_xlabel('(mm)', fontsize=25)
                #ax[r, c].set_ylabel('absorption density (W/m$^3$)', fontsize=20, rotation=360)
plot_absorbDensity_and_cem_profiles(arc_nc_filepath, cem_xy_boundaries)  
          

#BELOW WAS A TEST FOR A SINGLE FILE
#filepath = '/Users/charliejeynes/Projects/dia/sim_data/absorption_dens_4mm_GNR.nc'
# source = '/Users/charliejeynes/Projects/dia/dia/output/mcrt/hits.nc'

# def read_nc_file(filepath):
#     ''' 
#     reads in SINGLE data from .nc file output from ARCTORUS and returns a 2D cross 
#     section of the datacube
#     '''
#     nc = netCDF4.Dataset(filepath)
#     print(nc.variables.keys()) # get all variable names
#     data = nc.variables['data']  # access variable
#     print(data) 
#     section = data[:,100, :]
#     section = section.data
#     nc.close()
#     return  section

# section = read_nc_file(filepath)


# sectionLog = np.log10(section) # log10 transform data
# def image_2D_section(section, sectionlog):
# # view the native data and log10 data
#     fig, (ax0, ax1) = plt.subplots(figsize=(13,5), ncols=2)
#     img = ax0.imshow(section, cmap='jet', vmin=1e5, vmax=1e7)
#     fig.colorbar(img, ax=ax0)
#     img = ax1.imshow(sectionLog,cmap='jet') # logged image
#     fig.colorbar(img, ax=ax1)
#     plt.show()
# image_2D_section(section, sectionLog)




 