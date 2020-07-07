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
import pandas as pd
import netCDF4
import sys
import os
import scipy.io as spio



filepath = '/Users/charliejeynes/Projects/dia/sim_data/move_back_tumour/power1W/'

#filepath = '/Users/charliejeynes/Projects/git_NP_PTT/'
 


def list_nc_files(filepath):
    '''
    lists the file path for all .nc simulation files in the subdirectory
    '''
    files_list = [os.path.join(filepath, file) for file in os.listdir(filepath)]
    return files_list



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


def read_in_heat_simulation(data):
    
    return spio.loadmat(data, squeeze_me=True)



def get_data_from_spio_dic(data, name_of_data):
    
    ''' this passes in the 'data' which is a 3d matrix from matlab, and the name 
    of the data as a string , so that it can look it up in the dict
    
    '''

    return data[name_of_data]



def convert_3Dmatrix_to_lst(heat_sims_matrix):
    '''
    
    takes a 3D matlab array and makes the 3rd dimension a list in an n*m array. 
    This is then in the right format for the rest of the functions in this
    scripts

    '''
    
    lst_of_heat_sims = []
    for third_dimension in range(heat_sims_matrix.shape[-1]):
        
        #profile = np.mean(cem43[:, 90:110, i], axis=0)
        heat_sim = heat_sims_matrix[:, :, third_dimension]
        lst_of_heat_sims.append(heat_sim)
            
    return lst_of_heat_sims   
    


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
        profile = section[130, :]
        line_profiles.append(profile)
        plt.plot(profile)
            
    return line_profiles   
 
 



def plot_image_profile(images, profiles):
    
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
            else:
                ax[r, c].semilogy(profiles[r])
                ax[r, c].axis(ymin=1,ymax=3e5)
                
                
def plot_hirsch_data():
    
    depth_from_skin = [0,  1,    2,  3,   4,   5,   6] # in mm 
    control_1min =    [5,  5,    4,  3,  2.5,  2,   1] # heat rise in C
    control_6min =    [12, 11.5, 9,  9,  8.0,  7.5, 6] # heat rise in C


    plt.plot(depth_from_skin, control_6min, scaley=False)

def plot_arc_nc_files():
    files_list = list_nc_files(filepath)
    cross_sections = get_sections_as_lst(files_list)
    log_cross_sections = log_data(cross_sections)   
    image_2D_sections(cross_sections, log_cross_sections)
    line_profiles = get_line_profile(cross_sections) 
    plot_image_profile(log_cross_sections, line_profiles)


#plot_arc_nc_files()


filepath = get_full_file_path(filepath, 'cem43.mat')
cem43 = read_in_heat_simulation(filepath)
mat_cem43 = get_data_from_spio_dic(cem43, 'cem43')
lst_cem43 = convert_3Dmatrix_to_lst(mat_cem43)
log_cem43 = log_data(lst_cem43)
#image_2D_sections(lst_cem43, log_cem43)
line_profiles = get_line_profile(lst_cem43) 
plot_image_profile(log_cem43, line_profiles)


plot_hirsch_data()






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




 