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
 

filepath = '/Users/charliejeynes/Projects/dia/sim_data/'

def list_nc_files(filepath):
    '''
    lists the file path for all .nc simulation files in the subdirectory
    '''
    files_list = [os.path.join(filepath, file) for file in os.listdir(filepath)]
    return files_list
files_list = list_nc_files(filepath)


def read_multiple_nc(files_list):
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
    cross_sections = []
    for file in files_list:
        nc = netCDF4.Dataset(file)
        print(nc.variables.keys()) # get all variable names
        data = nc.variables['data']  # access variable
        print(data) 
        section = data[:,100, :]
        section = section.data
        cross_sections.append(section)
        nc.close()
    return  cross_sections
cross_sections = read_multiple_nc(files_list)


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
log_cross_sections = log_data(cross_sections)   


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
image_2D_sections(cross_sections, log_cross_sections)

def get_line_profile(cross_sections):
    line_profiles = []
    for section in cross_sections:
        profile = np.mean(section[90:110, :], axis=0)
        line_profiles.append(profile)
        plt.plot(profile)
    return line_profiles
line_profiles = get_line_profile(cross_sections)    

def plot_image_profile(images, profiles):
    fig, ax = plt.subplots(3, 2, figsize=(20,13))
    fig.tight_layout()
    # axes are in a two-dimensional array, indexed by [row, col]
    for r in range(3):
        for c in range(2):
            if c%2 == 0:
                img = ax[r, c].imshow(images[r], cmap='jet', vmin=1, vmax=7)
                fig.colorbar(img, ax=ax[r,c])
            else:
                ax[r, c].semilogy(profiles[r])
                ax[r, c].axis(ymin=1,ymax=1e7)
plot_image_profile(log_cross_sections, line_profiles)

def main():
    '''
    read in file
    extract the datacube as an 3D matrix
    

    Returns
    -------
    None.

    '''











#filepath = '/Users/charliejeynes/Projects/dia/sim_data/absorption_dens_4mm_GNR.nc'
# source = '/Users/charliejeynes/Projects/dia/dia/output/mcrt/hits.nc'

def read_nc_file(filepath):
    ''' 
    reads in SINGLE data from .nc file output from ARCTORUS and returns a 2D cross 
    section of the datacube
    '''
    nc = netCDF4.Dataset(filepath)
    print(nc.variables.keys()) # get all variable names
    data = nc.variables['data']  # access variable
    print(data) 
    section = data[:,100, :]
    section = section.data
    nc.close()
    return  section

section = read_nc_file(filepath)


sectionLog = np.log10(section) # log10 transform data
def image_2D_section(section, sectionlog):
# view the native data and log10 data
    fig, (ax0, ax1) = plt.subplots(figsize=(13,5), ncols=2)
    img = ax0.imshow(section, cmap='jet', vmin=1e5, vmax=1e7)
    fig.colorbar(img, ax=ax0)
    img = ax1.imshow(sectionLog,cmap='jet') # logged image
    fig.colorbar(img, ax=ax1)
    plt.show()
image_2D_section(section, sectionLog)




 