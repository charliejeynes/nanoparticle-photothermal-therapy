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

filepath = '/Users/charliejeynes/Projects/dia/sim_data/absorption_dens_4mm_GNR.nc'
# source = '/Users/charliejeynes/Projects/dia/dia/output/mcrt/hits.nc'

def read_nc_file(filepath):
    ''' 
    reads in data from .nc file output from ARCTORUS and returns a 2D cross 
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


 