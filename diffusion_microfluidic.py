#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 11:08:12 2019

@author: charliejeynes

This script simulates diffusion of 40 kDa Dextran, using known parameters, and compares the simulation 
to timelapse microscopy data (in particular an image taken at 2 minutes after diffusion of the dextran from 
vascular channel into the tissue chamber)

"""

import numpy as np
import matplotlib.pyplot as plt
import skimage.io as skio
from skimage.color import rgb2gray
from skimage.transform import rescale, resize, downscale_local_mean


# this sets up the diffusion grid to correspond to pixels in the 'synvivo_chamber.bmp' image
w = 172 #this is in pixels
h = 281

# intervals in x-, y- directions, pixels
dx = dy = 1


# scale the pixels to microns ( 10 smaller, so divideby 10 to convert drug speed in um to pixels )
OCchannelwidthPixels = 18 # pixels
scalefactor = 200/OCchannelwidthPixels # the channel width in the microfludic is 200 um and 18 pixels in the image 

#  diffusivity of drug, microns^2.s-1 / pixel scalar factor (i.e. about 10)
D = 40 / scalefactor  # 40 microns^2.s-1 is diffusivity for FITC dextran, scalefactor is for direct comparison with the image data

Tcool, Thot = 0, 1

nx, ny = int(w/dx), int(h/dy)

dx2, dy2 = dx*dx, dy*dy
dt = dx2 * dy2 / (2 * D * (dx2 + dy2))

u0 = Tcool * np.ones((nx, ny))
u = np.empty((nx, ny))

#superimpose the image here 
im1 = skio.imread("synvivo_chamber.bmp")
imgray = rgb2gray(im1)
maskCh = imgray > 0.86
u0[maskCh ==True] = Thot # this creates the initial conditions for the simluation
#plt.imshow(u0[maskCh ==True])
plt.imshow(maskCh)


# this gets get the fluorescence in the tissue chamber bit for the real real
# Initial conditions - ring of inner radius r, width dr centred at (cx,cy) (mm)
circlemask = np.zeros((nx, ny))
r, cx, cy = 80, 105, 142
r2 = r**2
for i in range(nx):
    for j in range(ny):
        p2 = (i*dx-cx)**2 + (j*dy-cy)**2
        if p2 < r2:
            circlemask[i,j] = 1
imgray1 = imgray.copy()
circlemask1 = circlemask == 1 # makes it a mask 
imgray1[~circlemask1 == 1] = 0
plt.imshow(imgray1)  
#plt.imshow(imgray)    
#newimage =  np.empty((nx, ny)) 
#imgray[circlemask1 == True]
#plt.imshow(newimage)
# 
#plt.imshow(imgray(circlemask1 == 1))   

def do_timestep(u0, u):
    # Propagate with forward-difference in time, central-difference in space
    u[1:-1, 1:-1] = u0[1:-1, 1:-1] + D * dt * (
          (u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/dx2
          + (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/dy2 )

    u0 = u.copy()
    u0[maskCh ==True] = Thot # added this bit so the vascular channel 
    return u0, u

# Number of timesteps # 
nsteps = 2001
end = nsteps - 1
# Output 4 figures at these timesteps
mfig = [0, 1, 5, end]
fignum = 0
all_times = np.empty((nsteps, w, h))
fig = plt.figure()
for m in range(nsteps):
    u0, u = do_timestep(u0, u)
    all_times[m, : ,:] = u.copy()# figute this bit out!!!
    if m in mfig:
        fignum += 1
        print(m, fignum)
        ax = fig.add_subplot(220 + fignum)
        im = ax.imshow(u.copy(), cmap=plt.get_cmap('hot'), vmin=Tcool,vmax=Thot)
        ax.set_axis_off()
#        ax.set_title('{:.1f} ms'.format(m*dt*1000))
        ax.set_title('{:.1f} s'.format(m*dt))
#        ax.set_title('{:.1f} mins'.format(m*dt/60))
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
cbar_ax.set_xlabel('$T$ / K', labelpad=20)
fig.colorbar(im, cax=cbar_ax)
plt.show()

#compare the images 
fig, (ax1, ax2) = plt.subplots(1, 2)
#fig.suptitle('Horizontally stacked subplots')
ax1.imshow(imgray)
ax2.imshow(all_times[end, :, :])

# put a mask on the simulation image to deliniate the channels  
maskTC = imgray >0.4 #this is repeated above
plt.imshow(maskTC)
sim_image0 = all_times[end, :, :].copy()
maskTC[96:220, 65:218] = 1; maskTC[60:99, 75:200] = 1
sim_image0[~maskTC == 1] = 0
plt.imshow(sim_image0)
plt.colorbar()

#fig.subplots_adjust(right=0.85)
#cbar_ax = fig.add_axes([0.7, 0.15, 0.03, 0.5])
#im = ax.imshow(sim_image0, cmap=plt.get_cmap('hot'), vmin=Tcool,vmax=Thot)
#cbar_ax.set_xlabel('$T$ / K', labelpad=20)
#fig.colorbar(im, cax=cbar_ax)
#plt.show()

# this gets get the flurescnce in the tissue chamber for the simulated image
sim_image = all_times[end, :, :].copy()
sim_image[~circlemask1 == 1] = 0
plt.imshow(sim_image) 

#compare just the tissue chmaber bit
fig, (ax1, ax2) = plt.subplots(1, 2)
#fig.suptitle('Horizontally stacked subplots')
ax1.imshow(imgray1)
ax2.imshow(sim_image)

#binarised both the images to compare 
imgray1_bw = (imgray1 > 0.2) *1
sim_image_bw = (sim_image > 0.1) *1

#overaly the binarised images summing  
overlay_bw_im =  imgray1_bw + sim_image_bw
plt.imshow(overlay_bw_im)
plt.colorbar

#simplyt subtract 1 image from the other 
imgray2 = imgray1 - 0.1
substracted_im = imgray2 - sim_image
substracted_im[substracted_im < 0] = 0

# do structural similarity index on images
# mssim, S = ssim(imgray2, sim_image, full=True)

#THIS IS FIGURE that is in the MRC-CDA grant
#compare the exp with model with the good mask on it
#fig.suptitle('Horizontally stacked subplots')  
imgray_norm  = (imgray - np.min(imgray)) / (np.max(imgray) - np.min(imgray))

fig, (ax1, ax2) = plt.subplots(figsize=(13, 3), ncols=2)
exp = ax1.imshow(imgray_norm, extent=[0,2.81,0,1.72]) 
ax1.set_title("a. experimental", fontsize=20)
ax1.set_xlabel("mm")
ax1.set_ylabel("mm")
fig.colorbar(exp, ax=ax1)


np.savetxt('imgray_norm.csv', imgray_norm, delimiter=',')
np.savetxt('sim_image0.csv', sim_image0, delimiter=',')


sim = ax2.imshow(sim_image0, extent=[0,2.81,0,1.72])
ax2.set_title("b. simulation", fontsize=20)
ax2.set_xlabel("mm")
ax2.set_ylabel("mm")
fig.colorbar(sim, ax=ax2)

# diff = ax3.imshow(substracted_im, extent=[0,2.81,0,1.72], cmap='hot', vmin=Tcool, vmax=Thot,
#                             interpolation='none')

# # SSI = ax3.imshow(S, extent=[0,2.81,0,1.72], cmap='hot', vmin=Tcool, vmax=Thot,
# #                              interpolation='none')
# ax3.set_xlabel("mm")
# ax3.set_ylabel("mm")
# ax3.set_title("c. difference" ,fontsize=20)
# fig.colorbar(diff, ax=ax3)

plt.show()




#compute the dice score to compare model with experiment
#Y = pdist(X, 'dice')
def dice_coefficient(contour1, contour2):
    '''
    Calculate the dice score given 2 consensus metrics
    :param contour1: the first contour - citizen
    :param contour2: the second contour - the expert
    :return: The dice score as a value between 0 and 1 (float)
    '''
    return np.sum(contour1[contour2==1])*2.0 / (np.sum(contour1) + np.sum(contour2))

dice  = dice_coefficient(imgray1_bw, sim_image_bw)

# calculate the mean squared error
mse0 = np.square(np.subtract(imgray1[circlemask1 == 1], sim_image[circlemask1 == 1])).mean()

a = sim_image[circlemask1 == 1]
b = imgray1[circlemask1 == 1]
