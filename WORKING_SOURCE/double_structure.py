#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 08:31:26 2025

@author: coletamburri
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import skimage
import scipy
import tol_colors as tc
import matplotlib.pylab as pl

from scipy.signal import convolve2d
from scipy.signal import convolve
from scipy.ndimage import gaussian_filter

# Function definitions for Gaussian fitting
def Gauss_func(x,A,mu,sigma,m,b):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))+ m*x + b

def double_gaussian( x, c1, mu1, sigma1, c2, mu2, sigma2 ,m,b):
    res =   (c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) )) \
          + (c2 * np.exp( - (x - mu2)**2.0 / (2.0 * sigma2**2.0) )) \
          + (m * x + b)
    return res


directory = '/Users/coletamburri/Desktop/double_loop_frame0_pre_destretch/'
time = '2024-08-08T20:12:32.333333'

blur = 0
#Determine mu
d = 151.68e9 # distance to the sun on 8 August source: https://theskylive.com/planetarium?objects=sun-moon-mercury-venus-mars-jupiter-saturn-uranus-neptune-pluto&localdata=40.01499%7C-105.27055%7CBoulder%20CO%20(US)%7CAmerica%2FDenver%7C0&obj=sun&h=14&m=01&date=2024-08-08#ra|9.242130505796545|dec|15.985314118209057|fov|50
solrad = 695700000

# Coordinates from DKIST are not correct, but define them anyways as a starting
# point.  Will co-align later in routine.

hpc1_arcsec = -175
hpc2_arcsec= 375

# image center
x_center = d*np.cos(hpc1_arcsec/206265)*np.sin(hpc2_arcsec/206265) # m
y_center = d*np.sin(hpc1_arcsec/206265) # m
z = solrad - d*np.cos(hpc1_arcsec/206265)*np.cos(hpc2_arcsec/206265) # m

# to mu value
rho = np.sqrt(x_center**2+y_center**2)
mu = np.sqrt(1-(rho/solrad)**2)

# Constants
spatial_samp = 0.017 # for vbi red at 656nm
arcsec_to_km = 727 # approximate arcsec to km conversion

# Arrays for coordinates of start and end
startx = []
starty = []

endx = []
endy = []

l = 0 # initialization of cut number

# for npz loading
path = '/Users/coletamburri/Desktop/VBI_Destretching/'
folder_vbi = 'AXXJL/' # 8 August X-class flare decay phase
filename = 'AXXJLselection_predestretch.npz'
array = np.load(path+folder_vbi+filename)['first50'] #first50 or brightening
numareas = 1

#frame to work with
frame = array[0,:,:]

halpha_samp = 0.017 #arcsec/pixel
resolution_aia_var = .6**2 #arcsec spatial sampling of sdo/aia
resolution_trace_var = .5**2 #spatial sampling of rhessi/trace
resolution_vbi_var = 0.017**2 #spatial sampling of DKIST/VBI in H-alpha filter
bbsogst_var = 0.034**2 #spatial sampling of BBSO/GST at H-alpha


pixels_psf_sig = round((np.sqrt(bbsogst_var-resolution_vbi_var))/halpha_samp)
convolved = gaussian_filter(np.asarray(frame),pixels_psf_sig)
frameblur = convolved



# X and Y coordinates of frame
xarr = np.arange(np.shape(frame)[0])
yarr = np.arange(np.shape(frame)[1])

# X and Y coordinates, in KM
xarr_km = xarr*spatial_samp
yarr_km = yarr*spatial_samp

# Meshgrid for plotting
XKM,YKM =np.meshgrid(xarr_km,yarr_km)

# Plot first frame
fig,ax=plt.subplots(dpi=400,figsize=(10,10))
ax.pcolormesh(frame,cmap='grey')
ax.set_aspect('equal')
ax.invert_yaxis()


cc = plt.ginput(numareas*2,timeout = 120)

ylo, yhi, xlo, xhi = int(cc[0][0]), int(cc[1][0]), int(cc[0][1]),\
    int(cc[1][1])

smallframe = frame[xlo:xhi,ylo:yhi]

smallframeblur = frameblur[xlo:xhi,ylo:yhi]

ncuts=31
lenfr = np.shape(smallframe)[0]
inds = np.linspace(160,190,ncuts)

colors = plt.cm.berlin(np.linspace(0,1,ncuts))

fig,ax=plt.subplots(1,2,dpi=200)

def normalize(data):
    normarr=(data-np.nanmin(data))/(np.nanmax(data)-np.nanmin(data)) 
    return normarr

for i in range(len(inds)):
    ind = round(inds[i])
    if ind > 165 and ind < 183:
        ax[0].plot(smallframe[ind,63:86]+4500*i,'-',color=colors[i],linewidth=1.5)
        ax[1].plot(smallframeblur[ind,63:86]+4500*i,'-',color=colors[i],linewidth=1.5)
    else:
        ax[0].plot(smallframe[ind,63:86]+4500*i,'--',color=colors[i],linewidth=1)
        ax[1].plot(smallframeblur[ind,63:86]+4500*i,'--',color=colors[i],linewidth=1)
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    
fig,ax=plt.subplots(1,2,dpi=200)
ax[0].imshow(smallframe[130:220,60:90],cmap='grey')
ax[0].set_xticks([])
ax[0].set_yticks([])

ax[1].imshow(smallframeblur[130:220,60:90],cmap='grey')
ax[1].set_xticks([])
ax[1].set_yticks([])

for i in np.arange(0,len(inds),2):
    ind = round(inds[i])
    if ind > 165 and ind < 183:
        ax[0].axhline(ind-130,linestyle='-',c=colors[i],linewidth=1)
        ax[0].axhline(ind-130,linestyle='-',c=colors[i],linewidth=1)
        ax[0].axhline(ind-130,linestyle='-',c=colors[i],linewidth=1)
        ax[1].axhline(ind-130,linestyle='-',c=colors[i],linewidth=1)
        ax[1].axhline(ind-130,linestyle='-',c=colors[i],linewidth=1)
        ax[1].axhline(ind-130,linestyle='-',c=colors[i],linewidth=1)
    else:
        ax[0].axhline(ind-130,linestyle='--',c=colors[i],linewidth=.5)
        ax[0].axhline(ind-130,linestyle='--',c=colors[i],linewidth=.5)
        ax[0].axhline(ind-130,linestyle='--',c=colors[i],linewidth=.5)
        ax[1].axhline(ind-130,linestyle='--',c=colors[i],linewidth=.5)
        ax[1].axhline(ind-130,linestyle='--',c=colors[i],linewidth=.5)
        ax[1].axhline(ind-130,linestyle='--',c=colors[i],linewidth=.5)














