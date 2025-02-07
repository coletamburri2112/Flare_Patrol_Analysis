#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 14:39:00 2025

@author: coletamburri
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib
from datetime import date, datetime, timedelta, time

from scipy.signal import convolve2d
from scipy.signal import convolve
from scipy.ndimage import gaussian_filter
import utilvbi
import matplotlib.patches as patches

from skimage.exposure import match_histograms


def perdelta(start, end, delta):
    curr = start
    while curr <= end:
        yield curr
        curr += delta
        
def normalize_3d_array(arr):
    min_val = np.nanmin(arr)
    max_val = np.nanmax(arr)

    return (arr - min_val) / (max_val - min_val)

#make times
sttime = datetime(2024,8,8,20,12,32,333333)
endtime = datetime(2024,8,8,21,5,7,0)

# plt.rcParams['text.usetex']=True
# plt.rcParams['font.family']='sans-serif'
# plt.rcParams['font.sans-serif'] = ['Tahoma']
# plt.rcParams['axes.labelsize'] = 25
# plt.rcParams['lines.linewidth'] = 2
# matplotlib.rc('xtick', labelsize=20) 
# matplotlib.rc('ytick', labelsize=20) 


stack=[]
for result in perdelta(sttime , endtime, timedelta(seconds=2.666666)):
    stack.append(str(result))
    
timeshhmmss = []

# timesdt = []

# for i in stack:
#     timesdt.append(datetime.strptime(i, '%Y-%m-%dT%H:%M:%S.%f'))
    
# timesdthr = []

# for i in timesdt:
#     timesdthr.append(datetime.strftime(i, 
#                                  "%H:%M:%S"))

# times in correct format for plotting
for i in range(len(stack)):
    timeshhmmss.append(stack[i][-15:-7])
    
def gaussian_psf(x, fwhm):
	#x = wavelength [A]
	# fwhm in [A]
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))  # Calculate sigma from FWHM
    tr = np.exp(-(x)**2 / (2 * (sigma**2)))
    tr /= tr.sum()
    return tr
    
# load fits file

#this file is indices 150 to 250
data1 = np.load('/Users/coletamburri/Desktop/August_2024_DKIST_Flares/VBI_X_class/brighteningframes.npz')

data = data1['brightening']
# gbref = data[0]
# timegb,yy,xx=data.shape

# dataCubeTrackedhist=[]
# for idata in range(timegb):
#     image=data[idata,:,:]
#     matched = match_histograms(image, gbref)
#     dataCubeTrackedhist.append(matched)


# datahist = dataCubeTrackedhist

halpha_samp = 0.017 #arcsec/pixel
resolution_aia = .6 #arcsec spatial sampling of sdo/aia
resolution_trace = .5 #spatial sampling of rhessi/trace

pixels_psf_sig = round(resolution_trace/halpha_samp)

l=0
fig,ax=plt.subplots(3,3,dpi=300);
for i in np.arange(0,90,10):
    convolved= gaussian_filter(np.asarray(data[i]),pixels_psf_sig)
    
    
    ys=np.arange(1400,2300,1)
    xs=np.arange(2000,3200,1)
    
    X,Y = np.meshgrid(xs,ys)
    
    ys2=np.arange(1200,2500,1)
    xs2=np.arange(1800,3400,1)
    
    X2,Y2 = np.meshgrid(xs2,ys2)
    ax.flatten()[l].pcolormesh(X2,Y2,data[i,1200:2500,1800:3400],cmap='grey')
    ax.flatten()[l].contour(X,Y,convolved[1400:2300,2000:3200],\
               levels=np.linspace(Z.min(), Z.max(), 7), cmap='hot',width=7)
    rect = patches.Rectangle((2000,1400), 1200, 900, linewidth=3, edgecolor='k', facecolor='none')
    
    # Add the patch to the Axes
    ax.flatten()[l].add_patch(rect)
    ax.flatten()[l].set_ylim([2500,1200])
    ax.flatten()[l].set_xlim([1800,3400])
    ax.flatten()[l].set_axis_off()
    ax.flatten()[l].set_aspect('equal')
    l+=1
    
fig.show()

def running_difference(data):
    """
    Calculates the running difference of a list.

    Args:
        data: A list of numbers.

    Returns:
        A list containing the running difference.
    """
    result = np.zeros(np.shape(data))
    for i in range(len(data)):

        if i == 0:
            result[i,:,:]= np.zeros(np.shape(data)[1:]) # The first element has no previous element
        else:
            result[i,:,:] = data[i] - data[i - 1]
    return result

running_diff = running_difference(data)

hdu = fits.PrimaryHDU(running_diff)
hdu.writeto('runningdifference_Xclassbrightening.fits', overwrite=True) # overwrite=True will overwrite the file if it exists


