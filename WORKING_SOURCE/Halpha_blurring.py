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
data1 = np.load('/Users/coletamburri/Desktop/VBI_Destretching/AXXJL/brighteningframes.npz')

data = normalize_3d_array(data1['brightening'])

halpha_samp = 0.017 #arcsec/pixel
resolution_aia = .6 #arcsec spectral sampling of aia

pixels_psf_sig = round(resolution_aia/halpha_samp)


convolved= gaussian_filter(data[40,:,:],pixels_psf_sig)

fig,ax=plt.subplots();ax.imshow(convolved,cmap='grey');fig.show()


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
            result[i,:,:] = data[i,:,:] - data[i - 1,:,:]
    return result

running_diff = running_difference(data)

dataCube=np.asarray(running_diff)

#save pre-destretch data cube of the first 100 (sub-)frames
fits.writeto('/Users/coletamburri/Desktop/runningdifference.fits',dataCube,overwrite=True)

#restore datacube
dataCube=fits.open('/Users/coletamburri/Desktop/runningdifference.fits')[0].data

#make a movie 
utilvbi.storeSequence(dataCube,'/Users/coletamburri/Desktop/runningdifference.mp4', dpi=300, write=True)

Video('/Users/coletamburri/Desktop/runningdifference.mp4', embed=True, width=600, height=600)

