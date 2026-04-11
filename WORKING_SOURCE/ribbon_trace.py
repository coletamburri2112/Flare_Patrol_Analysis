#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 11:04:40 2026

@author: coletamburri
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import cv2
from skimage.measure import profile_line
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from pathlib import Path
from datetime import datetime
import matplotlib.dates as mdates ## Import required library
import sys
from scipy.optimize import curve_fit
sys.path.append('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/WORKING_SOURCE/')
import dkistpkg_ct as DKISTanalysis
from scipy.ndimage import map_coordinates


# Gaussian + linear background
def gaussian_linear(x, A, mu, sigma, B, C):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2)) + B*x + C

# get times and Fried parameter for entire series
vbi_filenames=[]

path = '/Volumes/ViSP_External/pid_2_11_VBI/MBVIDS'
# Specify the directory
folder = Path(path)

# List all files and folders
for item in folder.iterdir():
    if item.is_file():
        vbi_filenames.append(item.name)

vbi_filenames.sort()
del vbi_filenames[0:3]
vbi_filenames

timesvbi=[]
friedvbi=[]

for i in range(len(vbi_filenames)-1):
    timesvbi.append(fits.open(path+'/'+vbi_filenames[i])[1].header['DATE-BEG'])
    if fits.open(path+'/'+vbi_filenames[i])[1].header['AO_LOCK']==True:
        friedvbi.append(fits.open(path+'/'+vbi_filenames[i])[1].header['ATMOS_R0']*100)    
    else:
        friedvbi.append(float(0))
        
#convert to datetime
datetimevbi=[]
for i in range(len(timesvbi)):
    date_str = timesvbi[i]
    # Format: Year-Month-Day Hour:Minute:Second
    datetimevbi.append(datetime.strptime(date_str, "%Y-%m-%dT%H:%M:%S.%f"))
    
#strip day
onlytimevbi = []
for i in range(len(timesvbi)):
    date_str = timesvbi[i]
    # Format: Year-Month-Day Hour:Minute:Second
    onlytimevbi.append(date_str[-15:-7])
    
onlytimevbi_samp = onlytimevbi[180:350]
    
file = '/Users/coletamburri/Desktop/DKIST_Code/VBI_Destretching/MBVIDS/postdestretch_dataCube_Halpha_C_class_impulsive_phase_Halpha_180_350.fits'

destretch = fits.open(file)

ch=45

image = destretch[0].data[ch,:,:]

npoints=20

fig,ax=plt.subplots(dpi=200)
ax.imshow(image,cmap='sdoaia304')
ax.set_xlim([1700,2400])
ax.set_ylim([2550,1000])

fig.show()

cc = plt.ginput(npoints,timeout = 120)

xs = []
ys = []

for i in range(npoints):
    xs.append(cc[i][0])
    ys.append(cc[i][1])

def intensity_along_polyline(image, points):
    """
    Sample image intensity along straight-line segments connecting points.

    Parameters:
        image (2D ndarray): grayscale image
        points (list of (x, y)): ordered points

    Returns:
        coords (ndarray): (N, 2) sampled (x, y) coordinates
        intensities (ndarray): sampled intensity values
    """
    points = np.array(points)
    
    all_coords = []

    for i in range(len(points) - 1):
        x0, y0 = points[i]
        x1, y1 = points[i + 1]

        # Number of samples based on segment length
        length = int(np.hypot(x1 - x0, y1 - y0))
        if length == 0:
            continue

        # Evenly spaced points along the segment
        t = np.linspace(0, 1, length)
        x = x0 + t * (x1 - x0)
        y = y0 + t * (y1 - y0)

        all_coords.append(np.vstack([x, y]).T)

    coords = np.vstack(all_coords)

    # map_coordinates expects (row, col) = (y, x)
    sample_coords = np.vstack([coords[:, 1], coords[:, 0]])

    intensities = map_coordinates(image, sample_coords, order=1, mode='nearest')

    return coords, intensities

intensities_all = []
coords_all=[]

for i in range(0,140,1):
    image = destretch[0].data[i,:,:]
    coords, intensities = intensity_along_polyline(image, cc)
    intensities_all.append(intensities)
    coords_all.append(coords)
    
fig,ax=plt.subplots(dpi=200,figsize=(10,20));
ax.pcolormesh(np.transpose(intensities_all),cmap='Reds');
ax.invert_yaxis()
fig.show()

#overplot fried parameter on on this map

#compute FFT

l=len(intensities)
n=len(intensities)
dx = l / n      # Spatial sampling interval (meters)
space_x = np.linspace(0, l, n, endpoint=False)

all_freq = []
all_fft = []

for i in range(len(intensities_all)):
    chintensities = intensities_all[i]
    # Compute the 1D FFT
    fft_result = np.fft.fft(chintensities/np.max(chintensities))
    
    # Get the frequencies for the result
    freqs = np.fft.fftfreq(n,d=dx)
    
    all_freq.append(freqs)
    
    all_fft.append(fft_result)
    
all_fftarr = np.array(all_fft)
all_freqarr = np.array(all_freq)

pix_to_arcsec = 0.017
arcsec_to_km=727
km_to_Mm = 0.001
    
fig,ax=plt.subplots(2,1,dpi=200,figsize=(10,20));
ax.flatten()[0].pcolormesh(np.arange(140),np.arange(l)*pix_to_arcsec*arcsec_to_km*km_to_Mm,np.transpose(intensities_all),cmap='Reds');
ax.flatten()[0].invert_yaxis()

ax.flatten()[0].set_ylabel('Distance along trace [Mm]')
ax.flatten()[0].set_xlabel('Timestep')
ax.flatten()[1].pcolormesh(np.arange(140),km_to_Mm*arcsec_to_km*pix_to_arcsec*1/all_freqarr[0,1:n//2],np.transpose(np.abs(all_fftarr)[:,1:n//2]),cmap='magma',vmax=22);
#ax.flatten()[1].invert_yaxis()
ax.flatten()[1].axhline(0.050) # loop size
ax.flatten()[1].axhline(0.4) # bead size
ax.flatten()[1].set_ylim([.034,2])

ax.flatten()[1].set_yscale('log')
ax.flatten()[1].set_ylabel('Spatial scale [Mm]')
ax.flatten()[1].set_xlabel('Timestep')


fig.show()

from scipy.signal import savgol_filter



fig,ax=plt.subplots();
for i in range(0,120,20):
    ax.plot(km_to_Mm*arcsec_to_km*pix_to_arcsec*1/all_freqarr[0,1:n//2],savgol_filter(np.abs(all_fftarr)[i,1:n//2],window_length=20,polyorder=3),label=str(i))
ax.set_xscale('log')
ax.set_xlim([0,2])
ax.legend()
fig.show()






























