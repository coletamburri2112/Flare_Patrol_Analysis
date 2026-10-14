#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 06:18:46 2024

@author: coletamburri
"""

import numpy as np
from astropy.io import fits
from Clusterer import SpectralClusterer


## DEFINE NUMBER OF CLUSTERS (determined by inspection for now, but this should be updated)
k_values = [10,25]

## Want to output metrics?
metricstest = 0

## Define mask cutoff (median intensity)
cutoff0 = 1.1 # this is best option when taking mean across spectral window for mask
cutoff0 = 1.3 # testing options when taking core intensity for mask

## KEYWORDS FOR RUN
full_scan_qs = 0 # if =1, subtract the pre-flare sun pixel-by-pixel
                 # if 0, subtract averaged "non-flare" from all pixels
adjust='no' # if 'no', sorts clusters by weighted mean (or other choice); otherwise
            # 'byhand' means that some cluster order is switched for clarity
nsteps = 148 # number of slit steps per ViSP scan
start = 0 # where does the interesting bit of the ViSP data start in loaded datacube?
n_init = 10 # number of times to initialize the k-means clustering; 
            # 10 by default (though 1 is probably ok if using k-means++ as initializer)
manyscan = 1 # if =0, only one scan of the ViSP; if =1, many
nframes = 10 # number of scans (if manyscan)
c = 299792458 # speed of light in m/s

## SPECTRAL AXIS PIXELS FOR EACH LINE
linelow = 380
linecent = 445
linehigh = 520

# #for fe I
# linelow = 145
# linecent = 167
# linehigh = 185

## LINE CENTER IN NANOMETERS
cent = 854.21


## DEFINE SAVED FILE
#filename = '/Users/coletamburri/Desktop/Aaron_Jul2026_DKIST_dataset/window_1746_arm2_Fe_I_(630.25_nm)_aligned_roi.fits'
filename = '/Users/coletamburri/Desktop/Aaron_Jul2026_DKIST_dataset/window_1746_arm1_Ca_II_(854.21_nm)_aligned_roi.fits'

flare_file = fits.open(filename)

# datacube shape - [map repeat, Stokes param, wavelength, ROI x (slit direction), ROI y (scan direction)]
flare_arr0 = flare_file[0].data[:,0,:,:,:]

# spatial coords
yarr = flare_file[2].data[0,0,:]
xarr = flare_file[2].data[0,:,0]

# Wavelengths
wave = flare_file[1].data/10 # convert to nm

# Times for array - placeholder, will need to get these from Aaron
times = np.arange(np.shape(flare_arr0)[0]*np.shape(flare_arr0)[2])

## DEFINE SPATIAL LIMITS TO INCLUDE IN MASKING
startspace = 0
endspace = -1

## SELECT WHICH WAVELENGTHS FOR SPECTRAL AXIS
selwls = wave[linelow:linehigh]

# Clustering object
clusterer = SpectralClusterer(flare_arr0,pid='4.17',wave=wave,times=times)

# Prep data - make sure to define frame-line-metric for the map!
clusterer.prep_data(cutoff=cutoff0,cent=cent, line_low=linelow,line_high=linehigh,\
                    frame_line_metric = 'core', pid='4.17',\
                        nsteps=nsteps,nframes=nframes,coreind=linecent)
    
clusterer.clustering(k_values=k_values, metricstest=metricstest, n_init=n_init)

# Extract the output variables - this will be indexed based on the k value
results = clusterer.results

# decide this after looking at output of above:
k_int = 10

# a number that goes into k_int which makes things nice to plot
#... probably a way to automate this...
div = 5

# Extract info and sort for this run
clusterer.sort(k_choice=k_int,method='WM')


# Plot results
clusterer.plot_maps(k_choice=k_int,startspace=startspace,endspace=endspace,\
                    xarr=xarr,yarr=yarr,nframes=nframes,nsteps=nsteps, \
                        pid = '4.17',div=div)

# Plot clustering metrics
clusterer.plot_metrics(kstart=4,kend=50,klow=k_values[0])


