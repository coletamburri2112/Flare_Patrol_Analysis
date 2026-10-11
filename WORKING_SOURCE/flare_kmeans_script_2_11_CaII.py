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
k_values = np.arange(2,50,1)
k_int = 15

## Want to output metrics?
metricstest = 0

## Define mask cutoff (median intensity)
cutoff0=1.5

## KEYWORDS FOR RUN
full_scan_qs = 0 # if =1, subtract the pre-flare sun pixel-by-pixel
                 # if 0, subtract averaged "non-flare" from all pixels
adjust='no' # if 'no', sorts clusters by weighted mean (or other choice); otherwise
            # 'byhand' means that some cluster order is switched for clarity
nsteps = 91 # number of slit steps per ViSP scan
start = 57 # where does the interesting bit of the ViSP data start in loaded datacube?
n_init = 1 # number of times to initialize the k-means clustering; 
            # 10 by default (though 1 is probably ok if using k-means++ as initializer)
manyscan = 1 # if =0, only one scan of the ViSP; if =1, many
nframes = 10 # number of scans (if manyscan)
c = 299792458 # speed of light in m/s
pid = '2.11'

#for Ca II H
linelow = 480
linecent = 595
linehigh = 690

## LINE CENTER IN NANOMETERS
cent = 396.85

## DEFINE SAVED FILE
filename = '/Volumes/ViSP_External/CaII_11Aug2024_Cclass_newcalib.npz'
flare_file = np.load(filename)

## READ-IN SAVED CO-ALIGNED DKIST COORDINATION
#dkist_coord_file = '/Users/coletamburri/Desktop/DKIST_Flares/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/ViSPcoords_newcalib.npz'
dkist_coord_file = '/Users/coletamburri/Desktop/DKIST_August2024_Flare_RefResponse/ViSPcoords_newcalib.npz'
dkist_coords = np.load(dkist_coord_file)

## STORE SPATIAL COORDINATES
xarr_caII = dkist_coords['xarr_caII']
yarr_caII = dkist_coords['yarr_caII']

xarr_hbeta = dkist_coords['xarr_hbeta']
yarr_hbeta = dkist_coords['yarr_hbeta']

# datacube shape - [ROI X, wave, ROI Y]
flare_arr0 = flare_file['flare']
wave = flare_file['wl']
times = flare_file['time']

## DEFINE SPATIAL LIMITS TO INCLUDE IN MASKING
startspace = 300
endspace = 1700

## SELECT WHICH WAVELENGTHS FOR SPECTRAL AXIS
selwls = wave[linelow:linehigh]

# Clustering object
clusterer = SpectralClusterer(flare_arr0,pid,wave,times)

# Prep data - make sure to define frame-line-metric for the map!
clusterer.prep_data(cutoff=cutoff0,cent=cent, line_low=linelow,line_high=linehigh,\
                    frame_line_metric = 'mean', coreind=linecent,startspace=startspace,\
                        endspace=endspace, start=start,nsteps=nsteps,nframes=nframes)

# Perform the clustering
clusterer.clustering(k_values=k_values, metricstest=metricstest, n_init=n_init)

# Extract the output variables - this will be indexed based on the k value
results = clusterer.results

# Extract info and sort for this run
clusterer.sort(k_choice=k_int,method='WM')

# Plot results
clusterer.plot_maps(k_choice=k_int,startspace=startspace,endspace=endspace,
                    xarr=xarr_caII,yarr=yarr_caII,nframes=nframes,nsteps=nsteps)
