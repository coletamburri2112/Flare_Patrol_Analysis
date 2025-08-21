
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2 December 2023
Last updated: 7 March 2024
Author: Cole Tamburri, University of Colorado Boulder, National Solar 
Observatory, Laboratory for Atmospheric and Space Physics

Description of script: 
    Used to process and perform preparatory analysis on DKIST ViSP L1 data, 
    .fits file format.  Intensity calibration via quiet Sun also necessary as 
    input, performed via another script.  Here, observations are calibrated 
    using the factors determined through that external process, and adjustments
    made for limb darkening, etc., if necessary. Line widths and strengths are 
    calculated. Emission line modeling is performed in order to track flare 
    dynamics.

"""
# shift of wavelength range by inspection
end=0

# package initialize
import dkistpkg_ct as DKISTanalysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update(matplotlib.rcParamsDefault)

# color scheme for plotting
muted = DKISTanalysis.color_muted2()

# path and file ID for ViSP data
path = '/Users/coletamburri/Desktop/August_2022_Observations_Paper/'
folder1 = 'All_Stokes_pid_1_84'
#folder2 = 'ANYDJ' # for QS calibration - data from 22UT on 19 August
#folder2 = 'AJDZR' # alternative for QS calibration 

# list of files in directory for DKIST/ViSP
dir_list2 = DKISTanalysis.pathdef(path,folder1,stokesi=0) #flaretime
#dir_list3 = DKISTanalysis.pathdef(path,folder2) #qs

# Stonyhurst lon/lat position of the AR from JHv
lon = 59 #degrees
lat = -31 #degrees

# Stonyhurst lon/lat position of the AR from JHv for QS, ANYDJ
lonqs = 59 #degrees
latqs = -28 #degrees

# Stonyhurst lon/lat position of the AR from JHv for QS, BGL


wl = 396.847 # central wavelength, Ca II H

# spatial coordinates
hpc1_arcsec, hpc2_arcsec, x_center, y_center, z, rho, mu, doppshrel,\
    doppshnonrel = \
    DKISTanalysis.spatialinit(path,folder1,dir_list2,lon,lat,wl,flag=1)
    
# # spatial coordinates, qs    
# hpc1_arcsecqs, hpc2_arcsecqs, x_centerqs, y_centerqs, zqs, rhoqs, muqs, doppshrelqs,\
#     doppshnonrelqs = \
#     DKISTanalysis.spatialinit(path,folder2,dir_list3,lonqs,latqs,wl,flag=2)

# get limb darkening coefficient 
clv_corr = DKISTanalysis.limbdarkening(wl, mu=mu, nm=True)
    # for Ca II H (require mu value for determination, be sure to specify
    # correct wl units)
    
# # limb darkening coefficient, qs
# clv_corrqs = DKISTanalysis.limbdarkening(wl, mu=muqs, nm=True)
#     # for Ca II H (require mu value for determination, be sure to specify
#     # correct wl units)
    
# time step start for chosen QS observations
#startstepqs = 250 # for longer file
startstepqs = 0 # for shorter file


# # process multi-step raster - for qs time
# image_data_arr_arr_qs, rasterpos_qs, times_qs = \
#     DKISTanalysis.multistepprocess(path,folder2,dir_list3,div=10,
#                                    startstep=startstepqs)
    
# for pid_1_84 - process four-step raster
image_data_arr_arr,i_file_raster1, for_scale, times_raster1,times = \
    DKISTanalysis.fourstepprocess(path,folder1,dir_list2,fullstokes=1)