#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

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
end=2

# package initialize
import dkistpkg_ct as DKISTanalysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update(matplotlib.rcParamsDefault)

# color scheme for plotting
muted = DKISTanalysis.color_muted2()

# path and file ID for ViSP data
path = '/Volumes/ViSP_External/pid_2_11/'
path2 = '/Volumes/ViSP_External/pid_2_11/'
folder1 = 'ANQZR'
folder2 = 'APQVP'  #need to add data for QS for 11 august...

# list of files in directory for DKIST/ViSP
dir_list2 = DKISTanalysis.pathdef(path,folder1) #flaretime
dir_list3 = DKISTanalysis.pathdef(path2,folder2) #qs

hbeta_low = 480
hbeta_high = 655


# Stonyhurst lon/lat position of the AR from JHv
lon = 55 #degrees
lat = -13 #degrees

lonqs = 0
latqs = 0


wl = 486.135 # central wavelength, H-beta

# spatial coordi nates
hpc1_arcsec, hpc2_arcsec, x_center, y_center, z, rho, mu, doppshrel,\
    doppshnonrel = \
    DKISTanalysis.spatialinit(path,folder1,dir_list2,lon,lat,wl)
    
# spatial coordinates, qs    
hpc1_arcsecqs, hpc2_arcsecqs, x_centerqs, y_centerqs, zqs, rhoqs, muqs, doppshrelqs,\
    doppshnonrelqs = \
    DKISTanalysis.spatialinit(path2,folder2,dir_list3,lon,lat,wl)

# get limb darkening coefficient 
clv_corr = DKISTanalysis.limbdarkening(wl, mu=mu, nm=True)
    # for Ca II H (require mu value for determination, be sure to specify
    # correct wl units)
    
# limb darkening coefficient, qs
clv_corrqs = DKISTanalysis.limbdarkening(wl, mu=muqs, nm=True)
    # for Ca II H (require mu value for determination, be sure to specify
    # correct wl units)
    
# time step start for chosen QS observations
startstepqs = 0
endstepqs=100
# # for earliest part of flare
# startstep=0#where does interesting bit begin?
# endstep=1400#where does interesting bit end?

# for later in flare
startstep=10000
endstep = 11400
print('here2')
# process multi-step raster - for qs time
image_data_arr_arr, rasterpos, times = \
    DKISTanalysis.multistepprocess(path,folder1,dir_list2,startstep=startstep,div=1,endstep=endstep)
    
# process multi-step raster - for qs time
image_data_arr_arr_qs, rasterpos_qs, times_qs = \
    DKISTanalysis.multistepprocess(path2,folder2,dir_list3,div=1,
                                   startstep=startstepqs,endstep=endstepqs)
    

#odd shift in the data for eid_2_11 - added logic to adjust for    
shift = .101

# spatial and dispersion axes for single observation (single slit step)
spatial_range, dispersion_range = DKISTanalysis.spatialaxis(path,folder1,
                                                            dir_list2,line='Hbeta',
                                                            pid='2_11',shift=shift)

# # old code, when basing QS on 15 August 2022 disk-center observations
# #only for 19 August observations, really - the QS will be different for others
# nonflare_average = np.load('/Users/coletamburri/Desktop/'+\
#                            'DKIST_Data/bolow_nonflare_average.npy')
# nonflare_stdevs = np.load('/Users/coletamburri/Desktop/'+\
#                           'DKIST_Data/bolow_nonflare_stdevs.npy')
# nonflare_fitvals = np.load('/Users/coletamburri/Desktop/'+\
#                            'DKIST_Data/bolow_nonflare_fit_vals.npy')
# nonflare_multfact = np.load('/Users/coletamburri/Desktop/'+\
#                             'DKIST_Data/bolow_nonflare_mult_fact.npy')
    
# Begin calibration based on QS

# Load Kurucz FTS Atlas
wlsel, ilamsel = DKISTanalysis.load_fts(dispersion_range)
wlsel=wlsel/10

# Average the QS data for space and time, for selected ranges
space_and_time_averaged_qs = \
    DKISTanalysis.comp_fts_to_qs(wlsel,ilamsel,dispersion_range, 
                                 image_data_arr_arr_qs)

# telluric lines for comparison (or other absorption lines if telluric not 
# available, as is the case for the Ca II H window).  Most of the next steps
# are not used for pid_1_38, but include in any case to model the use of lines
line1 = 485.974 #fe I
line2 = 486.261 # VI

# Definition of "telluric" or other absorption lines for calibration.
# The QS spectrum is already averaged in space and time, so any absorption lines
# (as is necessary for the 19 August 2022 observations) are ok

# indices of telluric lines in spectrum - lower
lowinds = [396,666]

# indices of telluric lines in spectrum - upper
highinds = [443,691]

# define the multiplication factor (polynomial), new dispersion range, fit values
# to scale the quiet sun to FTS atlas
cont_mult_facts,fit_vals,new_dispersion_range,dispersion_range_fin,rat=\
    DKISTanalysis.get_calibration_poly(dispersion_range,
                                       space_and_time_averaged_qs,wlsel,ilamsel,
                                       DKISTanalysis.find_nearest,line1,line2,
                                       lowinds,highinds,limbdark_fact=clv_corrqs,
                                       noqs_flag=2,cont_vals = [485.64,485.76,485.86,486.26,486.30,
                                                486.38,3486.424])

# calibrate the quiet sun 
calibrated_qs=fit_vals*space_and_time_averaged_qs/clv_corrqs
nonflare_average_avg = calibrated_qs
nonflare_multfact = fit_vals

# full width of half max of PSF to convolve with atlas to match instrument
fwhm = 0.05

# number of points to interpolate Atlas to in PSF convolve to match instrument
ntw = 45

# perform PSF convolution.  Result 'yconv' is Atlas*PSF
yconv=DKISTanalysis.psf_adjust(wlsel,ilamsel,fwhm,new_dispersion_range,
                               calibrated_qs,clv_corrqs,ntw,
                               DKISTanalysis.gaussian_psf)

#show comparison of atlas to qs
fig,ax=plt.subplots()
ax.plot(dispersion_range_fin,calibrated_qs,label='visp')
ax.plot(dispersion_range_fin,yconv*clv_corrqs,label='convolved')
ax.plot(wlsel,clv_corrqs*ilamsel,label='raw')
#ax.set_xlim([396.6,397.2]);ax.set_ylim([0,0.6e6])
ax.legend();plt.show()

#do another iteration of the calibration step after peforming the PSF conv.
cont_mult_facts,fit_vals,\
    new_dispersion_range,dispersion_range_fin,rat=\
        DKISTanalysis.get_calibration_poly(dispersion_range,
                                           space_and_time_averaged_qs,
                                           new_dispersion_range,yconv,
                                           DKISTanalysis.find_nearest,
                                           line1,line2,lowinds,highinds,
                                           limbdark_fact=clv_corrqs,noqs_flag=2)
        
# calibrated quiet sun, again, using updated fit values
calibrated_qs=fit_vals*space_and_time_averaged_qs/clv_corrqs
nonflare_average_avg = calibrated_qs
nonflare_multfact = fit_vals

# intensity calibration, background subtraction for flare-time                            
scaled_flare_time, bkgd_subtract_flaretime = \
    DKISTanalysis.scaling(image_data_arr_arr, nonflare_multfact,clv_corr,
                          nonflare_average_avg,end=end)


#indices to search within to find flare kernel center
spacelow = 0
spacehigh = -1

#indices of max intensity in each frame
maxindices = DKISTanalysis.maxintind(dispersion_range,bkgd_subtract_flaretime,
                                     hbeta_low,hbeta_high,
                                     spacelow,spacehigh)

### remove if desire tracking of central position ###
#testing with just the initial brightest point in the ribbon
# firstmaxonly = []
# for i in range(len(maxindices)):
#     firstmaxonly.append(maxindices[0])
# maxindices = firstmaxonly
####

### remove if desire tracking of central position ###
#testing with edge of flare ribbon
# edgeind = []
# for i in range(len(maxindices)):
#     edgeind.append(maxindices[i]+40)
# maxindices = edgeind
####


# plot intensity calibrated, background-subtracted spectra
DKISTanalysis.pltsubtract(dispersion_range_fin,nonflare_average_avg,
                          scaled_flare_time,muted,maxindices,end=end,pid='pid_2_11')

# variation in intensity value corresponding to wavelengths; PTE; to test
# for variations in pseudo-continuum.  If PTE high, cannot be explained by the
# solar activity difference between quiet sun and flare-time

########
## Comment out this block if don't need PTE (shouldn't, usually)
# stdevs_flaretime, ptes_flaretime = \
    # DKISTanalysis.deviations(bkgd_subtract_flaretime,nonflare_average,
                             # nonflare_stdevs)

# DKISTanalysis.pltptes(ptes_flaretime,image_data_arr_arr_raster1)

########
