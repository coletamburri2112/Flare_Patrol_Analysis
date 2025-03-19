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
path = '/Volumes/ViSP_External/pid_1_84/'
folder1 = 'AZVXV'
folder2 = 'ANYDJ' # for QS calibration - data from 22UT on 19 August

# list of files in directory for DKIST/ViSP
dir_list2 = DKISTanalysis.pathdef(path,folder1) #flaretime
dir_list3 = DKISTanalysis.pathdef(path,folder2) #qs

# Stonyhurst lon/lat position of the AR from JHv
lon = 58.57 #degrees
lat = -29.14 #degrees

# Stonyhurst lon/lat position of the AR from JHv for QS
lon = 60 #degrees
lat = -28 #degrees


wl = 396.847 # central wavelength, Ca II H

# spatial coordinates
hpc1_arcsec, hpc2_arcsec, x_center, y_center, z, rho, mu, doppshrel,\
    doppshnonrel = \
    DKISTanalysis.spatialinit(path,folder1,dir_list2,lon,lat,wl)
    
# spatial coordinates, qs    
hpc1_arcsecqs, hpc2_arcsecqs, x_centerqs, y_centerqs, zqs, rhoqs, muqs, doppshrelqs,\
    doppshnonrelqs = \
    DKISTanalysis.spatialinit(path,folder2,dir_list3,lon,lat,wl)

# get limb darkening coefficient 
clv_corr = DKISTanalysis.limbdarkening(wl, mu=mu, nm=True)
    # for Ca II H (require mu value for determination, be sure to specify
    # correct wl units)
    
# limb darkening coefficient, qs
clv_corrqs = DKISTanalysis.limbdarkening(wl, mu=muqs, nm=True)
    # for Ca II H (require mu value for determination, be sure to specify
    # correct wl units)
    
# time step start for chosen QS observations
startstepqs = 250

# process multi-step raster - for qs time
image_data_arr_arr_qs, rasterpos_qs, times_qs = \
    DKISTanalysis.multistepprocess(path,folder2,dir_list3,div=10,
                                   startstep=startstepqs)
    
# for pid_1_84 - process four-step raster
image_data_arr_arr,i_file_raster1, for_scale, times_raster1, times = \
    DKISTanalysis.fourstepprocess(path,folder1,dir_list2)
    
# spatial and dispersion axes for single observation (single slit step)
spatial_range, dispersion_range = DKISTanalysis.spatialaxis(path,folder1,
                                                            dir_list2,line='Ca II H')

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
                                 image_data_arr_arr_qs, lowint = -500,highint=-1,
                                 timelow=20,timehigh=46)

# telluric lines for comparison (or other absorption lines if telluric not 
# available, as is the case for the Ca II H window).  Most of the next steps
# are not used for pid_1_38, but include in any case to model the use of lines
line1 = 396.7432
line2 = 396.926

# Definition of "telluric" or other absorption lines for calibration.
# The QS spectrum is already averaged in space and time, so any absorption lines
# (as is necessary for the 19 August 2022 observations) are ok

# indices of telluric lines in spectrum - lower
lowinds = [400,625]

# indices of telluric lines in spectrum - upper
highinds = [450,690]

# define the multiplication factor (polynomial), new dispersion range, fit values
# to scale the quiet sun to FTS atlas
cont_mult_facts,fit_vals,new_dispersion_range,dispersion_range_fin,rat=\
    DKISTanalysis.get_calibration_poly(dispersion_range,
                                       space_and_time_averaged_qs,wlsel,ilamsel,
                                       DKISTanalysis.find_nearest,line1,line2,
                                       lowinds,highinds,limbdark_fact=clv_corrqs,
                                       noqs_flag=2)

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
ax.set_xlim([396.6,397.2]);ax.set_ylim([0,0.6e6])
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
    DKISTanalysis.scaling(for_scale, nonflare_multfact,clv_corr,
                          nonflare_average_avg,end=end)
  
# definition of index bounds for Ca II H lines
caII_low_foravg = 525
caII_high_foravg = 600

#indices to search within to find flare kernel center
spacelow = 1000
spacehigh = -1

scaled_flare_time_ar, bkgd_subtract_flaretime_ar = \
    DKISTanalysis.scaling(image_data_arr_arr, nonflare_multfact,clv_corr,
                          nonflare_average_avg,end=end)
    
maxindices_ar = DKISTanalysis.maxintind(dispersion_range,
                                                bkgd_subtract_flaretime_ar,
                                                caII_low_foravg,caII_high_foravg,
                                                spacelow,spacehigh)



#indices of max intensity in each frame
maxindices = DKISTanalysis.maxintind(dispersion_range,bkgd_subtract_flaretime,
                                     caII_low_foravg,caII_high_foravg,
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
                          scaled_flare_time,muted,maxindices,end=end)

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

# equivalent widths, effective widths, widths
# indices defining line locations
caII_low = 480
caII_high = 650
hep_low = 700
hep_high = 850

#first profile, to use in testing
sample_flaretime = bkgd_subtract_flaretime[0,:,maxindices[0]]

########
# perform following block only if need to calculate continuum window 
# independently of width determiation; 
# exists within width determination script as well

#nolines, cont_int_array, cont_int_wave_array = \
    # DKISTanalysis.contwind(sample_flaretime,dispersion_range,maxinds,avgs,
                           # low,high)
# line widths, strengths - initialize arrays
ew_CaII_all_fs_ar = np.zeros((len(scaled_flare_time_ar)-0,
                           np.shape(bkgd_subtract_flaretime_ar)[2]))
ew_hep_all_fs_ar = np.zeros((len(scaled_flare_time_ar)-0,
                          np.shape(bkgd_subtract_flaretime_ar)[2]))
eqw_CaII_all_fs_ar = np.zeros((len(scaled_flare_time_ar)-0,
                            np.shape(bkgd_subtract_flaretime_ar)[2]))
eqw_hep_all_fs_ar = np.zeros((len(scaled_flare_time_ar)-0,
                           np.shape(bkgd_subtract_flaretime_ar)[2]))
width_CaII_all_fs_ar = np.zeros((len(scaled_flare_time_ar)-0,
                              np.shape(bkgd_subtract_flaretime_ar)[2]))
width_hep_all_fs_ar = np.zeros((len(scaled_flare_time_ar)-0,
                             np.shape(bkgd_subtract_flaretime_ar)[2]))

int_CaII_all_fs_ar = np.zeros((len(scaled_flare_time_ar)-0,
                             np.shape(bkgd_subtract_flaretime_ar)[2]))
int_hep_all_fs_ar = np.zeros((len(scaled_flare_time_ar)-0,
                             np.shape(bkgd_subtract_flaretime_ar)[2]))


# line widths, strength determination
ew_CaII_all_fs_ar, ew_hep_all_fs_ar, eqw_CaII_all_fs,\
    eqw_hep_all_fs, width_CaII_all_fs, width_hep_all_fs = \
        DKISTanalysis.widths_strengths(ew_CaII_all_fs_ar,eqw_CaII_all_fs_ar,
                                       width_CaII_all_fs_ar,ew_hep_all_fs_ar,
                                       eqw_hep_all_fs_ar,width_hep_all_fs_ar,caII_low,
                                       caII_high,hep_low,hep_high,
                                       scaled_flare_time_ar,
                                       bkgd_subtract_flaretime_ar, 
                                       dispersion_range_fin)
        
ew_CaII_all_fs_arr, eqw_CaII_all_fs, width_CaII_all_fs_arr, int_CaII_all_fs_ar=\
    DKISTanalysis.widths_strengths_oneline(ew_CaII_all_fs_ar,eqw_CaII_all_fs_ar,
                                           width_CaII_all_fs_ar,int_CaII_all_fs_ar,
                                           caII_low,caII_high,scaled_flare_time_ar,
                                           bkgd_subtract_flaretime_ar,dispersion_range,
                                           maxindices_ar,alt=1)
    
ew_hep_all_fs_arr, eqw_hep_all_fs, width_hep_all_fs_arr, int_hep_all_fs_ar=\
    DKISTanalysis.widths_strengths_oneline(ew_hep_all_fs_ar,eqw_hep_all_fs_ar,
                                           width_hep_all_fs_ar,int_hep_all_fs_ar,
                                           hep_low,caII_high,scaled_flare_time_ar,
                                           bkgd_subtract_flaretime_ar,dispersion_range,
                                           maxindices_ar,alt=1)
    
DKISTanalysis.plt_line_characteristics(ew_CaII_all_fs_ar,eqw_CaII_all_fs_ar,
                                       width_CaII_all_fs_ar,int_CaII_all_fs_ar,
                                       ew_hep_all_fs_ar,eqw_hep_all_fs_ar,
                                       width_hep_all_fs_ar,int_hep_all_fs_ar,
                                       maxindices_ar,times,muted)
        