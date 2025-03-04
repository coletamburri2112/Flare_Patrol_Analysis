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
end=5

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
folder2 = 'APQVP'  #QS - but change to hbeta specific
folder3hbeta = 'ARYEE'

# list of files in directory for DKIST/ViSP
dir_list2hbeta = DKISTanalysis.pathdef(path,folder3hbeta) #flaretime
dir_list3 = DKISTanalysis.pathdef(path2,folder2) #qs

# Stonyhurst lon/lat position of the AR from JHv
lon = 57.99 #degrees
lat = -15 #degrees

lonqs = 0
latqs = 0



wlhbeta = 486.2 # central wavelength, Hbeta


    
# spatial coordinates hbeta
hpc1_arcsechbeta, hpc2_arcsechbeta, x_centerhbeta, y_centerhbeta, zhbeta, rhohbeta, muhbeta, doppshrelhbeta,\
    doppshnonrelhbeta = \
    DKISTanalysis.spatialinit(path,folder3hbeta,dir_list2hbeta,lon,lat,wlhbeta)
    
# spatial coordinates, qs    
hpc1_arcsecqs, hpc2_arcsecqs, x_centerqs, y_centerqs, zqs, rhoqs, muqs, doppshrelqs,\
    doppshnonrelqs = \
    DKISTanalysis.spatialinit(path2,folder2,dir_list3,lon,lat,wlhbeta)


    
# limb darkening coefficient, qs
clv_corrqs = DKISTanalysis.limbdarkening(wlhbeta, mu=muqs, nm=True)
    # for H-beta (require mu value for determination, be sure to specify
    # correct wl units)
    
clv_corrhbeta = DKISTanalysis.limbdarkening(wlhbeta, mu=muhbeta, nm=True)
    # for H-beta (require mu value for determination, be sure to specify
    # correct wl units)
# time step start for chosen QS observations
startstepqs = 0
endstepqs=100
startstep=2000 #where does interesting bit begin?
endstep=2500#where does interesting bit end?

    
# process multi-step raster - for hbeta
image_data_arr_arrhbeta, rasterposhbeta, timeshbeta = \
    DKISTanalysis.multistepprocess(path,folder3hbeta,dir_list2hbeta,startstep=startstep,div=1,endstep=endstep)

#array processed
# process multi-step raster - for qs time
image_data_arr_arr_qs, rasterpos_qs, times_qs = \
    DKISTanalysis.multistepprocess(path2,folder2,dir_list3,div=1,
                                   startstep=startstepqs,endstep=endstepqs)
    
    
    

# # spatial and dispersion axes for single observation (single slit step)

spatial_range , dispersion_range = DKISTanalysis.spatialaxis(path,folder3hbeta,
                                                            dir_list2hbeta,line='H-beta',
                                                            pid='2_11')

# # # old code, when basing QS on 15 August 2022 disk-center observations
# # #only for 19 August observations, really - the QS will be different for others
# # nonflare_average = np.load('/Users/coletamburri/Desktop/'+\
# #                            'DKIST_Data/bolow_nonflare_average.npy')
# # nonflare_stdevs = np.load('/Users/coletamburri/Desktop/'+\
# #                           'DKIST_Data/bolow_nonflare_stdevs.npy')
# # nonflare_fitvals = np.load('/Users/coletamburri/Desktop/'+\
# #                            'DKIST_Data/bolow_nonflare_fit_vals.npy')
# # nonflare_multfact = np.load('/Users/coletamburri/Desktop/'+\
# #                             'DKIST_Data/bolow_nonflare_mult_fact.npy')
    
# Begin calibration based on QS

# Load Kurucz FTS Atlas
#wlsel, ilamsel = DKISTanalysis.load_fts(dispersion_range-.01)
#wlsel=wlsel/10

# # Load Kurucz FTS Atlas
wlsel, ilamsel = DKISTanalysis.load_fts(dispersion_range-.101)
wlsel=wlsel/10


# # Average the QS data for space and time, for selected ranges
space_and_time_averaged_qs = \
    DKISTanalysis.comp_fts_to_qs(wlsel,ilamsel,dispersion_range, 
                                  image_data_arr_arr_qs)

# # telluric lines for comparison (or other absorption lines if telluric not 
# # available, as is the case for the Ca II H window).  Most of the next steps
# # are not used for pid_1_38, but include in any case to model the use of lines
line1 = 485.97411 #fe I
line2 = 486.2610 #v1

# # Definition of "telluric" or other absorption lines for calibration.
# # The QS spectrum is already averaged in space and time, so any absorption lines
# # (as is necessary for the 19 August 2022 observations) are ok

# # indices of telluric lines in spectrum - lower
lowinds = [390,666]

# # indices of telluric lines in spectrum - upper
highinds = [455,700]

cont_vals = [485.636,485.666,485.762,485.796,485.866,485.942,486.014,486.061,486.213,486.286,486.341,486.453,486.454,486.494]

# # define the multiplication factor (polynomial), new dispersion range, fit values
# # to scale the quiet sun to FTS atlas
cont_mult_facts,fit_vals,new_dispersion_range,dispersion_range_fin,rat=\
    DKISTanalysis.get_calibration_poly(dispersion_range,
                                        space_and_time_averaged_qs,wlsel,ilamsel,
                                        DKISTanalysis.find_nearest,line1,line2,
                                        lowinds,highinds,limbdark_fact=clv_corrqs,
                                        cont_vals=cont_vals,order=1)

# # calibrate the quiet sun 
calibrated_qs=fit_vals*space_and_time_averaged_qs/clv_corrqs
nonflare_average_avg = calibrated_qs
nonflare_multfact = fit_vals

# # full width of half max of PSF to convolve with atlas to match instrument
fwhm = 0.05

# # number of points to interpolate Atlas to in PSF convolve to match instrument
ntw = 45

# # perform PSF convolution.  Result 'yconv' is Atlas*PSF
yconv=DKISTanalysis.psf_adjust(wlsel,ilamsel,fwhm,new_dispersion_range,
                                calibrated_qs,clv_corrqs,ntw,
                                DKISTanalysis.gaussian_psf)



#do another iteration of the calibration step after peforming the PSF conv.
cont_mult_facts,fit_vals,\
    new_dispersion_range,dispersion_range_fin,rat=\
        DKISTanalysis.get_calibration_poly(dispersion_range,
                                            space_and_time_averaged_qs,
                                            new_dispersion_range,yconv,
                                            DKISTanalysis.find_nearest,
                                            line1,line2,lowinds,highinds,
                                            limbdark_fact=clv_corrqs,noqs_flag=2,
                                            cont_vals=cont_vals)
        
        

        
# calibrated quiet sun, again, using updated fit values
calibrated_qs=fit_vals*space_and_time_averaged_qs/clv_corrqs
nonflare_average_avg = calibrated_qs
nonflare_multfact = fit_vals

#show comparison of atlas to qs
fig,ax=plt.subplots()
ax.plot(dispersion_range_fin,calibrated_qs,label='visp')
ax.plot(dispersion_range_fin,yconv*clv_corrqs,label='convolved')
ax.plot(wlsel,clv_corrqs*ilamsel,label='raw')
ax.set_xlim([485.5,487]);
ax.legend();plt.show()

# intensity calibration, background subtraction for flare-time                            
scaled_flare_time, bkgd_subtract_flaretime = \
    DKISTanalysis.scaling(image_data_arr_arrhbeta, nonflare_multfact,clv_corrhbeta,
                          nonflare_average_avg,end=end)


low_foravg = 500
high_foravg = 650

# # definition of index bounds for Ca II H lines
# caII_low_foravg = 525
# caII_high_foravg = 600

# #indices to search within to find flare kernel center
# spacelow = 1000
# spacehigh = -1

# #indices of max intensity in each frame
# maxindices = DKISTanalysis.maxintind(dispersion_range,bkgd_subtract_flaretime,
#                                      caII_low_foravg,caII_high_foravg,
#                                      spacelow,spacehigh)

# ### remove if desire tracking of central position ###
# #testing with just the initial brightest point in the ribbon
# # firstmaxonly = []
# # for i in range(len(maxindices)):
# #     firstmaxonly.append(maxindices[0])
# # maxindices = firstmaxonly
# ####

# ### remove if desire tracking of central position ###
# #testing with edge of flare ribbon
# # edgeind = []
# # for i in range(len(maxindices)):
# #     edgeind.append(maxindices[i]+40)
# # maxindices = edgeind
# ####


# # plot intensity calibrated, background-subtracted spectra
# DKISTanalysis.pltsubtract(dispersion_range_fin,nonflare_average_avg,
#                           scaled_flare_time,muted,maxindices,end=end,pid='pid_1_50')

# # variation in intensity value corresponding to wavelengths; PTE; to test
# # for variations in pseudo-continuum.  If PTE high, cannot be explained by the
# # solar activity difference between quiet sun and flare-time

# ########
# ## Comment out this block if don't need PTE (shouldn't, usually)
# # stdevs_flaretime, ptes_flaretime = \
#     # DKISTanalysis.deviations(bkgd_subtract_flaretime,nonflare_average,
#                              # nonflare_stdevs)

# # DKISTanalysis.pltptes(ptes_flaretime,image_data_arr_arr_raster1)

# ########

# # equivalent widths, effective widths, widths
# # indices defining line locations
# caII_low = 480
# caII_high = 650
# hep_low = 700
# hep_high = 850

# #first profile, to use in testing
# sample_flaretime = bkgd_subtract_flaretime[0,:,maxindices[0]]

# ########
# # perform following block only if need to calculate continuum window 
# # independently of width determiation; 
# # exists within width determination script as well

# #nolines, cont_int_array, cont_int_wave_array = \
#     # DKISTanalysis.contwind(sample_flaretime,dispersion_range,maxinds,avgs,
#                            # low,high)
                           
# ########

# # line widths, strengths - initialize arrays
# ew_CaII_all_fs = np.zeros((len(scaled_flare_time)-5,
#                            np.shape(bkgd_subtract_flaretime)[2]))
# ew_hep_all_fs = np.zeros((len(scaled_flare_time)-5,
#                           np.shape(bkgd_subtract_flaretime)[2]))
# eqw_CaII_all_fs = np.zeros((len(scaled_flare_time)-5,
#                             np.shape(bkgd_subtract_flaretime)[2]))
# eqw_hep_all_fs = np.zeros((len(scaled_flare_time)-5,
#                            np.shape(bkgd_subtract_flaretime)[2]))
# width_CaII_all_fs = np.zeros((len(scaled_flare_time)-5,
#                               np.shape(bkgd_subtract_flaretime)[2]))
# width_hep_all_fs = np.zeros((len(scaled_flare_time)-5,
#                              np.shape(bkgd_subtract_flaretime)[2]))


# # line widths, strength determination
# ew_CaII_all_fs, ew_hep_all_fs, eqw_CaII_all_fs,\
#     eqw_hep_all_fs, width_CaII_all_fs, width_hep_all_fs = \
#         DKISTanalysis.widths_strengths(ew_CaII_all_fs,eqw_CaII_all_fs,
#                                        width_CaII_all_fs,ew_hep_all_fs,
#                                        eqw_hep_all_fs,width_hep_all_fs,caII_low,
#                                        caII_high,hep_low,hep_high,
#                                        scaled_flare_time,
#                                        bkgd_subtract_flaretime, 
#                                        dispersion_range_fin)
        
# # Gaussian fitting
# # automate for all timesteps
# storeamp1 = []
# storemu1 = []
# storesig1 = []
# storeamp2 = []
# storemu2 = []
# storesig2 = []

# # spatial index corresponding to part of observation of interest
# sliceind = maxindices[0]
# sel = bkgd_subtract_flaretime[0,caII_low:caII_high,sliceind]-\
#     min(bkgd_subtract_flaretime[0,caII_low:caII_high,sliceind])
# selwl = new_dispersion_range[caII_low:caII_high]
        
# # parameter guesses
# parameters=[2e6,396.82,0.015,6e6,396.86,0.015]

# # fitting of spectral lines, e.g. double-Gaussian (for Ca II H)
# storeamp1,storeamp2,storesig1,storesig2,storemu1,storemu2 = \
#     DKISTanalysis.gauss2fit(storeamp1,storemu1,storesig1,storeamp2,storemu2,
#                             storesig2,bkgd_subtract_flaretime,dispersion_range_fin,
#                             DKISTanalysis.double_gaussian_fit,maxindices,times,
#                             caII_low,caII_high,DKISTanalysis.double_gaussian,
#                             DKISTanalysis.gaussian,selwl,sel,
#                             parameters = parameters,pid='pid_1_50')

# # width determination
# store_ten_width = []
# store_quarter_width = []
# store_half_width = []

# store_ten_width, store_quarter_width, store_half_width = \
#     DKISTanalysis.perclevels(bkgd_subtract_flaretime,dispersion_range_fin,caII_low,
#                              caII_high,store_ten_width,store_quarter_width,
#                              store_half_width)

# storeamp1_2 = []
# storemu1_2 = []
# storesig1_2 = []
# storeamp2_2 = []
# storemu2_2 = []
# storesig2_2 = []

# # output fit parameters
# fits_1g,fits_2g,fits_2gneg,params2gaussnew,stopind,storeamp1_2,\
#     storemu1_2,storesig1_2,storeamp2_2,storemu2_2,storesig2_2,selwl,sel = \
#     DKISTanalysis.fittingroutines(bkgd_subtract_flaretime,dispersion_range_fin,
#                                   times, caII_low, caII_high,
#                                   DKISTanalysis.double_gaussian, 
#                                   DKISTanalysis.gaussian, 
#                                   selwl,sel,[4e6,396.85,0.02],
#                                   parameters,
#                                   [.5e6,396.85,0.015,-1e6,396.85,0.015],
#                                   maxindices,storeamp1_2,storemu1_2,
#                                   storesig1_2,storeamp2_2,storemu2_2,
#                                   storesig2_2,pid='pid_2_11', date = '8/8/2024',
#                                   line = 'Ca II H',nimg = 8,
#                                   inds=[380,390,400,410,450,480,647,700,820,850,900],deg=7)

# vel1,vel2 = DKISTanalysis.conv_to_vel(storemu1_2,storemu2_2,mu)
# # plot results of Gaussian fitting

# # note for this particular go-around
# note = ', adjust with mu angle'

# #plot results of fitting
# DKISTanalysis.pltfitresults(bkgd_subtract_flaretime,dispersion_range_fin,
#                             DKISTanalysis.double_gaussian,
#                             DKISTanalysis.gaussian,times,muted,
#                             caII_low,caII_high,fits_1g,fits_2g,fits_2gneg,maxindices,
#                             mu, pid='pid_1_84', date = '08092022',line = 'Ca II H',
#                             nimg = 7, nrow=2,ncol=4,lamb0=wl,note=note,yhigh=7.5e6,
#                             inds=[380,390,400,410,450,480,647,700,820,850,900],deg=7)


# #processing of all raster slices
# scaled_flare_time_allrasters, bkgd_subtract_flaretime_allrasters = \
#     DKISTanalysis.scaling(image_data_arr_arr, nonflare_multfact,clv_corr,
#                           nonflare_average_avg,end=end)
    
# maxindices_allrasters = DKISTanalysis.maxintind(dispersion_range,
#                                                 bkgd_subtract_flaretime_allrasters,
#                                                 caII_low_foravg,caII_high_foravg,
#                                                 spacelow,spacehigh)
   




