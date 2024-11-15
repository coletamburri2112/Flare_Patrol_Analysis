#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
4 December 2023
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
    
    Applied to pid_1_38, observations from 20 April 2022, Ca II 8542

"""

# package initialize
import dkistpkg_ct as DKISTanalysis
import numpy as np
import matplotlib

matplotlib.rcParams.update(matplotlib.rcParamsDefault)

# establish plotting methods and formats
# plt.rcParams['text.usetex']=True
# plt.rcParams['font.family']='sans-serif'
# plt.rcParams['font.sans-serif'] = ['Helvetica']
# plt.rcParams['axes.labelsize'] = 25
# plt.rcParams['lines.linewidth'] = 2
# matplotlib.rc('xtick', labelsize=20) 
# matplotlib.rc('ytick', labelsize=20) 

# color scheme for plotting
muted = DKISTanalysis.color_muted2()

# path and file ID for ViSP data
path = '/Volumes/ViSP_External/pid_1_38/'
folder1 = ''

# list of files in directory for DKIST/ViSP
dir_list2 = DKISTanalysis.pathdef(path,folder1)

# Stonyhurst lon/lat position of the AR from JHv
lon = -27.69 #degrees
lat = 18.29 #degrees

wl = 854.209 # central wavelength, Ca II 854.209, in nm

# spatial coordinates
hpc1_arcsec, hpc2_arcsec, x_center, y_center, z, rho, mu, \
    doppshnonrel, doppshrel = \
    DKISTanalysis.spatialinit(path,folder1,dir_list2,lon,lat,wl)

# get limb darkening coefficient 
clv_corr = DKISTanalysis.limbdarkening(wl, mu=mu, nm=True)
    # for Ca II H (require mu value for determination, be sure to specify
    # correct wl units)
    
# two different starting coefficients - one for QS observations, one for 
# flare-time observations
startstepqs = 0
startstepflare = 2720

# process multi-step raster - for qs time
image_data_arr_arr_qs, rasterpos_qs, times_qs = \
    DKISTanalysis.multistepprocess(path,folder1,dir_list2,div=1000,
                                   startstep=startstepqs)
    
# process multi-step raster - flaretime
image_data_arr_arr, rasterpos, times = \
    DKISTanalysis.multistepprocess(path,folder1,dir_list2,div=20,
                                   startstep=startstepflare)
    
# spatial and dispersion axes for single observation (single slit step)
spatial_range, dispersion_range = \
    DKISTanalysis.spatialaxis(path,folder1,dir_list2,line='Ca II 8542')
    
#qs calibration
wlsel, ilamsel = DKISTanalysis.load_fts(dispersion_range)

wlsel = wlsel/10
space_and_time_averaged_qs = \
    DKISTanalysis.comp_fts_to_qs(wlsel,ilamsel,dispersion_range, 
                                 image_data_arr_arr_qs, lowint = -100,highint=-1,
                                 timelow=0,timehigh=5)

# telluric lines for comparison (or other absorption lines if telluric not 
# available, as is the case for the Ca II H window).  Most of the next steps
# are not used for pid_1_38, but include in any case to model the use of lines
line1 = 853.801
line2 = 854.804

# indices of telluric lines in spectrum - lower
lowinds = [277,737]

# indices of telluric lines in spectrum - upper
highinds = [297,754]

# dispersion range and multiplicative ratio from intensity calibration
new_dispersion_range, rat = \
    DKISTanalysis.calib_qs_shift(wlsel,ilamsel,dispersion_range,
                                 space_and_time_averaged_qs, 
                                 line1,line2,lowinds,highinds)

# final calibration - the output of above function is good as "initial guess"
# to find best chi-squared fit with this function, below 

# but for pid_1_38, observations are not in a quiet sun region, so this function
# really just takes a single index to compare the atlas and the observations
# and uses that as the multiplicative factor
calibration, calibrated_qs,new_dispersion_range2= \
    DKISTanalysis.get_calibration_singleval(dispersion_range,\
                                            space_and_time_averaged_qs,\
                                                wlsel,ilamsel,\
                                                    limbdark_fact=clv_corr,
                                                    noqs_flag = 1)

# plotting intensity calibration results
DKISTanalysis.plot_calibration(new_dispersion_range2,calibrated_qs,wlsel,ilamsel,
                               pid='pid_1_38')

# load QS intensity calibration results here!
nonflare_average = calibrated_qs
nonflare_multfact = np.full(len(dispersion_range), calibration[0])

# intensity calibration, background subtraction                            
scaled_flare_time, bkgd_subtract_flaretime = \
    DKISTanalysis.scaling(image_data_arr_arr, nonflare_multfact, clv_corr,
                          nonflare_average,end=1)


# equivalent widths, effective widths, widths
caII_8542_low_foravg = 535
caII_8542_high_foravg = 565

caII_8542_low = 500
caII_8542_high = 600

caII_8542_low_for_fit = 490
caII_8542_high_for_fit = 610

spacelow = 1200
spacehigh = 2000

#indices of max intensity - presumably spatial indices of flare kernel
maxindices = DKISTanalysis.maxintind(new_dispersion_range2,image_data_arr_arr,
                                     caII_8542_low_foravg,caII_8542_high_foravg,
                                     spacelow,spacehigh)

# plot intensity calibrated, background-subtracted spectra
DKISTanalysis.pltsubtract(new_dispersion_range2,nonflare_average,
                          scaled_flare_time,muted,maxindices,pid='pid_1_38')


# put bkgd_subtract_flaretime here, when ready
sample_flaretime = bkgd_subtract_flaretime[0,:,maxindices[0]]

# perform following line only if need to calculate continuum window 
# independently of width determiation; 
# exists within width determination script as well

#nolines, cont_int_array, cont_int_wave_array = \
    # DKISTanalysis.contwind(sample_flaretime,dispersion_range,maxindices,avgs,
                           # low,high)

end = 1
# line widths, strengths initial arrays
ew_CaII_all_fs = np.zeros((len(scaled_flare_time)-end,
                            np.shape(bkgd_subtract_flaretime)[2]))
eqw_CaII_all_fs = np.zeros((len(scaled_flare_time)-end,
                            np.shape(bkgd_subtract_flaretime)[2]))
width_CaII_all_fs = np.zeros((len(scaled_flare_time)-end,
                              np.shape(bkgd_subtract_flaretime)[2]))
int_CaII_all_fs = np.zeros((len(scaled_flare_time)-end,
                              np.shape(bkgd_subtract_flaretime)[2]))

altinds=[410,460,510,520,590,600,610,620,630,640,650,660,670,680,690,700,
      720]

# # line widths, strength determination
ew_CaII_all_fs, eqw_CaII_all_fs,\
    width_CaII_all_fs, int_CaII_all_fs= \
        DKISTanalysis.widths_strengths_oneline(ew_CaII_all_fs,eqw_CaII_all_fs,
                                                width_CaII_all_fs,int_CaII_all_fs,
                                                caII_8542_low,caII_8542_high,
                                                scaled_flare_time,
                                                bkgd_subtract_flaretime,
                                                dispersion_range,maxindices,
                                                deg=7,low0=105,
                                                high0=134,low1=185,high1=202,
                                                low2=290,high2=306,low3=357,
                                                high3=370,low4=446,high4=464,
                                                low5=616,high5=645,low6=692,
                                                high6=715,alt=1,altinds=altinds)
        
# plot results of line characteristics - widths, strengths, fluxes
DKISTanalysis.plt_line_characteristics(ew_CaII_all_fs,eqw_CaII_all_fs,
                                       width_CaII_all_fs, int_CaII_all_fs,
                                       maxindices,times,muted,pid='pid_1_38',
                                       nslitpos=45,raster=0,nframes=10,
                                       line='Ca II 8542')
        
# Gaussian fitting
# automate for all timesteps
storeamp1 = []
storemu1 = []
storesig1 = []
storeamp2 = []
storemu2 = []
storesig2 = []


sel = bkgd_subtract_flaretime[0,caII_8542_low:caII_8542_high,maxindices[0]]-\
    min(bkgd_subtract_flaretime[0,caII_8542_low:caII_8542_high,maxindices[0]])
selwl = dispersion_range[caII_8542_low:caII_8542_high]
        
# # width determination - percentage levels - untested
# store_ten_width = []
# store_quarter_width = []
# store_half_width = []

# store_ten_width, store_quarter_width, store_half_width = \
#     DKISTanalysis.perclevels(bkgd_subtract_flaretime,new_dispersion_range2,
#                              caII_8542_low,caII_8542_high,store_ten_width,
#                              store_quarter_width,store_half_width)
    
nimg = 100
# output fit parameters
fits_1g,fits_2g,fits_2gneg,params2gaussnew,stopind= \
    DKISTanalysis.fittingroutines(bkgd_subtract_flaretime,new_dispersion_range2,
                                  times, caII_8542_low_for_fit, 
                                  caII_8542_high_for_fit,
                                  DKISTanalysis.double_gaussian, 
                                  DKISTanalysis.gaussian, 
                                  selwl,sel,[1.1,854.28,0.05],
                                  [.3e6,854.2,0.01,.5e6,854.22,0.015],
                                  [.5e6,396.85,0.015,-1e6,396.85,0.015],
                                  maxindices,pid='pid_1_38', date = '04/20/2022',
                                  line = 'Ca II 854.2',nimg = nimg)
    
# plot results of Gaussian fitting
# yhigh will set the maximum for 
DKISTanalysis.pltfitresults(bkgd_subtract_flaretime,new_dispersion_range2,
                            DKISTanalysis.double_gaussian,
                            DKISTanalysis.gaussian,times,muted,caII_8542_low,
                            caII_8542_high,fits_1g,fits_2g,fits_2gneg,
                            maxindices,pid='pid_1_38', date = '04202022',
                            line = 'Ca 854.2', nimg =nimg, nrow=4,ncol=5,
                            note='_only_flare_12_13_2023',lim=0.1,lamb0=wl)
    




