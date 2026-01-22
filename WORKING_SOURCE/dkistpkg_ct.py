# -*- coding: utf-8 -*-
"""
2 December 2023
Author: Cole Tamburri, University of Colorado Boulder, National Solar 
Observatory, Laboratory for Atmospheric and Space Physics

Description of script: 
    Main working functions for analysis package of DKIST data, applied first to 
    pid_1_84 ("flare patrol"), with PI Kowalski and Co-Is Cauzzi, Tristain, Notsu, 
    Kazachenko, and (unlisted) Tamburri.  Also applied to pid_2_11, with nearly 
    identical science objectives.  See READMe for details.  
    
    Includes intensity calibration, Gaussian fitting, and co-alignment routes 
    between ViSP and VBI, and, externally, SDO/HMI.  This code was developed with 
    the ViSP Ca II H and VBI blue continuum, TiO, and H-alpha channels as priority, 
    and using HMI in the continuum and 304 Angstrom bandpasses, but there is room 
    for application to other channels and science objectives.


"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import scipy.integrate as integrate
import sunpy
import sunpy.coordinates
#import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
import astropy.units as u
import astropy
from astropy.table import Table
from scipy.optimize import differential_evolution
from scipy.interpolate import interp1d
#from scipy.signal import convolve
import matplotlib.pylab as pl 
from datetime import datetime
from astropy.convolution import convolve, Gaussian1DKernel


# from sunpy.net import Fido, attrs as a
# import pandas as pd
# from astropy.utils.data import get_pkg_data_filename
# import shutil
# import fitsio
# import matplotlib.animation as animat
# import ffmpeg
# import latex
# import radx_ct
# import math as math
# import scipy.special as sp
# from scipy.stats import skewnorm
# from lmfit.models import SkewedGaussianModel
# import matplotlib
# from matplotlib import animation
# from lmfit import Model
# from pylab import *
# from astropy.coordinates import SkyCoord
# from astropy.time import Time
# from astropy.visualization import ImageNormalize, SqrtStretch



# define functions to be used for line fitting
def gaussian(x, c1, mu1, sigma1):
    res = c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) )
    return res

def gaussfit(params,selwl,sel):
    fit = gaussian( selwl, params )
    return (fit - sel)

def double_gaussian( x, c1, mu1, sigma1, c2, mu2, sigma2 ):
    res =   (c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) )) \
          + (c2 * np.exp( - (x - mu2)**2.0 / (2.0 * sigma2**2.0) ))
    return res

def double_gaussian_fit( params ,selwl,sel):
    fit = double_gaussian( selwl, params[0],params[1],params[2],params[3],
                          params[4],params[5] )
    return (fit - sel)

def numgreaterthan(array,value):
    count = 0
    for i in array:
        if i > value:
            count = count + 1
    return count

# for calibration
def gaussian_for_calib(x, A, mu, sigma, baseline):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))+ baseline

def fitline(xdata, ydata, plot='True'):
	#Gaussian fit
	params, covariance = curve_fit(gaussian_for_calib, xdata, ydata, p0=[-0.2, np.mean(xdata), 1,0.4])
	a,b,c,d = params
	xfit  = np.linspace(xdata[0],xdata[-1], 1000)
	yfit = gaussian_for_calib(xfit, a,b,c,d)

	w0 = xfit[np.argmin(yfit)]
	print(' MINIMA loc--->',w0)
	
	if plot=='True':
		plt.plot(xdata, ydata, '.')
		plt.plot(xfit, yfit)
		plt.plot([w0,w0], [yfit.min()-0.1, yfit.max()], color='gray')
		plt.ylim(ydata.min()-0.03, ydata.max()+0.03)

	return xfit, yfit, w0


def color_muted2():
    #define colors for plotting
    
    #  0= indigo
    # 1 = cyan
    # 2 = teal
    # 3 = green
    # 4 = olive
    # 5= sand
    # 6 = rose
    # 7 = wine
    # 8 = purple
    # 9=grey
    
    muted =['#332288', '#88CCEE', '#44AA99','#117733','#999933','#DDCC77', '#CC6677','#882255','#AA4499','#DDDDDD']

    return muted

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx


def pathdef(path,folder1,stokesi=1):
    
    # define path for DKIST data; assumes L1 .fits files; 

    dir_list = os.listdir(path+folder1)
    
    #redefine list of directory contents
    dir_list2 = []
    for i in range(len(dir_list)):
        filename = dir_list[i]
        
        # modify based on filename structure and Stokes file preference
        if stokesi == 1 and filename[-5:] == '.fits' and '_I_' in filename:
            dir_list2.append(filename)
        elif stokesi == 0 and filename[-5:] == '.fits':
            dir_list2.append(filename)
    
    dir_list2.sort()
    
    return dir_list2

def spatialinit(path,folder1,dir_list2,lon,lat,wl, flag=0):
    
    # initialize spatial parameters, mu angle for use in determining limb darkening

    i_file_raster1 = fits.open(path+folder1+'/'+dir_list2[0]) #first image frame
    d = 150.84e9 #average distance to sun, in meters
    solrad = sunpy.sun.constants.radius.value
    
    # Coordinates from DKIST are not correct, but define them anyways as a starting
    # point.  Will co-align later in routine.
    
    if flag == 1: #if pid_1_84
        hpc1_arcsec = i_file_raster1[1].header['CRVAL2'] # first axis coords # for pid_1_84
    else:
        hpc1_arcsec = i_file_raster1[1].header['CRVAL1'] # first axis coords # for pid_1_84
        
    hpc2_arcsec = i_file_raster1[1].header['CRVAL3'] # second axis corods
    

    if flag == 2: # if QS for 22:34 on 19 August; coordinates are wrong
        print('here')
        hpc1_arcsec = -445
        hpc2_arcsec= 711

    # image center
    x_center = d*np.cos(hpc1_arcsec/206265)*np.sin(hpc2_arcsec/206265) # m
    y_center = d*np.sin(hpc1_arcsec/206265) # m
    z = solrad - d*np.cos(hpc1_arcsec/206265)*np.cos(hpc2_arcsec/206265) # m
    
    # to mu value
    rho = np.sqrt(x_center**2+y_center**2)
    mu = np.sqrt(1-(rho/solrad)**2)
    
    rotrate = 13.2 # approx. rotation rate; degrees of longitude per day
    
    #convert to cylindrical coordinates
    rcyl = solrad*np.cos(-lat*2*np.pi/360) #radius of corresponding cylinder
    circyl = 2*np.pi*solrad
    dratecyl = rotrate*circyl/360
    dratecylms = dratecyl/24/3600 # meters/day to meters/s
    
    # potential redshift effects due to solar rotation
    # testing impact of this on spectra to account for perceived wl shift in
    # data? Seems to account for shift in pid_1_84 by chance, uncertain if this
    # is the solution to the issue
    redsh = dratecylms*np.cos((90-lon)*2*np.pi/360) #redshift value
    
    # doppler shift in line due to rotation
    doppshnonrel = wl*(1+(redsh/3e8)) - wl #non-relativistic
    doppshrel = wl*np.sqrt(1-redsh/3e8)/np.sqrt(1+redsh/3e8) - wl #relativistic
        
    return hpc1_arcsec, hpc2_arcsec, x_center, y_center, z, rho, mu, \
        doppshnonrel, doppshrel
        
def limbdarkening(wave, mu=1.0, nm=False, 
                  directory = '/Users/coletamburri/Desktop/DKIST_Data_Tools_misc/Rahul_ViSP_Cal/'):
    """
    Return limb-darkening factor given wavelength and viewing angle
    mu=cos(theta)
    Arguments:
        wave: scalar or 1D array with wavelength(s).
    Keyword arguments:
        mu: cosine of heliocentric viewing angle (default 1.0 -> disc centre)
        nm: input wavelength units are nanometers (default False)
    Returns:
        factor: scaling factor(s) to be applied for given input wavelengths. Has
        as many elements as `wave`.
    Example:
        factor = limbdarkening(630.25, mu=0.7, nm=True)
    Author:
        Gregal Vissers (ISP/SU 2020)
    """

    DATA_PATH = os.path.join(directory, "limbdarkening_Neckel_Labs_1994.fits")

    wave = np.atleast_1d(wave)  # Ensure input is iterable

    table = astropy.table.Table(fits.getdata(DATA_PATH))
    wavetable = np.array(table['wavelength'])
    if nm is False:
        wavetable *= 10.

    # Get table into 2D numpy array
    Atable = np.array([ table['A0'], table['A1'], table['A2'],
        table['A3'], table['A4'], table['A5'] ])
    
    factor = np.zeros((wave.size), dtype='float64')
    for ii in range(6):
      Aint = np.interp(wave, wavetable, Atable[ii,:])
      factor += Aint * mu**ii

    return factor[0]



def fourstepprocess(path,folder1,dir_list2,fullstokes=1):
    
    # Simplest initial processing of data, when only a four-step raster; will 
    # certainly need to be generalized to observations with more raster steps in
    # ViSP observations

    image_data_arrs_raster1 = []
    image_data_arrs_raster2 = []
    image_data_arrs_raster3 = []
    image_data_arrs_raster4 = []
    rasterpos1 = []
    rasterpos2 = []
    rasterpos3 = []
    rasterpos4 = []
    times_raster1 = []
    times_raster2 = []
    times_raster3 = []
    times_raster4 = []
    times = []

    image_data_arrs = []
    
    #four raster steps; make array for each
    if fullstokes==1:
        for i in range(0,len(dir_list2),4):
            i_file_raster1 = fits.open(path+folder1+'/'+dir_list2[i])
            q_file_raster1 = fits.open(path+folder1+'/'+dir_list2[i+1])
            u_file_raster1 = fits.open(path+folder1+'/'+dir_list2[i+2])
            v_file_raster1 = fits.open(path+folder1+'/'+dir_list2[i+3])
            
            times_raster1.append(i_file_raster1[1].header['DATE-BEG'])
        
            
            times.append(i_file_raster1[1].header['DATE-BEG'])
            times.append(q_file_raster1[1].header['DATE-BEG'])
            times.append(u_file_raster1[1].header['DATE-BEG'])
            times.append(v_file_raster1[1].header['DATE-BEG'])
            
            i_data_raster1 = i_file_raster1[1].data[0]
            q_data_raster1 = q_file_raster1[1].data[0]
            u_data_raster1 = u_file_raster1[1].data[0]
            v_data_raster1 = v_file_raster1[1].data[0]
            
            # # collect observations belonging to the same slit position
            # image_data_arrs_raster1.append(i_data_raster1)
            # image_data_arrs_raster2.append(i_data_raster2)
            # image_data_arrs_raster3.append(i_data_raster3)
            # image_data_arrs_raster4.append(i_data_raster4)
            
            # rasterpos1.append(i_file_raster1[1].header['CRPIX3'])
            # rasterpos2.append(i_file_raster2[1].header['CRPIX3'])
            # rasterpos3.append(i_file_raster3[1].header['CRPIX3'])
            # rasterpos4.append(i_file_raster4[1].header['CRPIX3'])
            
            #array including all raster positions
            image_data_arrs.append(i_data_raster1)
            image_data_arrs.append(q_data_raster1)
            image_data_arrs.append(u_data_raster1)
            image_data_arrs.append(v_data_raster1)
    else:
        for i in range(0,len(dir_list2),4):
            i_file_raster1 = fits.open(path+folder1+'/'+dir_list2[i])
            i_file_raster2 = fits.open(path+folder1+'/'+dir_list2[i+1])
            i_file_raster3 = fits.open(path+folder1+'/'+dir_list2[i+2])
            i_file_raster4 = fits.open(path+folder1+'/'+dir_list2[i+3])
            
            times_raster1.append(i_file_raster1[1].header['DATE-BEG'])
            times_raster2.append(i_file_raster2[1].header['DATE-BEG'])
            times_raster3.append(i_file_raster3[1].header['DATE-BEG'])
            times_raster4.append(i_file_raster4[1].header['DATE-BEG'])
            
            times.append(i_file_raster1[1].header['DATE-BEG'])
            times.append(i_file_raster2[1].header['DATE-BEG'])
            times.append(i_file_raster3[1].header['DATE-BEG'])
            times.append(i_file_raster4[1].header['DATE-BEG'])
            
            i_data_raster1 = i_file_raster1[1].data[0]
            i_data_raster2 = i_file_raster2[1].data[0]
            i_data_raster3 = i_file_raster3[1].data[0]
            i_data_raster4 = i_file_raster4[1].data[0]
            
            # collect observations belonging to the same slit position
            image_data_arrs_raster1.append(i_data_raster1)
            image_data_arrs_raster2.append(i_data_raster2)
            image_data_arrs_raster3.append(i_data_raster3)
            image_data_arrs_raster4.append(i_data_raster4)
            
            rasterpos1.append(i_file_raster1[1].header['CRPIX3'])
            rasterpos2.append(i_file_raster2[1].header['CRPIX3'])
            rasterpos3.append(i_file_raster3[1].header['CRPIX3'])
            rasterpos4.append(i_file_raster4[1].header['CRPIX3'])
            
            #array including all raster positions
            image_data_arrs.append(i_data_raster1)
            image_data_arrs.append(i_data_raster2)
            image_data_arrs.append(i_data_raster3)
            image_data_arrs.append(i_data_raster4)
        
            # array version of images corresponding to each slit position
            image_data_arr_arr_raster1 = np.array(image_data_arrs_raster1)
            image_data_arr_arr_raster2 = np.array(image_data_arrs_raster2)
            image_data_arr_arr_raster3 = np.array(image_data_arrs_raster3)
            image_data_arr_arr_raster4 = np.array(image_data_arrs_raster4)
    
    # all rasters
    image_data_arr_arr = np.array(image_data_arrs)
    
    # for intensity calibration purposes, only images from first raster pos
    if fullstokes == 0:
        for_scale = image_data_arr_arr_raster1[:,:,:]
    else:
        for_scale = image_data_arr_arr[0,:,:]
    # uncomment second and third lines to return individual raster arrays
    return image_data_arr_arr, i_file_raster1, for_scale, times_raster1,times#,\
        # image_data_arr_arr_raster1, image_data_arr_arr_raster2,\
        # image-data_arr-arr_raster3, image_data_arr_arr_raster3
        
def multistepprocess(path,folder,dir_list,div=10,startstep=0,endstep=-1):
    
    # Multi-step raster procesing of DKIST data; bare bones storage, no
    # separation based on slit step, output includes all slit step positions
    
    # div - dividing factor, proportion of directory to search through for data
    
    # startstep - number of files to skip in directory before begining storage

    image_data_arrs0 = []
    rasterpos = []
    times = []
    # extract relevant information
    for i in range(startstep,endstep-1,1):
        print(i)
        
        i_file = fits.open(path+folder+'/'+dir_list[i])
        times.append(i_file[1].header['DATE-BEG'])
        i_data = i_file[1].data[0]
        image_data_arrs0.append(i_data)
        rasterpos.append(i_file[1].header['CRPIX3'])
    if endstep==-1:
        for i in range(startstep,len(dir_list)-1,1):
            print(i)
            
            i_file = fits.open(path+folder+'/'+dir_list[i])
            times.append(i_file[1].header['DATE-BEG'])
            i_data = i_file[1].data[0]
            image_data_arrs0.append(i_data)
            rasterpos.append(i_file[1].header['CRPIX3'])
    
    # all rasters
    image_data_arr_arr = np.array(image_data_arrs0)

    
    return image_data_arr_arr, rasterpos, times


def spatialaxis(path,folder1,dir_list2,line='Ca II H',pid='84',shift=0):
    
    # find the axes of ViSP observations based on values given in L1 header;
    # spectral axis can be trusted as long as DKIST data set caveats have been
    # accounted for ( https://nso.atlassian.net/wiki/spaces/DDCHD/pages/1959985377/DKIST+Data+Set+Caveats+ViSP+VBI ).
    # Spatial coordiantes should not be trusted; only to be used for co-align
    # routines with SDO and between ViSP/VBI
    
    i_file_raster1 = fits.open(path+folder1+'/'+dir_list2[0])
    
    #crval1,cdelt1
    hdul1 = i_file_raster1
    centerlambda = hdul1[1].header['CRVAL1']
    

    deltlambda = hdul1[1].header['CDELT1'] 
    nlambda = hdul1[1].header['NAXIS2']
    
    if pid == '2_11':
        centerlambda = hdul1[1].header['CRVAL2']
        deltlambda = hdul1[1].header['CDELT2']       
    if pid == '50':
        centerlambda = hdul1[1].header['CRVAL2']
        deltlambda = hdul1[1].header['CDELT2']
    
    dispersion_range = np.linspace(centerlambda-deltlambda*(nlambda-1)/2,
                                   centerlambda+deltlambda*(nlambda-1)/2,nlambda)
    
    if line == 'Ca II 8542':
        dispersion_range = []
        for i in range(nlambda):
            dispersion_range.append(853.1749+1.88966e-3*i-7.506692e-9*i*i)
    
    centerspace = hdul1[1].header['CRVAL2']
    deltspace = hdul1[1].header['CDELT2'] #this actually is fine
    nspace = hdul1[1].header['NAXIS1']
    spatial_range = np.linspace(centerspace-deltspace*(nspace-1)/2,
                                centerspace+deltspace*(nspace-1)/2,nspace)
    
    if pid == '2_11':
        centerspace = hdul1[1].header['CRVAL1']
        deltspace = hdul1[1].header['CDELT1'] #this actually is fine
        nspace = hdul1[1].header['NAXIS1']
        spatial_range = np.linspace(centerspace-deltspace*(nspace-1)/2,
                                    centerspace+deltspace*(nspace-1)/2,nspace)   
        
    if shift > 0:
        dispersion_range = dispersion_range-shift
    return spatial_range, dispersion_range

       

def scaling(for_scale,nonflare_multfact,limbdarkening,nonflare_average,
            limbd = 1,end=5):
    # Scaling relative to QS values.  For this, require inputs of "nonflare" -
    # this can take the form of off-kernel observations.  In our case, was a disk
    # center observation, hence the allowance for limb darkening correction.  
    # Disk-center QS were compared to the Neckel-Hamburg disk center atlas, and a
    # multiplicative factor applied to give intensity values of observations. Based
    # on the spectral range being studied, may (1) have a single scaling factor, as
    # is found by e.g. Rahul Yadav, or (2) a variable scaling factor, as applied 
    # below, found for the Ca II H window by Cole Tamburri.  Solution to this 
    # theoretical dilemma not found as of 2 Dec 2023
    
    # set limbd to 0 if QS/nonflare values used for calibration were not at
    # disk center (or if limb darkening accounted for at a previoius time)
    
    # multiply dispersion dimension in each frame by fit_vals
    # gives us intensity calibrated spectra during flare time
    scaled_flare_time = np.zeros(np.shape(for_scale))
    
    if len(nonflare_multfact)>1:
        for i in range(np.shape(for_scale)[0]):
            for k in range(np.shape(for_scale)[2]):
                # intensity calibration factor
                if limbd == 1:
                    scaled_flare_time[i,:,k] = nonflare_multfact[:]*for_scale[i,:,k]/\
                        limbdarkening
                elif limbd == 0:
                    scaled_flare_time[i,:,k] = nonflare_multfact[:]*for_scale[i,:,k]          
    else:
        for i in range(np.shape(for_scale)[0]):
            for k in range(np.shape(for_scale)[2]):
                # intensity calibration factor
                if limbd == 1:
                    scaled_flare_time[i,:,k] = nonflare_multfact*for_scale[i,:,k]/\
                        limbdarkening
                elif limbd == 0:
                    scaled_flare_time[i,:,k] = nonflare_multfact*for_scale[i,:,k]          
            
    
    bkgd_subtract_flaretime = np.zeros(np.shape(for_scale))

    # Subtracted averaged, scaled data cube 
    # from each time step in scaled data cube
    for i in range(np.shape(scaled_flare_time)[0]):
        for j in range(np.shape(scaled_flare_time)[2]):
            if end ==0:
                bkgd_subtract_flaretime[i,:,j] = scaled_flare_time[i,:,j]-nonflare_average[:]
            else:
                bkgd_subtract_flaretime[i,end:,j] = scaled_flare_time[i,end:,j]-nonflare_average[:-end]

        
    return scaled_flare_time, bkgd_subtract_flaretime
    
                            
def pltsubtract(dispersion_range,nonflare_average,scaled_flare_time,muted,indexs,end=5,pid='pid_1_84'):
    
    # plotting routines to compare flare-time with non-flare spectra
    
    fig,ax = plt.subplots(figsize=(10,5))

    if end == 0:
        ax.plot(dispersion_range*10,nonflare_average,\
                color=muted[4],label='Non-Flare')
        ax.plot(dispersion_range*10,scaled_flare_time[0,:,indexs[0]],\
                color=muted[7],label='Flare-Time')
        ax.plot(dispersion_range*10,scaled_flare_time[0,:,indexs[0]]-\
                nonflare_average,color=muted[6],label='Flare-Only')
    else:
        ax.plot(dispersion_range[end:]*10,nonflare_average[:-end],\
                color=muted[4],label='Non-Flare')
        ax.plot(dispersion_range[end:]*10,scaled_flare_time[0,end:,indexs[0]],\
                color=muted[7],label='Flare-Time')
        ax.plot(dispersion_range[end:]*10,scaled_flare_time[0,end:,indexs[0]]-\
                nonflare_average[:-end],color=muted[6],label='Flare-Only')    
    ax.grid()
    #ax.set_ylim([0,5e6])
    ax.legend(loc=0)
    ax.set_title('Non-Flare Estimate vs. Flare-time ',fontsize=25)
    ax.set_xlabel(r'Wavelength [$\mathring A$]',fontsize=15)
    ax.set_ylabel(r'Intensity [$erg/s/cm^2/sr/\mathring A$]',fontsize=15)
    plt.show()
    
    fig.savefig('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/'+pid+\
                '/pltprofile.png')

    return None

def deviations(bkgd_subtract_flaretime,nonflare_average,nonflare_stdevs,end=5):
    
    # in pid_1_84, there is a strange difference between flare-time spectra and
    # nonflare, even outside of the main wings of Ca II H; testing for the
    # flare-time variations in these parts of the spectra and comparing to the
    # variation in the nonflare average in order to test if this variation is
    # due to flaring activity or something else in the spectra
    
    stdevs_flaretime = np.zeros([np.shape(nonflare_average)[0]-5,
                                 np.shape(nonflare_average)[1]])

    for i in range(np.shape(bkgd_subtract_flaretime)[1]-5):
        for j in range(np.shape(bkgd_subtract_flaretime)[2]):
            stdevs_flaretime[i,j] = np.std(bkgd_subtract_flaretime[:,i+end,j])
            
    # Compute PTE in dispersion dimension for all spatial locations
    ptes_flaretime = np.zeros([np.shape(nonflare_average)[0]-5,\
                               np.shape(nonflare_average)[1]]) # just one loc
    totnum = np.shape(ptes_flaretime)[1] # total number of spatial points
    for i in range(np.shape(nonflare_stdevs)[0]-5):
        for j in range(np.shape(nonflare_stdevs)[1]):
            stdev_inquestion_flaretime = stdevs_flaretime[i,j]
            stdev_inquestion_nonflare = nonflare_stdevs[i,:]
            num_gt = numgreaterthan(stdev_inquestion_nonflare,\
                                    stdev_inquestion_flaretime)
            pte = num_gt/totnum
            ptes_flaretime[i,j] = pte

    return stdevs_flaretime, ptes_flaretime


def pltptes(ptes_flaretime,image_data_arr_arr_raster1,pid='pid_1_84'):
    
    # plot probabity to exceed the variation in nonflare spectra at each
    # wavlength in flare-time - high probability to exceed means that variation
    # can be explained by the same variations seen in the quiet sun, low PTE 
    # means that variation can only be explained by something specific to the
    # flare data - either the flare itself, or something specific to that set, 
    # if different from the dataset used to determine QS
    
    fig,[ax,ax1] = plt.subplots(2,1,figsize=(7,12))
    im = ax.pcolormesh(np.log(ptes_flaretime))
    #cbar = plt.colorbar(im, shrink=0.9, pad = 0.05)
    #position=fig1.add_axes([0.93,0.1,0.02,0.35])
    #fig.colorbar(im, orientation="vertical",ax=ax,cax=position)
    ax1.set_xlabel('Spatial',fontsize=13)
    ax.set_ylabel('Dispersion',fontsize=13)
    ax.set_xlabel('Spatial',fontsize=13)
    ax.set_title('P.T.E. standard deviation in non-flare spectrum',fontsize=15)
    ax1.set_title('First Image Frame - Intensity',fontsize=15)
    ax1.pcolormesh(image_data_arr_arr_raster1[0,:,:])
    ax.axvline(1350)
    ax1.axvline(1350)
    ax.axvline(300)
    ax1.axvline(300)
    
    fig.savefig('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/'+pid+\
                '/pltpptes.png')

    return None

#definition of window for continuum

def contwind(sample_flaretime,dispersion_range,maxinds,scaled_flare_time,
             caII_low,caII_high,deg=8,low0=29,high0=29,low1=60,high1=150,
             low2=265,high2=290,low3=360,high3=400,low4=450,high4=480,
             low5=845,high5=870,low6=945,high6=965):
    
    # Define continuum (or "pseudo-continuum") window locations outside of the
    # main lines; this is used to isolate the line, to determine line flux,
    # strength, and for modeling.  In the Ca II window, this "pseudo-continuum"
    # is likely affected by the line itself - a problem to be discussed at a
    # later date
    
    # Default limits are for Ca II H in pid_1_84.  Select six (or fewer? then
    # comment extra lines) windows to fit a polynomial to as an estimate of
    # the pseudo-continuum with absorption (telluric or otherwise) and emission
    # lines removed
    
    avgs = []
    for i in range(len(scaled_flare_time)):
        snapshot = scaled_flare_time[i,:,:]
        average = np.mean(snapshot,axis=0)
        avgs.append(average)
    
    contwind0_1 = sample_flaretime[low0:high0]
    contwind0_1_wave = dispersion_range[low0:high0]
    contwind1 = np.mean(sample_flaretime[low1:high1])
    contwind1_wave = np.mean(dispersion_range[low1:high1])
    contwind2 = np.mean(sample_flaretime[low2:high2])
    contwind2_wave = np.mean(dispersion_range[low2:high2])
    contwind3 = np.mean(sample_flaretime[low3:high3])
    contwind3_wave = np.mean(dispersion_range[low3:high3])
    contwind4 = np.mean(sample_flaretime[low4:high4])
    contwind4_wave = np.mean(dispersion_range[low4:high4])
    contwind5 = np.mean(sample_flaretime[low5:high5])
    contwind5_wave = np.mean(dispersion_range[low5:high5])
    contwind6 = np.mean(sample_flaretime[low6:high6])
    contwind6_wave = np.mean(dispersion_range[low6:high6])
    
    cont_int_array = [contwind0_1,contwind1,contwind2,contwind3,contwind4,
                      contwind5,contwind6]
    cont_int_wave_array = [contwind0_1_wave,contwind1_wave,contwind2_wave,
                           contwind3_wave,contwind4_wave,contwind5_wave,contwind6_wave]
    
    # polynomial fit of degree deg; deg = 8 likely oversolves?
    p = np.poly1d(np.polyfit(cont_int_wave_array,cont_int_array,deg))
    nolines = p(dispersion_range)
    
    return nolines, cont_int_array, cont_int_wave_array

def widths_strengths(ew_CaII_all_fs,eqw_CaII_all_fs,width_CaII_all_fs,
                     ew_hep_all_fs,eqw_hep_all_fs,width_hep_all_fs,
                     caII_low,caII_high,hep_low,hep_high,
                     scaled_flare_time,bkgd_subtract_flaretime,
                     dispersion_range,deg=6,low0=29,high0=30,low1=60,high1=150,
                     low2=265,high2=290,low3=360,high3=400,low4=450,high4=480,
                     low5=845,high5=870,low6=945,high6=965):
    
    # determine equivanet widths, effective widths, line widths for Ca II line
    
    # Uses pseudo-continuum polynomial determination similar to that described
    # in function above; see that doc for description
    
    avgs = []
    for i in range(len(scaled_flare_time)):
        snapshot = scaled_flare_time[i,:,:]
        average = np.mean(snapshot,axis=0)
        avgs.append(average)
    
    # "eq" means for use in equivalent width determination - equivalent width
    # determination does not use the background-subtracted values, but the raw
    # intensity-calibrated spectra (see description of effective vs. equivalent
    # with for justification)
    for j in range(np.shape(bkgd_subtract_flaretime)[2]):
        for i in range(len(scaled_flare_time)-5):
            print(i)
            print(j)
            sample_flaretime = bkgd_subtract_flaretime[i,:,j]
            foreqw = scaled_flare_time[i,:,j]
            contwind0_1 = np.nanmean(sample_flaretime[low0:high0])
            contwind0_1eq = np.nanmean(foreqw[low0:high0])
            contwind0_1_wave = np.nanmean(dispersion_range[low0:high0])
            contwind1eq = np.mean(foreqw[low1:high1])
            contwind1 = np.mean(sample_flaretime[low1:high1])
            contwind1_wave = np.mean(dispersion_range[low1:high1])
            contwind2eq = np.mean(foreqw[low2:high2])
            contwind2 = np.mean(sample_flaretime[low2:high2])
            contwind2_wave = np.mean(dispersion_range[low2:high2])
            contwind3eq = np.mean(foreqw[low3:high3])
            contwind3 = np.mean(sample_flaretime[low3:high3])
            contwind3_wave = np.mean(dispersion_range[low3:high3])
            contwind4eq = np.mean(foreqw[low4:high4])
            contwind4 = np.mean(sample_flaretime[low4:high4])
            contwind4_wave = np.mean(dispersion_range[low4:high4])
            contwind5eq = np.mean(foreqw[low5:high5])
            contwind5 = np.mean(sample_flaretime[low5:high5])
            contwind5_wave = np.mean(dispersion_range[low5:high5])
            contwind6eq = np.mean(foreqw[low6:high6])
            contwind6 = np.mean(sample_flaretime[low6:high6])
            contwind6_wave = np.mean(dispersion_range[low6:high6])
    
            cont_int_arrayeqw = [contwind0_1eq,contwind1eq,contwind2eq,
                                 contwind3eq,contwind4eq,contwind5eq,contwind6eq]
    
            cont_int_array = [contwind0_1,contwind1,contwind2,contwind3,
                              contwind4,contwind5,contwind6]
            cont_int_wave_array = [contwind0_1_wave,contwind1_wave,
                                   contwind2_wave,contwind3_wave,
                                   contwind4_wave,contwind5_wave,
                                   contwind6_wave]
    
            deg = 6
    
            p = np.poly1d(np.polyfit(cont_int_wave_array,cont_int_array,deg))
    
            peq = np.poly1d(np.polyfit(cont_int_wave_array,cont_int_arrayeqw,
                                       deg))
    
            nolines = p(dispersion_range)
            nolineseq = peq(dispersion_range)
    
            maxind = np.argmax(sample_flaretime)
            maxint = sample_flaretime[maxind]
            maxcont = nolines[maxind]
    
            integrand = (sample_flaretime-nolines)/(maxint-maxcont)
            normflux = np.divide(foreqw,nolineseq)
            
            integrand2 = 1 - normflux
    
            #equivalent width determination, efefctive width determinatino
            ew_caII = integrate.trapezoid(integrand[caII_low:caII_high],
                                         dispersion_range[caII_low:caII_high])
                
            eqw_caII = integrate.trapezoid(integrand2[caII_low:caII_high],
                                          dispersion_range[caII_low:caII_high])
            
            maxind_Hep = np.argmax(sample_flaretime[hep_low:hep_high])
            maxint_Hep = sample_flaretime[maxind_Hep+hep_low]
            maxcont_Hep = nolines[maxind_Hep+hep_low]
    
            integrand_Hep = (sample_flaretime[hep_low:hep_high]-
                             nolines[hep_low:hep_high])/\
                (maxint_Hep-maxcont_Hep)
    
            ew_Hep = integrate.trapezoid(integrand_Hep,
                                        dispersion_range[hep_low:hep_high])
            eqw_Hep = integrate.trapezoid(integrand2[hep_low:hep_high],
                                         dispersion_range[hep_low:hep_high])
    
            ew_CaII_all_fs[i,j]=ew_caII
            ew_hep_all_fs[i,j]=ew_Hep
            eqw_CaII_all_fs[i,j]=eqw_caII
            eqw_hep_all_fs[i,j]=eqw_Hep
            
            caII_isolate = sample_flaretime[caII_low:caII_high]
            mincaIIH = min(caII_isolate)
            maxcaIIH = max(caII_isolate)
            meancaIIH = (maxcaIIH+mincaIIH)/2
    
            caIIHmidlow, caIImidlowindex = \
                find_nearest(caII_isolate[:round(len(caII_isolate)/2)],meancaIIH)
            caIIHmidhigh,caIImidhighindex = \
                find_nearest(caII_isolate[round(len(caII_isolate)/2):],meancaIIH)
    
            widthAng_caII = dispersion_range[caII_low+\
                                             round(len(caII_isolate)/2)+\
                                                 caIImidhighindex-1]-\
                dispersion_range[caII_low+caIImidlowindex-1] 
    
            hep_isolate = sample_flaretime[hep_low:hep_high]
            minhep = min(hep_isolate)
            maxhep = max(hep_isolate)
            meanhep = (maxhep+minhep)/2
    
            hepmidlow,hepmidlowindex = \
                find_nearest(hep_isolate[:round(len(hep_isolate)/2)],meanhep)
            hepmidhigh,hepmidhighindex = \
                find_nearest(hep_isolate[round(len(hep_isolate)/2):],meanhep)
    
            widthAng_hep = dispersion_range[hep_low+round(len(hep_isolate)/2)+\
                                            hepmidhighindex-1]-\
                dispersion_range[hep_low+hepmidlowindex-1] 
    
            width_CaII_all_fs[i,j]=widthAng_caII
            width_hep_all_fs[i,j]=widthAng_hep
        
    return ew_CaII_all_fs, ew_hep_all_fs, eqw_CaII_all_fs, eqw_hep_all_fs,\
        width_CaII_all_fs, width_hep_all_fs
        
def widths_strengths_oneline(ew_line_all_fs,eqw_line_all_fs,width_line_all_fs,
                             int_line_all_fs,line_low,line_high,scaled_flare_time,
                             bkgd_subtract_flaretime,dispersion_range,maxinds,
                             low0=29,high0=29,low1=60,high1=150,low2=265,
                             high2=290,low3=360,high3=400,low4=450,high4=480,
                             low5=845,high5=870,low6=945,high6=965,end=0,deg=7,
                             alt=0,
                             altinds = [380, 390, 400, 410, 450, 480, 647, 700,\
                                        820, 850, 900]):
    
    # determine equivanet widths, effective widths, line widths for single line
    
    # Uses pseudo-continuum polynomial determination similar to that described
    # in function above; see that doc for description
    
    # deg = 7 is best for pid_1_38, deg = 6 for pid_1_84
    
    avgs = []
    for i in range(len(scaled_flare_time)):
        snapshot = scaled_flare_time[i,:,:]
        average = np.mean(snapshot,axis=0)
        avgs.append(average)
    
    # "eq" means for use in equivalent width determination - equivalent width
    # determination does not use the background-subtracted values, but the raw
    # intensity-calibrated spectra (see description of effective vs. equivalent
    # with for justification)
    for j in range(np.shape(bkgd_subtract_flaretime)[2]):
        for i in range(len(scaled_flare_time)-end-1):
            
            # two methods - averaging continuum windows (pid_1_84) is ind = 0
            if alt == 0:
                sample_flaretime = bkgd_subtract_flaretime[i,:,j]
                foreqw = scaled_flare_time[i,:,j]
                contwind0_1 = np.mean(sample_flaretime[low0:high0])
                contwind0_1eq = np.mean(foreqw[low0:high0])
                contwind0_1_wave = np.mean(dispersion_range[low0:high0])
                contwind1eq = np.mean(foreqw[low1:high1])
                contwind1 = np.mean(sample_flaretime[low1:high1])
                contwind1_wave = np.mean(dispersion_range[low1:high1])
                contwind2eq = np.mean(foreqw[low2:high2])
                contwind2 = np.mean(sample_flaretime[low2:high2])
                contwind2_wave = np.mean(dispersion_range[low2:high2])
                contwind3eq = np.mean(foreqw[low3:high3])
                contwind3 = np.mean(sample_flaretime[low3:high3])
                contwind3_wave = np.mean(dispersion_range[low3:high3])
                contwind4eq = np.mean(foreqw[low4:high4])
                contwind4 = np.mean(sample_flaretime[low4:high4])
                contwind4_wave = np.mean(dispersion_range[low4:high4])
                contwind5eq = np.mean(foreqw[low5:high5])
                contwind5 = np.mean(sample_flaretime[low5:high5])
                contwind5_wave = np.mean(dispersion_range[low5:high5])
                contwind6eq = np.mean(foreqw[low6:high6])
                contwind6 = np.mean(sample_flaretime[low6:high6])
                contwind6_wave = np.mean(dispersion_range[low6:high6])
        
                cont_int_arrayeqw = [contwind0_1eq,contwind1eq,contwind2eq,
                                     contwind3eq,contwind4eq,contwind5eq,contwind6eq]
        
                cont_int_array = [contwind0_1,contwind1,contwind2,contwind3,
                                  contwind4,contwind5,contwind6]
                cont_int_wave_array = [contwind0_1_wave,contwind1_wave,
                                       contwind2_wave,contwind3_wave,
                                       contwind4_wave,contwind5_wave,
                                       contwind6_wave]
        
        
                p = np.poly1d(np.polyfit(cont_int_wave_array,cont_int_array,deg))
        
                peq = np.poly1d(np.polyfit(cont_int_wave_array,cont_int_arrayeqw,
                                           deg))
        
                nolines = p(dispersion_range)
                nolineseq = peq(dispersion_range)
        
                maxind = np.argmax(sample_flaretime)
                maxint = sample_flaretime[maxind]
                maxcont = nolines[maxind]
                selwl = dispersion_range[line_low:line_high]
                sel = bkgd_subtract_flaretime[i,line_low:line_high,j]-\
                    nolines[line_low:line_high]
                
            #second is just taking points from the continuum, in altinds
            elif alt == 1: 
                sample_flaretime = bkgd_subtract_flaretime[i,:,j]
                foreqw = scaled_flare_time[i,:,j]
                cont_int_array = bkgd_subtract_flaretime[i,altinds,j]
                cont_int_wave_array = [dispersion_range[i] for i in altinds]
                cont_int_arrayeqw = foreqw[altinds]
                p = np.poly1d(np.polyfit(cont_int_wave_array,cont_int_array,deg))
                peq = np.poly1d(np.polyfit(cont_int_wave_array,cont_int_arrayeqw,deg))
                nolines = p(dispersion_range)
                nolineseq = peq(dispersion_range)
                maxind = np.argmax(sample_flaretime)
                maxint = sample_flaretime[maxind]
                maxcont = nolines[maxind]
                selwl = dispersion_range[line_low:line_high]
                sel = bkgd_subtract_flaretime[i,line_low:line_high,j]-\
                    nolines[line_low:line_high]
                
    
            integrand = (sample_flaretime-nolines)/(maxint-maxcont)
            normflux = np.divide(foreqw,nolineseq)
    
            integrand2 = 1 - normflux
    
            #equivalent width determination, efefctive width determinatino
            ew_line = integrate.trapezoid(integrand[line_low:line_high],
                                         dispersion_range[line_low:line_high])\
                [-1]
            eqw_line = integrate.trapezoid(integrand2[line_low:line_high],
                                          dispersion_range[line_low:line_high])\
                [-1]
            int_line = integrate.trapezoid(sel,selwl)[-1]*10
            print(int_line)
            print(i)
            print(j)

            ew_line_all_fs[i,j]=ew_line
            eqw_line_all_fs[i,j]=eqw_line
            int_line_all_fs[i,j]=int_line
            
            line_isolate = sample_flaretime[line_low:line_high]
            minline = min(line_isolate)
            maxline = max(line_isolate)
            meanline = (maxline+minline)/2
    
            linemidlow, linemidlowindex = \
                find_nearest(line_isolate[:round(len(line_isolate)/2)],meanline)
            linemidhigh,linemidhighindex = \
                find_nearest(line_isolate[round(len(line_isolate)/2):],meanline)
    
            widthAng_line = dispersion_range[line_low+\
                                             round(len(line_isolate)/2)+\
                                                 linemidhighindex-1]-\
                dispersion_range[line_low+linemidlowindex-1] 
    
            width_line_all_fs[i,j]=widthAng_line
        
    return ew_line_all_fs, eqw_line_all_fs, width_line_all_fs, int_line_all_fs

# NOTE: Add plotting routine for widths?
def plt_line_characteristics(ew_line_all_fs,eqw_line_all_fs,width_line_all_fs,
                             int_line_all_fs,ew_line2_all_fs,eqw_line2_all_fs,width_line2_all_fs,
                             int_line2_all_fs,maxindices,times,muted,pid='pid_1_84',
                             nslitpos=4,raster=0,nframes=11,line='Ca II H',line2='H-Epsilon',end=16):
    
    # defined quantites are all for pid_1_84
    # raster is start position to begin tracking
    # nslitpos is the number of slit positions in a scan
    # nframes is number of raster scans loaded in the array 

    # arrays to populate with specific kernel
    kernindews = []
    kernindeqws = []
    kernindwidths = []
    kernindflxs = []
    
    # arrays to populate with specific kernel
    kernindews2 = []
    kernindeqws2 = []
    kernindwidths2 = []
    kernindflxs2 = []
    
    timeshhmmss = []
    
    timesdt = []
    
    for i in times:
        timesdt.append(datetime.strptime(i, '%Y-%m-%dT%H:%M:%S.%f'))
        
    timesdthr = []
    
    for i in timesdt:
        timesdthr.append(datetime.strftime(i, 
                                     "%H:%M:%S"))
    
    # times in correct format for plotting
    for i in range(len(times)):
        timeshhmmss.append(times[i][-12:-4])

   
    # append correct value depending on kernel location defined by input array
    for i in range(len(times)):
        kernindews.append(ew_line_all_fs[i,maxindices[i]])
        kernindeqws.append(eqw_line_all_fs[i,maxindices[i]])
        kernindwidths.append(width_line_all_fs[i,maxindices[i]])
        kernindflxs.append(int_line_all_fs[i,maxindices[i]])
        
        kernindews2.append(ew_line2_all_fs[i,maxindices[i]])
        kernindeqws2.append(eqw_line2_all_fs[i,maxindices[i]])
        kernindwidths2.append(width_line2_all_fs[i,maxindices[i]])
        kernindflxs2.append(int_line2_all_fs[i,maxindices[i]])
    
    
    print('length = '+str(len(kernindews)))
    # plotting routines
    #print(len(timeshhmmss))
    fig,((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8))=plt.subplots(4,2,figsize=(10,7))
    #fig.suptitle(line+' Line Characteristics, '+pid+'',fontsize=25)
    lns1 = ax1.plot(timesdt[0:-end:nslitpos],kernindews[:-end:nslitpos],'-o',color=muted[0],label='Ca II H')
    lns2 = ax2.plot(timesdt[1:-end:nslitpos],kernindews[1:-end:nslitpos],'-o',color=muted[1],label='Ca II H')
    lns3 = ax3.plot(timesdt[2:-end:nslitpos],kernindews[2:-end:nslitpos],'-o',color=muted[2],label='Ca II H')
    lns4 = ax4.plot(timesdt[3:-end:nslitpos],kernindews[3:-end:nslitpos],'-o',color=muted[3],label='Ca II H')
    
    lns5 = ax5.plot(timesdt[0:-end:nslitpos],kernindeqws[:-end:nslitpos],'-o',color=muted[0],label='Ca II H')
    lns6 = ax6.plot(timesdt[1:-end:nslitpos],kernindeqws[1:-end:nslitpos],'-o',color=muted[1],label='Ca II H')
    lns7 = ax7.plot(timesdt[2:-end:nslitpos],kernindeqws[2:-end:nslitpos],'-o',color=muted[2],label='Ca II H')
    lns8 = ax8.plot(timesdt[3:-end:nslitpos],kernindeqws[3:-end:nslitpos],'-o',color=muted[3],label='Ca II H')
    
    ax1_0 = ax1.twinx()
    ax2_0 = ax2.twinx()
    ax3_0 = ax3.twinx()
    ax4_0 = ax4.twinx()
    ax5_0 = ax5.twinx()
    ax6_0 = ax6.twinx()
    ax7_0 = ax7.twinx()
    ax8_0 = ax8.twinx()
    
    lns1_0 = ax1_0.plot(timesdt[0:-end:nslitpos],kernindews2[:-end:nslitpos],'-o',color=muted[4],label=r'$H\epsilon$')
    lns2_0 = ax2_0.plot(timesdt[1:-end:nslitpos],kernindews2[1:-end:nslitpos],'-o',color=muted[5],label=r'$H\epsilon$')
    lns3_0 =ax3_0.plot(timesdt[2:-end:nslitpos],kernindews2[2:-end:nslitpos],'-o',color=muted[6],label=r'$H\epsilon$')
    lns4_0 =ax4_0.plot(timesdt[3:-end:nslitpos],kernindews2[3:-end:nslitpos],'-o',color=muted[7],label=r'$H\epsilon$')
    
    lns5_0 =ax5_0.plot(timesdt[0:-end:nslitpos],kernindeqws2[:-end:nslitpos],'-o',color=muted[4],label=r'$H\epsilon$')
    lns6_0 =ax6_0.plot(timesdt[1:-end:nslitpos],kernindeqws2[1:-end:nslitpos],'-o',color=muted[5],label=r'$H\epsilon$')
    lns7_0 =ax7_0.plot(timesdt[2:-end:nslitpos],kernindeqws2[2:-end:nslitpos],'-o',color=muted[6],label=r'$H\epsilon$')
    lns8_0 =ax8_0.plot(timesdt[3:-end:nslitpos],kernindeqws2[3:-end:nslitpos],'-o',color=muted[7],label=r'$H\epsilon$')
        
    
    
    ax1.set_ylabel(r'$\Delta\lambda_{eff}$ [nm], Ca II H',fontsize=10,color=muted[0])
    ax1_0.set_ylabel(r'$\Delta\lambda_{eff}$ [nm], H$\epsilon$',fontsize=10,color=muted[4])
    ax2.set_ylabel(r'$\Delta\lambda_{eff}$ [nm], Ca II H',fontsize=10,color=muted[1])
    ax2_0.set_ylabel(r'$\Delta\lambda_{eff}$ [nm], H$\epsilon$',fontsize=10,color=muted[5])
    ax3.set_ylabel(r'$\Delta\lambda_{eff}$ [nm], Ca II H',fontsize=10,color=muted[2])
    ax3_0.set_ylabel(r'$\Delta\lambda_{eff}$ [nm], H$\epsilon$',fontsize=10,color=muted[6])
    ax4.set_ylabel(r'$\Delta\lambda_{eff}$ [nm], Ca II H',fontsize=10,color=muted[3])
    ax4_0.set_ylabel(r'$\Delta\lambda_{eff}$ [nm], H$\epsilon$',fontsize=10,color=muted[7])
    
    ax5.set_ylabel(r'$\Delta\lambda_{eq}$ [nm], Ca II H',fontsize=10,color=muted[0])
    ax5_0.set_ylabel(r'$\Delta\lambda_{eq}$ [nm], H$\epsilon$',fontsize=10,color=muted[4])
    ax6.set_ylabel(r'$\Delta\lambda_{eq}$ [nm], Ca II H',fontsize=10,color=muted[1])
    ax6_0.set_ylabel(r'$\Delta\lambda_{eq}$ [nm], H$\epsilon$',fontsize=10,color=muted[5])
    ax7.set_ylabel(r'$\Delta\lambda_{eq}$ [nm], Ca II H',fontsize=10,color=muted[2])
    ax7_0.set_ylabel(r'$\Delta\lambda_{eq}$ [nm], H$\epsilon$',fontsize=10,color=muted[6])
    ax8.set_ylabel(r'$\Delta\lambda_{eq}$ [nm], Ca II H',fontsize=10,color=muted[3])
    ax8_0.set_ylabel(r'$\Delta\lambda_{eq}$ [nm], H$\epsilon$',fontsize=10,color=muted[7])
    
    ax1.set_xticks(timesdt[0:-end:4],[])
    ax2.set_xticks(timesdt[0:-end:4],[])
    ax3.set_xticks(timesdt[0:-end:4],[])
    ax4.set_xticks(timesdt[0:-end:4],[])
    ax5.set_xticks(timesdt[0:-end:4],[])
    ax6.set_xticks(timesdt[0:-end:4],[])
    ax7.set_xticks(timesdt[0:-end:4],timeshhmmss[0:-end:4],rotation=45)
    ax8.set_xticks(timesdt[0:-end:4],timeshhmmss[0:-end:4],rotation=45)

    
    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax4.grid()
    ax5.grid()
    ax6.grid()
    ax7.grid()
    ax8.grid()
    timesdt0 = datetime(2022, 8, 19, 20, 42, 3, 5)
    
    ax1.set_xlim([timesdt0,timesdt[-end-1]])
    ax2.set_xlim([timesdt0,timesdt[-end-1]])
    ax3.set_xlim([timesdt0,timesdt[-end-1]])
    ax4.set_xlim([timesdt0,timesdt[-end-1]])
    ax5.set_xlim([timesdt0,timesdt[-end-1]])
    ax6.set_xlim([timesdt0,timesdt[-end-1]])
    ax7.set_xlim([timesdt0,timesdt[-end-1]])
    ax8.set_xlim([timesdt0,timesdt[-end-1]])
    
    ax1_0.set_xlim([timesdt0,timesdt[-end-1]])
    ax2_0.set_xlim([timesdt0,timesdt[-end-1]])
    ax3_0.set_xlim([timesdt0,timesdt[-end-1]])
    ax4_0.set_xlim([timesdt0,timesdt[-end-1]])
    ax5_0.set_xlim([timesdt0,timesdt[-end-1]])
    ax6_0.set_xlim([timesdt0,timesdt[-end-1]])
    ax7_0.set_xlim([timesdt0,timesdt[-end-1]])
    ax8_0.set_xlim([timesdt0,timesdt[-end-1]])
    
    #ax1.set_ylim([])
    # ax2.set_ylim([timesdt0,timesdt[-end-1]])
    ax3.set_ylim([.0455,0.050])
    #ax4.set_ylim([timesdt0,timesdt[-end-1]])
    ax5.set_ylim([-0.22,-0.07])
    # ax6.set_ylim([timesdt0,timesdt[-end-1]])
    # ax7.set_ylim([timesdt0,timesdt[-end-1]])
    #ax8.set_ylim([timesdt0,timesdt[-end-1]])
    
    ax1_0.set_ylim([.006 ,.011])
    #ax2_0.set_ylim([timesdt0,timesdt[-end-1]])
    #ax3_0.set_ylim([timesdt0,timesdt[-end-1]])
    ax4_0.set_ylim([.0085,0.0105])
    ax5_0.set_ylim([-0.03,-0.004])
    #ax6_0.set_ylim([timesdt0,timesdt[-end-1]])
    #ax7_0.set_ylim([timesdt0,timesdt[-end-1]])
    #ax8_0.set_ylim([timesdt0,timesdt[-end-1]])

    lns = lns1+lns1_0
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs)
    

    lns = lns2+lns2_0
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs)


    lns = lns3+lns3_0
    labs = [l.get_label() for l in lns]
    ax3.legend(lns, labs)
    
    
    lns = lns4+lns4_0
    labs = [l.get_label() for l in lns]
    ax4.legend(lns, labs)
    
    
    lns = lns5+lns5_0
    labs = [l.get_label() for l in lns]
    ax5.legend(lns, labs)
    
    
    lns = lns6+lns6_0
    labs = [l.get_label() for l in lns]
    ax6.legend(lns, labs)
    
    
    lns = lns7+lns7_0
    labs = [l.get_label() for l in lns]
    ax7.legend(lns, labs)
    
    
    lns = lns8+lns8_0
    labs = [l.get_label() for l in lns]
    ax8.legend(lns, labs)

    
    # ax1.get_xaxis().set_visible(False)
    # ax1_0.get_xaxis().set_visible(False)

    # ax2.get_xaxis().set_visible(False)
    # ax2_0.get_xaxis().set_visible(False)
    
    # ax3.get_xaxis().set_visible(False)
    # ax3_0.get_xaxis().set_visible(False)
    
    # ax4.get_xaxis().set_visible(False)
    # ax4_0.get_xaxis().set_visible(False)
    
    # ax5.get_xaxis().set_visible(False)
    # ax5_0.get_xaxis().set_visible(False)
    
    # ax6.get_xaxis().set_visible(False)
    # ax6_0.get_xaxis().set_visible(False)
    

    fig.tight_layout(h_pad=0,w_pad=4)
    
    plt.show()
    
    #fig.savefig('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/'+pid+\
    #            '/linecharacteristics.png')
    
    return None
    
    

def gauss2fit(storeamp1,storemu1,storesig1,storeamp2,storemu2,storesig2,
              bkgd_subtract_flaretime,dispersion_range, double_gaussian_fit,
              maxinds,times_raster1,caII_low,caII_high,double_gaussian,gaussian,selwl,sel,
              pid='pid_1_84',parameters = [2e6,396.82,0.01,2e6,396.86,0.015]):
    fig, ax = plt.subplots(3,4,figsize=(30,30))


    
    # Original script for double-Gaussian fitting and plotting; no error metrics,
    # just visualization, limited room for model selection.  Not for
    # current use; for better option, look to functions "fittingroutines"
    # and "pltfitresults"

    fig.suptitle('Ca II H line evolution, 19-Aug-2022, Raster 1, Kernel Center'
                 ,fontsize=35)
    
    for i in range(11):
        selwl = dispersion_range[caII_low:caII_high]
        sel = bkgd_subtract_flaretime[i,caII_low:caII_high,maxinds[i]]-\
            min(bkgd_subtract_flaretime[i,caII_low:caII_high,maxinds[i]])
        fit = leastsq(double_gaussian_fit,parameters,(selwl,sel))
        [c1,mu1,sigma1,c2,mu2,sigma2] = fit[0]
        print(fit[0][0])
        storeamp1.append(c1)
        storemu1.append(mu1)
        storesig1.append(sigma1)
        storeamp2.append(c2)
        storemu2.append(mu2)
        storesig2.append(sigma2)
        ax.flatten()[i].plot(selwl,sel)
        ax.flatten()[i].plot(selwl, double_gaussian( selwl,fit[0][0],fit[0][1],fit[0][2],fit[0][3],fit[0][4],fit[0][5]))
        ax.flatten()[i].plot(selwl,gaussian(selwl,fit[0][0],fit[0][1],fit[0][2]),c='g')
        ax.flatten()[i].plot(selwl,gaussian(selwl,fit[0][3],fit[0][4],fit[0][5]),c='g')
        ax.flatten()[i].axis(ymin=0,ymax=7.5e6)
        ax.flatten()[i].set_title(times_raster1[i],fontsize=30)
        ax.flatten()[i].grid()
    ax.flatten()[-1].axis('off')
    
    fig.savefig('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/'+pid+\
                '/gaussianfits.png')
    
    #save
    
    return storeamp1, storeamp2, storesig1, storesig2, storemu1, storemu2

def fittingroutines(bkgd_subtract_flaretime,dispersion_range,
                    times_raster1, line_low, line_high,
                    double_gaussian, gaussian, selwl,sel,paramsgauss,
                    params2gauss,params2gaussneg,maxinds,storeamp1,storemu1,
                    storesig1,storeamp2,storemu2,storesig2,pid='pid_1_84',
                    date = '08/09/2022',line = 'Ca II H',nimg = 7,
                    inds=[410,460,510,520,590,600,610,620,630,640,650,660,670,680,690,700,
                          720],deg=7,start=0):
    # More flexible line fitting routines; currently for Ca II H as observed
    # in pid_1_84, but flexible for application to other lines.  Currently also
    # only includes functinoality for single Gaussian and double Gaussian fits;
    # additions welcome (skewed Gauss? Lorentz? Voigt? etc.)
    
    # Returns, via scipy.optimize.curve_fit, both the fit parameters and the 
    # error metrics for each model
    
    # defaults currently for pid_1_38
    
    fits_1g = []
    fits_2g = []
    fits_2gneg = []
    
    l=0
    j=0
    stopind=0
    # first iteration. Use this as guide for the fit parameters to follow
    for j in range(nimg):
        kernind = maxinds[start+j]
        if l == 0:
             
            cont_int_array = bkgd_subtract_flaretime[start+j,inds,kernind]
            cont_int_wave_array = dispersion_range[inds]
            deg = deg
            p = np.poly1d(np.polyfit(cont_int_wave_array,cont_int_array,deg))
            nolines = p(dispersion_range)
            selwl = dispersion_range[line_low:line_high]
            sel = bkgd_subtract_flaretime[start+j,line_low:line_high,kernind]-\
                nolines[line_low:line_high]
                
            
            
            if j == 0:
                fig,ax=plt.subplots()
                ax.plot(selwl,sel)
                ax.plot(dispersion_range,nolines)
                ax.plot(dispersion_range,bkgd_subtract_flaretime[start+j,:,kernind])
                ax.set_xlim([dispersion_range[0],dispersion_range[-1]])
                plt.show()
                
                
            try:
                
                #initial attempt at a fit
                fit2g, fit2gcov = curve_fit(double_gaussian,selwl,sel, 
                                            p0=params2gauss,
                                            maxfev=5000)
                
                if fit2g[0]/fit2g[3] > 1 or np.abs(fit2g[4]-fit2g[1])>0.04:
                    continue
                else:
                    # appropriate fit found!
                    l=1
                
            except RuntimeError:
                continue
        elif l==1:
            # ID the first time the fit was good
            stopind = start+j
            break
        
            
    
    # redefine gaussian parameters based on result of initial fit
    params2gauss = fit2g
    params2gaussnew = params2gauss
    paramsgauss = paramsgauss # should be fine based on initial gauss, with 1g
    
    
    # print(params2gauss)
    for i in range(nimg):
        
        kernind = maxinds[i+start]
        cont_int_array = bkgd_subtract_flaretime[i+start,inds,kernind]
        cont_int_wave_array = dispersion_range[inds]
        deg =deg
        p = np.poly1d(np.polyfit(cont_int_wave_array,cont_int_array,deg))
        nolines = p(dispersion_range)
        selwl = dispersion_range[line_low:line_high]
        sel = bkgd_subtract_flaretime[i+start,line_low:line_high,kernind]-\
            nolines[line_low:line_high]
                
        try:
            
            fit1g, fit1gcov = curve_fit(gaussian,selwl,sel,p0=paramsgauss)
            fit2g, fit2gcov = curve_fit(double_gaussian,selwl,sel, 
                                        p0=params2gaussnew,
                                        maxfev=10000)
            
            #not iterative
            fit2g_neg, fit2gcov_neg = curve_fit(double_gaussian,selwl,sel,
                                                p0 = params2gaussneg,
                                                maxfev=10000)
            #fit2gneg, fit2gnegcov = curve_fit(double_gaussian,selwl,\ 
                #sel,p0=params2gaussneg,maxfev=5000)
                
            fits_1g.append([fit1g,fit1gcov])
            fits_2g.append([fit2g,fit2gcov])
            fits_2gneg.append([fit2g_neg,fit2gcov_neg])
            
            params2gauss = fit2g
            storeamp1.append(params2gauss[0])
            storemu1.append(params2gauss[1])
            storesig1.append(params2gauss[2])
            storeamp2.append(params2gauss[3])
            storemu2.append(params2gauss[4])
            storesig2.append(params2gauss[5])
            
            
        except RuntimeError:
            print('Runtime Error! Worth a check at index '+str(i))
            fits_1g.append([fit1g,fit1gcov])
            fits_2g.append([fit2g,fit2gcov])
            fits_2gneg.append(['NaN','NaN'])
            continue
        
        #fits_2gneg.append([fit2gneg,fit2gnegcov])

    return fits_1g, fits_2g, fits_2gneg, params2gaussnew,stopind,storeamp1,\
        storemu1,storesig1,storeamp2,storemu2,storesig2,selwl,sel
    

def pltfitresults(bkgd_subtract_flaretime,dispersion_range,double_gaussian,
                  gaussian,times,muted,
                  line_low,line_high,fits_1g,fits_2g,fits_2gneg,maxinds,mu,
                  pid='pid_1_84',
                  date = '08092022',line = 'Ca II H',nimg = 7,nrow=2,ncol=4,
                  lamb0 = 396.847,c=2.99e5,
                  note='',lim=0.3e6,
                  yhigh=1.5e6,inds=[410,460,510,520,590,600,610,620,630,640,650,
                                    660,670,680,690,700,720],deg=7,start=0):
    
    # plotting of the output of "fittingroutines"; can expand to beyond first
    # few image frames.  Tested 1 Dec 2023 for pid_1_84 Ca II H but not beyond
    # this.
    # "inds" is an array containing values for continuum or pseudo-continuum 
    # locations in the images provided from DKIST.  Will vary based on line in
    # question; default is for pid_1_84
    
    # define path
    path = '/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/'+\
                pid+'/'
    if os.path.isdir(path) == False:
        os.mkdir(path)
        
    # initialize plot
    fig, ax = plt.subplots(nrow,ncol,figsize=(30,15))
    fig.suptitle(line+' evolution w/ fitting, '+date+note,fontsize=20)
    
    # select wavelengths of interest
    selwl = dispersion_range[line_low:line_high]
    
    # change type
    if type(selwl) == list:
        selwl = np.array(selwl)
    
    # different x-axes (doppler shift, for example, or deviation from central)
    selwlshift = selwl-lamb0
    
    selwlvel = (selwl/lamb0-1)*c
    
    def veltrans(x):
        return ((((x+lamb0)/lamb0)-1)*c)/mu
        #return ((((x+lamb0)/lamb0)-1)*c) #old, no mu angle difference
    
    def wltrans(x):
        return ((((x/c)+1)*lamb0)-lamb0)
    
    for i in range(nimg):
        
        # define continuum for subtraction, similar to fitting routines above
        kernind = maxinds[i+start]
        cont_int_array = bkgd_subtract_flaretime[i+start,inds,kernind]
        cont_int_wave_array = dispersion_range[inds]
        deg = deg
        p = np.poly1d(np.polyfit(cont_int_wave_array,cont_int_array,deg))
        nolines = p(dispersion_range)
        selwl = dispersion_range[line_low:line_high]
        sel = bkgd_subtract_flaretime[i+start,line_low:line_high,kernind]-\
            nolines[line_low:line_high]
        
        if i == 0:
            maxprofile = max(sel)
        else:
            maxprofilenow = max(sel)
            if maxprofilenow > maxprofile:
                maxprofile = maxprofilenow
            
        # fits in question
        fit1g = fits_1g[i][0]
        fit2g = fits_2g[i][0]
        #fit2gneg = fits_2gneg[i][0]
        
        # plotting routines
        if i > 0 and i % (nrow*ncol) == 0:
            secaxx = ax.flatten()[0].secondary_xaxis('top', functions=(veltrans,
                                                                       wltrans))
            secaxx.set_xlabel(r'Velocity $[km\; s^{-1}]$')
            ax.flatten()[0].set_xlabel(r' $\lambda$ - $\lambda_0$ [nm]')
            ax.flatten()[0].set_ylabel(r'Intensity (- $I_{min}$) $[W\; cm^{-2} sr^{-1} \AA^{-1}]$')


            
            #plt.subplots_adjust(wspace=0.4,hspace=1.0)
            plt.tight_layout(pad = 4)
            fig.savefig('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/'+\
                        pid+'/fits'+date+note+'_endatstep_'+str(i-1)+'.png')
            fig, ax = plt.subplots(nrow,ncol,figsize=(30,15))
            fig.suptitle(line+' evolution w/ fitting, '+date+note,fontsize=20)
            
            
        if len(fit2g) == 6:
            l = (i) % (nrow*ncol)
            gaussfity = gaussian(selwl,fit1g[0],fit1g[1],fit1g[2])
            gauss2fity = double_gaussian(selwl,fit2g[0],fit2g[1],fit2g[2],\
                                         fit2g[3],fit2g[4],fit2g[5])
                
            comp1fity = gaussian(selwl,fit2g[0],fit2g[1],fit2g[2])
            comp2fity = gaussian(selwl,fit2g[3],fit2g[4],fit2g[5])
            #gauss2negfity = double_gaussian(selwl,fit2gneg[0],fit2gneg[1],\
                                            # fit2gneg[2],fit2gneg[3],fit2gneg[4],
                                            # fit2gneg[5])
            

            ax.flatten()[l].plot(selwlshift,sel/1e6,label='data')
            #ax.flatten()[i].plot(selwlshift,gaussfity,label='G1')
            ax.flatten()[l].plot(selwlshift,gauss2fity/1e6,label='G2',
                                 color=muted[2])
            ax.flatten()[l].plot(selwlshift,comp1fity/1e6,label='G2,C1',
                                 color=muted[4])
            ax.flatten()[l].plot(selwlshift,comp2fity/1e6,label='G2,C2',
                                 color=muted[6])
            #ax.flatten()[i].plot(selwl,gauss2negfity,label='Gauss2neg')
            #ax.flatten()[i].legend()
            ax.flatten()[l].axis(ymin=0,ymax=(maxprofile+lim)/1e6)
            ax.flatten()[l].axvline(0,linestyle='dotted')
            secaxx = ax.flatten()[l].secondary_xaxis('top', 
                                                     functions=(veltrans,wltrans))
            ax.flatten()[l].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))  
            ax.flatten()[l].set_title(times[start+i][-12:-4]+' UT')
            ax.flatten()[l].set_ylim([0,yhigh/1e6])
        else: 
            l = (i) % (nrow*ncol)
            ax.flatten()[l].plot(selwlshift,sel/1e6,label='data')
            ax.flatten()[l].axis(ymin=0,ymax=(maxprofile+lim)/1e6)
            ax.flatten()[l].axvline(0,linestyle='dotted')
            secaxx = ax.flatten()[l].secondary_xaxis('top', 
                                                     functions=(veltrans,wltrans))
            ax.flatten()[l].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))  
            ax.flatten()[l].set_title(times[start+i][-12:-4]+' UT')
            ax.flatten()[l].set_ylim([0,yhigh/1e6])
            
    secaxx = ax.flatten()[0].secondary_xaxis('top', functions=(veltrans,wltrans))
    secaxx.set_xlabel(r'Velocity $[km\; s^{-1}]$')
    ax.flatten()[0].set_xlabel(r' $\lambda$ - $\lambda_0$ [nm]')
    ax.flatten()[0].set_ylabel(r'Intensity (- $I_{min}$) $[10^6\; W\; cm^{-2} sr^{-1} \AA^{-1}]$')

    plt.tight_layout(pad = 4)
    plt.show()
    
    fig.savefig('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/'+\
                pid+'/fits'+date+note+'_endatstep_'+str(i)+'.png')
    
    return None


def perclevels(bkgd_subtract_flaretime,dispersion_range,caII_low,caII_high,
               store_ten_width,store_quarter_width,store_half_width):
    
    # Initial code, yet to be expanded to other wavelength ranges or proposals,
    # for determining the line width at different heights along emission line
    # profile; useful in tracking different types of broadening (or different
    # atmospheric heights) - see Graham and Cauzzi, 2015, 2020
    
    for i in range(np.shape(bkgd_subtract_flaretime)[0]):
        sel = bkgd_subtract_flaretime[i,caII_low:caII_high,1350]-\
            min(bkgd_subtract_flaretime[i,caII_low:caII_high,1350])
        selwl = dispersion_range[caII_low:caII_high]
        
        maxsel = np.max(sel)
        minsel = np.min(sel)
        
        ten_perc_level = 0.1*(maxsel-minsel)
        quarter_perc_level = 0.25*(maxsel-minsel)
        half_level = 0.5*(maxsel-minsel)
        
        lensel=len(sel)
        
        ten_lev_low, ten_ind_low = find_nearest(sel[0:round(lensel/2)],
                                                ten_perc_level)
        ten_lev_high, ten_ind_high = find_nearest(sel[round(lensel/2):],
                                                  ten_perc_level)
        quarter_lev_low, quarter_ind_low = find_nearest(sel[0:round(lensel/2)],
                                                        quarter_perc_level)
        quarter_lev_high, quarter_ind_high = find_nearest(sel[round(lensel/2):],
                                                          quarter_perc_level)
        half_lev_low, half_ind_low = find_nearest(sel[0:round(lensel/2)],
                                                  half_level)
        half_lev_high, half_ind_high = find_nearest(sel[round(lensel/2):],
                                                    half_level)
        
        ten_perc_width = dispersion_range[round(lensel/2)+ten_ind_high]-\
            dispersion_range[ten_ind_low]
        quarter_perc_width = dispersion_range[round(lensel/2)+quarter_ind_high]\
            -dispersion_range[quarter_ind_low]
        half_perc_width = dispersion_range[round(lensel/2)+half_ind_high]\
            -dispersion_range[half_ind_low]
        
        store_ten_width.append(ten_perc_width)
        store_quarter_width.append(quarter_perc_width)
        store_half_width.append(half_perc_width)
        
    return store_ten_width, store_quarter_width, store_half_width

# Co-alignment routines

def space_range(hdul1,pid='1_84'):
    
    # Define spatial range for co-alignment given in L1 headers
    
    x_cent = hdul1[1].header['CRVAL2']
    y_cent = hdul1[1].header['CRVAL3']
    
    x_delt = hdul1[1].header['CDELT2']
    y_delt = hdul1[1].header['CDELT3']
    
    if pid=='2_11':
        x_cent = hdul1[1].header['CRVAL1']
        y_cent = hdul1[1].header['CRVAL3']
        
        x_delt = hdul1[1].header['CDELT1']
        y_delt = hdul1[1].header['CDELT3']       
    
    nspace = hdul1[1].header['NAXIS1']
    
    #x_range is helioprojective latitude position along slit
    #y_range is helioprojective longitude position of raster step
    x_range = np.linspace(x_cent-x_delt*(nspace-1)/2,x_cent+x_delt*(nspace-1)/2,
                          nspace)
    y_range = np.linspace(y_cent-y_delt*(nspace-1)/2,y_cent+y_delt*(nspace-1)/2,
                          nspace)
    
    arcsec_slit = np.linspace(0,nspace*x_delt,nspace)
    return x_cent, y_cent, x_delt, y_delt, x_range, y_range, arcsec_slit, nspace

def vispranges(hdul1,spatial_range,nslitpos=4,pid='1_84'):
    
    # Define spatial and wavelength ranges for ViSP; this takes all 
    # slit positions in a single raster scan and uses that as a second spatial
    # axis for the "ViSP image" which will be used to co-align with VBI
    
    if pid == '2_11':
        slitlen = hdul1[1].header['CDELT1']*len(spatial_range) #in arcsec
        rastersize = hdul1[1].header['CDELT3']*nslitpos
    else:
        slitlen = hdul1[1].header['CDELT2']*len(spatial_range) #in arcsec
        rastersize = hdul1[1].header['CDELT3']*nslitpos
    
    raster_range = np.linspace(0,hdul1[1].header['CDELT3']*nslitpos,nslitpos+1)
    
    spatial_range2 = np.insert(spatial_range,0,spatial_range[0]-
                               (spatial_range[1]-spatial_range[0]))
    
    return spatial_range2, raster_range, slitlen, rastersize

def imgprep(path,folder1,dir_list2,startstep,endstep,pid='1_84'):
    
    # Prepare initial image in the ViSP set, for comparison to VBI. Could be 
    # any, but make sure the correct timestamp
    
    image_data_arrs0 = []
    rasterpos = []
    times = []
    
    image_data_arrs = []
    image_date = []

    for i in range(startstep,endstep,1):
        i_file = fits.open(path+folder1+'/'+dir_list2[i])
        
        times.append(i_file[1].header['DATE-BEG'])
        
        lammin = i_file[1].header['WAVEMIN']
        lamcen = i_file[1].header['CRVAL1']
        
        if pid=='2_11':
            lamcen = i_file[1].header['CRVAL2']
        
        i_data = i_file[1].data[0]
        
        image_data_arrs0.append(i_data)
    
        rasterpos.append(i_file[1].header['CRPIX3'])
        
    return image_data_arrs0

def line_avg(image_data_arr_arr,lowind,highind,nslit,nwave):
    caiiavgs = np.zeros((nslit,nwave))
    
    # define the boundaries (in dispersion direction) for the line and get an 
    # "average" intensity; could also use line flux for all positions along slit?
    # This simply gives an idea, when averaged lines from all slit positions are
    # plotted side-by-side, for the kernel locations.  Should use a low and high
    # index which approximately straddles the main line; Ca II H in the case of
    # the code originally developed for pid_1_84
    
    for i in range(nslit):
        caiiavgs[i,:] = np.mean(image_data_arr_arr[i][lowind:highind,:],0)


            
    return caiiavgs

def pltraster(caiiavgs,raster_range,spatial_range2,pid='pid_1_84'):
    
    # plot the intensity images for the ViSP scan
    
    X,Y = np.meshgrid(raster_range,spatial_range2)
    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,np.transpose(caiiavgs),cmap='gray')
    ax.set_aspect('equal')
    
    plt.show()
    fig.savefig('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/'+pid+
                '/initslit.png')
    
    return None


def vbi_process(path_vbi,folder1_vbi,timestep):
    
    # Process VBI data similarly to ViSP above; only intensity files
    
    dir_list_vbi = os.listdir(path_vbi+folder1_vbi)

    dir_list2_vbi = []
    
    for i in range(len(dir_list_vbi)):
        filename = dir_list_vbi[i]
        if filename[-5:] == '.fits' and '_I_' in filename:
            dir_list2_vbi.append(filename)
    
    dir_list2_vbi.sort()
    
    dir_list2_vbi
    
    hdul1_vbi = fits.open(path_vbi+folder1_vbi+'/'+dir_list2_vbi[timestep])
    dat0_vbi = hdul1_vbi[1].data[0,:,:]
    
    xcent = hdul1_vbi[1].header['CRVAL1']
    xnum = hdul1_vbi[1].header['NAXIS1']
    xdelt = hdul1_vbi[1].header['CDELT1']
    
    ycent = hdul1_vbi[1].header['CRVAL2']
    ynum = hdul1_vbi[1].header['NAXIS2']
    ydelt = hdul1_vbi[1].header['CDELT2']
    
    xarr = np.linspace(((xcent-xdelt/2)-((xnum-1)/2)*xdelt),
                       ((xcent-xdelt/2)+((xnum-1)/2)*xdelt),xnum)
    yarr = np.linspace(((ycent-ydelt/2)-((ynum-1)/2)*ydelt),
                       ((ycent-ydelt/2)+((ynum-1)/2)*ydelt),ynum)
    
    vbi_X,vbi_Y = np.meshgrid(xarr,np.flip(yarr))
    
    dat0_vbi = hdul1_vbi[1].data[0,:,:]
    
    return vbi_X,vbi_Y,hdul1_vbi, dat0_vbi

def plt_precoalign(vbi_X, vbi_Y, hdul1_vbi, visp_X, visp_Y, vispimg,matplotlib, 
                   dat0_vbi,pid='pid_1_84'):
    
    # Plot VBI and ViSP, prior to co-alignment, together.  Use ginput to ID 
    # points for similar structures.  The result of the process which this
    # starts will be ViSP axes transformed into the VBI coordinate system; 
    # which is still not correct, but we can use to simultaneously get ViSP
    # and VBI in the SDO image frfname
    
    # VBI is much higher-resolution, obviously, so easier/more effective in 
    # comparison to appropriate SDO bandpass.  
    
    # Minor errors in comparison of features lead to big errors in
    # transformation matrix, particularly farther from the features; so (1) 
    # choose features wisely (localized, bright); (2) choose bandpasses wisely
    # (similar atmospheric height, contribution function); (3) select points
    # in uncaffeinated state (jitters)


    gs = matplotlib.gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[3, 1])
    fig=plt.figure(figsize=(10,5))
    ax=fig.add_subplot(gs[0,0])
    ax1=fig.add_subplot(gs[0,1])
    ax.pcolormesh(vbi_X,vbi_Y,hdul1_vbi[1].data[0],cmap='grey')
    #ax[0].set_aspect('equal')
    ax.grid()
    ax1.pcolormesh(visp_X,visp_Y,np.transpose(vispimg),cmap='grey')
    ax1.grid()
    
    
    #plt.tight_layout()
    
    plt.show()
    
    # ginput - press on vbi first, then visp, and so on, 
    # until three pairs collected
    # be VERY careful in selection of points! Require basis vectors
    # for each coordinate system; points should not be
    # colinear
    
    #matplotlib.use('macosx')
    aa = plt.ginput(6,timeout = 120)
    
    fig.savefig('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/'+pid+
                'pre_coalign.png')
    
    return aa

def vbi_visp_transformation(aa, visp_X,visp_Y,matplotlib,nslit=91,vbi_X=[],vbi_Y=[],
                            npos=2556,vbiband='H-alpha',vispimg=[],dat0_vbi=[],
                            vispband='CaII H 396.8 nm',pid='pid_1_84',plot=0,d1=0):
    
    # Simple transformation matrix between ViSP and VBI using output of ginput
    # process in function above
    
    #d1 should be 1 if the input vectors are 1d (like in the masking scheme)
    
    # visp points
    A1 = np.array(aa[1])
    B1 = np.array(aa[3])
    C1 = np.array(aa[5])
    
    # vbi points
    A2 = np.array(aa[0])
    B2 = np.array(aa[2])
    C2 = np.array(aa[4])
    
    ViSP1 = B1 - A1
    ViSP2 = C1 - A1
    VBI1 = B2 - A2
    VBI2 = C2 - A2
    
    # insert test for linear combination
    
    # create basis matrices
    
    VBI_base = np.column_stack((VBI1,VBI2))
    
    ViSP_base = np.column_stack((ViSP1, ViSP2))
    
    # change of basis matrix
    
    COB = np.matmul(VBI_base, np.linalg.inv(ViSP_base))
    
    ViSP_points = [visp_X,visp_Y]
    
    new_ViSP = np.zeros(np.shape(ViSP_points))
       
    if d1 == 1:
        for i in range(len(visp_X)):
            point_x = visp_X[i]
            point_y = visp_Y[i]
            point = [point_x,point_y]
            ViSPvec = point - A1
            A2_1 = np.matmul(COB,ViSPvec)+A2
            
            new_ViSP[:,i] = A2_1
            
        visp_X_new = new_ViSP[0,:]
        visp_Y_new = new_ViSP[1,:]
    
    else:
        for i in range(npos):
            for j in range(nslit):
                point_x = visp_X[i,j]
                point_y = visp_Y[i,j]
                
                point = [point_x,point_y]
                ViSPvec = point - A1
                
                A2_1 = np.matmul(COB,ViSPvec)+A2
                
                new_ViSP[:,i,j] = A2_1
            
        visp_X_new = new_ViSP[0,:,:]
        visp_Y_new = new_ViSP[1,:,:]
    
    # plot new axes
    
    if plot == 1:
        
        fig,ax=plt.subplots(1,2,figsize=(10,5),sharey=True)
        ax[0].pcolormesh(vbi_X,vbi_Y,dat0_vbi,cmap='grey')
        ax[0].set_aspect('equal')
        #ax[0].invert_xaxis()
        ax[0].grid()
        ax[1].pcolormesh(visp_X_new,visp_Y_new,np.transpose(vispimg),cmap='hot')
        ax[1].set_aspect('equal')
        ax[1].grid()
        #ax[1].invert_xaxis()
        # custom_xlim = (-455,-410)
        # custom_ylim = (285,332)
        
        # # Setting the values for all axes.
        # plt.setp(ax, xlim=custom_xlim, ylim=custom_ylim)
        
        plt.tight_layout()
        
        plt.show()
        
        fig.savefig('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/'+pid+
                    'postcalib.png')
        
        # plot overlay
        
        verts=np.array([[visp_X_new[-1,-1],visp_Y_new[-1,-1]],
                        [visp_X_new[-1,0],visp_Y_new[-1,0]],
                        [visp_X_new[0,0],visp_Y_new[0,0]],
                        [visp_X_new[0,-1],visp_Y_new[0,-1]]])
        
        fig,ax = plt.subplots(1,1,figsize = (5,5))
        ax.pcolormesh(vbi_X,vbi_Y,dat0_vbi,cmap='gray')
        ax.pcolormesh(visp_X_new,visp_Y_new,np.transpose(vispimg),cmap='hot',
                      alpha = 0.3)
        boxy = matplotlib.patches.Polygon(verts,fill=False,edgecolor='black',lw=3)
        ax.add_patch(boxy)
        ax.set_title('VBI '+vbiband+' and ViSP '+vispband,fontsize=15)
        ax.set_xlabel('Heliocentric Longitude',fontsize=10)
        ax.set_ylabel('Heliocentric Latitude',fontsize=10)
        
        plt.show()
        
        fig.savefig('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/'+pid+
                    'postcalib_overlay.png')
    
    
    return visp_X_new, visp_Y_new

def query_sdo(start_time, email, cutout, matplotlib, 
              lowerx,upperx,lowery,uppery,
              wavelength = 304,
              timesamp=2,passband = '304'):
    
    # Use sunpy Fido to query different SDO bandpasses; on 2 Dec 2023 only 304 
    # and HMI continuum functionality, but can be expanded with a little work
    
    if passband == '304':
        query = Fido.search(
            a.Time(start_time - 0.01*u.h, start_time + .1*u.h),
            a.Wavelength(wavelength*u.angstrom),
            a.Sample(timesamp*u.h),
            a.jsoc.Series.aia_lev1_euv_12s,
            a.jsoc.Notify(email),
            a.jsoc.Segment.image,
            cutout,
        )

    if passband == 'cont':
        query = Fido.search(
            a.Time(start_time - 0.1*u.h, start_time + .09*u.h),
            a.Sample(timesamp*u.h),
            a.jsoc.Series.hmi_ic_45s,
            a.jsoc.Notify(email),
            a.jsoc.Segment.continuum,
            cutout,
        )
        
            
    files = Fido.fetch(query)
    
    sdohdul = fits.open(files.data[0])
    
    xr = np.linspace(lowerx,upperx,np.shape(sdohdul[1].data)[1])
    yr = np.linspace(lowery,uppery,np.shape(sdohdul[1].data)[0])
    
    X3,Y3 = np.meshgrid(xr,yr)
    
    # three points (remember which features - choose same features and same 
    # order when using points_vbi later!)

    fig,ax = plt.subplots()
    
    ax.pcolormesh(X3,Y3,sdohdul[1].data)

    plt.show()
    
    #matplotlib.use('macosx')

    bb = plt.ginput(3,timeout = 60)
    
    return query, bb

def points_vbi(vbi_X,vbi_Y,dat0_vbi,matplotlib):
    
    # VBI coordinates, same features as in query_sdo ginput 
    
    fig,ax1 = plt.subplots()
    ax1.pcolormesh(vbi_X,vbi_Y,dat0_vbi,cmap='gray')
    #ax1.set_xticklabels([])
    #ax1.set_yticklabels([])
    ax1.set_aspect('equal')
    ax1.invert_xaxis()

    #matplotlib.use('macosx')

    cc = plt.ginput(3,timeout=60)
    
    return cc

def vbi_to_sdo(bb, cc, vbi_X, vbi_Y):
    
    # Similar transformation, using basis vectors, as in visp to vbi 
    # transformation above
    
    A1 = np.array(bb[0])
    B1 = np.array(bb[1])
    C1 = np.array(bb[2])
    
    A2 = np.array(cc[0])
    B2 = np.array(cc[1])
    C2 = np.array(cc[2])
    
    #basis vectors
    SDO1 = B1 - A1
    SDO2 = C1 - A1
    VBI1 = B2 - A2
    VBI2 = C2 - A2
    
    SDO_base = np.column_stack((SDO1,SDO2))
    VBI_base = np.column_stack((VBI1,VBI2))

    COB2 = np.matmul(SDO_base,np.linalg.inv(VBI_base))
    
    VBI_points = [vbi_X,vbi_Y]
    new_VBI = np.zeros(np.shape(VBI_points))

    for i in range(np.shape(vbi_X)[0]):
        for j in range(np.shape(vbi_X)[1]):
            point_x = vbi_X[i,j]
            point_y = vbi_Y[i,j]
            
            point = [point_x,point_y]
            VBIvec = point - A2
            
            A2_1 = np.matmul(COB2,VBIvec)+A1
            
            new_VBI[:,i,j] = A2_1
            
    vbi_X_new = new_VBI[0,:,:]
    vbi_Y_new = new_VBI[1,:,:]
    
    return vbi_X_new, vbi_Y_new, COB2, A2, A1

def visp_sdo_trans(visp_X_new,visp_Y_new, COB2, A2, A1, nspace = 2544, nwave=4):
    
    # Finally, ViSP into the coordinate system defined by SDO

    ViSP_points = [visp_X_new,visp_Y_new]
    print(np.shape(ViSP_points))

    new_ViSP = np.zeros(np.shape(ViSP_points))
    
    for i in range(nspace+1):
        for j in range(nwave+1):
            point_x = visp_X_new[i,j]
            point_y = visp_Y_new[i,j]
            
            point = [point_x,point_y]
            ViSPvec = point-A2
            
            A2_1 = np.matmul(COB2,ViSPvec)+A1
            
            new_ViSP[:,i,j] = A2_1

    visp_X_new2 = new_ViSP[0,:,:]
    visp_Y_new2 = new_ViSP[1,:,:]
    
    return visp_X_new2, visp_Y_new2

def plt_final_coalign(vbi_X_new, vbi_Y_new, dat0_vbi2, 
                      visp_X_new2, visp_Y_new2, vispdat,
                      dat0_vbi, VBIpass1 = 'TiO',VBIpass2 = 
                      'H-alpha',ViSPpass = 'Ca II H',
                      obstimestr = '19 August 2022 20:42 UT',pid='pid_1_84'):
    
    # Plot final co-alignment
    
    fig,ax = plt.subplots(1,2,figsize = (10,5))
    ax[0].pcolormesh(vbi_X_new,vbi_Y_new,dat0_vbi2,cmap='gray')
    ax[0].pcolormesh(visp_X_new2,visp_Y_new2,np.transpose(vispdat),cmap='hot',alpha = 0.3)
    ax[0].patch.set_edgecolor('black')  
    ax[1].pcolormesh(vbi_X_new,vbi_Y_new,dat0_vbi,cmap='gray')
    ax[1].pcolormesh(visp_X_new2,visp_Y_new2,np.transpose(vispdat),cmap='hot',alpha = 0.3)
    ax[1].patch.set_edgecolor('black')  
    ax[0].set_title('VBI '+VBIpass1+ ' and ViSP '+ViSPpass,fontsize=15)
    ax[1].set_title('VBI '+VBIpass2+ ' and ViSP '+ViSPpass,fontsize=15)
    ax[0].set_xlabel('Helioprojective Longitude [arcsec]',fontsize=12)
    ax[0].set_ylabel('Helioprojective Latitude [arcsec]',fontsize=12)
    ax[1].set_xlabel('Helioprojective Longitude [arcsec]',fontsize=12)
    ax[1].set_ylabel('Helioprojective Latitude [arcsec]',fontsize=12)
    ax[0].set_box_aspect(1)
    ax[1].set_box_aspect(1)
    ax[0].grid()
    ax[1].grid()
    fig.suptitle(obstimestr+ ' Co-aligned',fontsize=20)
    
    fig.tight_layout()

    plt.show()
    
    fig.savefig('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/'+pid+'finalcoalign.png')
    
    return None

def load_fts(dispersion_range):
    # Load disk-center, quiet sun profile from Neckel and Hamburg
    # atlast
    
    path = '/Users/coletamburri/Desktop/DKIST_Data_Tools_misc/speclab-python/cal_data/'
    filename = 'neckel.hamburg.atlas.disk_center_intensity.cgs.ecsv'
    
    ecsv_content = path+filename
    
    disk_center_ilam = Table.read(ecsv_content,format='ascii.ecsv')
    
    wl = disk_center_ilam['wl_ang']
    ilam = disk_center_ilam['ilam_cgs']
    
    sel = np.where((wl>np.min(dispersion_range)*10) & \
                   (wl<np.max(dispersion_range)*10)) # in A
        
    wlsel = np.take(wl,sel)[0]
    ilamsel = np.take(ilam,sel)[0]
        
    return wlsel, ilamsel

def comp_fts_to_qs(wlsel, ilamsel, dispersion_range, qs_obs,lowint=0, highint=-1,
                   timelow = 0, timehigh = -1):
    
    # qs_obs is the observations for quiet sun - average this over
    # the time dimension, and then the spatial dimension, to get an averaged
    # estimate of the quiet sun
    
    # timelow and timehigh are indices where we are sure that observations
    # are quiet sun; lowint and highint are the spatial dimension of this
    
    # first restrict observations to region which is definitely quiet sun
    
    qs_sel = qs_obs[timelow:timehigh,:,lowint:highint]
    
    time_averaged_qs = np.mean(qs_sel,axis=0)
    
    space_and_time_averaged_qs = np.mean(time_averaged_qs,axis=1)
    
    fig,ax = plt.subplots()
    
    lns1 = ax.plot(dispersion_range,space_and_time_averaged_qs,color='red',label='DKIST')
    ax0 = ax.twinx()
    lns2 = ax0.plot(wlsel,ilamsel,label='Atlas')
    ax.legend()
    ax.grid()
    ax.set_xlabel('Wavelength [nm]',fontsize=13)
    ax.set_ylabel('Non-Flare Estimate [cts]',fontsize=13)
    ax0.set_ylabel(r'Disk Center Intensity [$erg/s/cm^2/sr/\mathring A$]',fontsize=13)
    ax.set_title('Observations and Flare Atlas',fontsize=20)
    lns = lns1+lns2
    
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)
    
    plt.show()
    
    return space_and_time_averaged_qs

def calib_qs_shift(wlsel, ilamsel, dispersion_range, space_and_time_averaged_qs,
             absline1, absline2, indmins, indmaxs):
    
    """
    Calibration of shift in spectrum using quiet sun and comparison to solar 
    atlas.  Performed after loading reference spectrum and quiet Sun, using 
    two absorption lines (with indices in observation spectrum, and reference
    wavelengths from atlas) in order to determine plate scale.
    
    This is outdated as of 7 March 2024.  The whole calibration process is now
    included in the get_calibration_poly function, applied to 19 August 2024
    """
    

    dispersion_range = np.array(dispersion_range)
    
    # define locations of reference lines
    shiftsel_obs1 = np.where((dispersion_range>dispersion_range[indmins[0]]) & 
                             (dispersion_range<dispersion_range[indmaxs[0]]))
    shiftsel_obs2 = np.where((dispersion_range>dispersion_range[indmins[1]]) & 
                             (dispersion_range<dispersion_range[indmaxs[1]]))
    
    wavesel_obs1 = np.take(dispersion_range,shiftsel_obs1)[0]
    fluxsel_obs1 = np.take(space_and_time_averaged_qs,shiftsel_obs1)[0]
    
    wavesel_obs2 = np.take(dispersion_range,shiftsel_obs2)[0]
    fluxsel_obs2 = np.take(space_and_time_averaged_qs,shiftsel_obs2)[0]

    pixsel_obs1 = np.arange(0, wavesel_obs1.shape[0])
    pixsel_obs2 = np.arange(0, wavesel_obs2.shape[0])
    
    xfit1, yfit1, w0_1 = fitline(pixsel_obs1,fluxsel_obs1)
    xfit2, yfit2, w0_2 = fitline(pixsel_obs2,fluxsel_obs2)
    
    # shift in lines
    dw1 = w0_1-absline1
    dw2 = w0_2-absline2
    
    dt = absline2 - absline1

    # correction in wavelengths of reference lines
    corr1 = w0_1 + indmins[0]
    corr2 = w0_2 + indmins[1]
    
    obsdiff = corr2 - corr1
    
    #change per nm
    rat = dt/obsdiff
    
    # definition of new dispersion range
    new_dispersion_range = (np.arange(0,len(dispersion_range))-corr1)*rat + \
        absline1
    
    return new_dispersion_range, rat

def find_specsamp(line1,line2,wlsel,wave_obs,spec_obs,indmins,indmaxs):

    line1fts = find_nearest(wlsel,line1)
    line2fts = find_nearest(wlsel,line2)

    dispersion_range = np.array(wave_obs)

    # define locations of reference lines
    shiftsel_obs1 = np.where((dispersion_range>dispersion_range[indmins[0]]) & 
                             (dispersion_range<dispersion_range[indmaxs[0]]))
    shiftsel_obs2 = np.where((dispersion_range>dispersion_range[indmins[1]]) & 
                             (dispersion_range<dispersion_range[indmaxs[1]]))

    wavesel_obs1 = np.take(dispersion_range,shiftsel_obs1)[0]
    fluxsel_obs1 = np.take(spec_obs,shiftsel_obs1)[0]

    wavesel_obs2 = np.take(dispersion_range,shiftsel_obs2)[0]
    fluxsel_obs2 = np.take(spec_obs,shiftsel_obs2)[0]

    pixsel_obs1 = np.arange(0, wavesel_obs1.shape[0])
    pixsel_obs2 = np.arange(0, wavesel_obs2.shape[0])

    xfit1, yfit1, w0_1 = fitline(pixsel_obs1,fluxsel_obs1)
    xfit2, yfit2, w0_2 = fitline(pixsel_obs2,fluxsel_obs2)

    # shift in lines
    dw1 = w0_1-line1
    dw2 = w0_2-line2

    dt = line2 - line1 #dt is difference in wavelength between two lines
    
    # correction in wavelengths of reference lines
    corr1 = w0_1 + indmins[0]
    corr2 = w0_2 + indmins[1]
    
    obsdiff = corr2 - corr1
    
    #change per nm
    rat = dt/obsdiff
    
    # definition of new dispersion range
    new_dispersion_range = (np.arange(0,len(dispersion_range))-corr1)*rat + \
        line1
        
    line1obs = find_nearest(new_dispersion_range,line1)
    line2obs = find_nearest(new_dispersion_range,line2)
    
    specsamp_fts = (line2fts[0]-line1fts[0])/(line2fts[1]-line1fts[1])
    specsamp_obs = (line2obs[0]-line1obs[0])/(line2obs[1]-line1obs[1])
    
    return specsamp_fts,specsamp_obs

def get_calibration_singleval(wave_obs, spec_obs, wave_atlas, spec_atlas, 
                              limbdark_fact = 1.0, wave_idx=None, extra_weight=20., bounds=None,
                              noqs_flag = 0,noqs_ind = 20):
    """
    Get calibration offsets from fitting `spec_obs` to `spec_atlas`, assuming
    wavelength grids `wave_obs` and `wave_atlas`
    Arguments:
        wave_obs: 1D array with observed wavelengths. Must be of same size as
            `spec_obs`.
        spec_obs: 1D array with observed intensities.
        wave_atlas: 1D array with wavelengths corresponding to `spec_atlas`.
        spec_atlas: 1D array with atlas intensity profile (e.g. from
            `ISpy.spec.atlas`)
    Keyword arguments:
        limbdark_fact: default 1.0 (i.e. disk center), mulitplicative factor to
            divide observations by to account for limb darkening
        wave_idx: wavelength indices that will get `extra_weight` during while
            fitting the intensity profile (default None -> all wavelengths get
            equal weight)
        extra_weight: amount of extra weight to give selected wavelength
            positions as specified by `wave_idx` (default 20)
        bounds: list of tuples [(ifact_low, ifact_upp), (woff_low, woff_upp)]
            suggesting lower and upper bounds for fitting the intensity factor
            and wavelength offset (defaults to 1/50th and 50 times the fraction
            of `spec_atlas` to `spec_obs` for ifact, and ±0.3 for woff)
    Returns:
        calibration: 2-element array [ifact, woff] with multiplication factor
        and wavelength offset to be applied to `spec_obs` and `wave_obs`
        respectively
    Author: Carlos Diaz Baso, Gregal Vissers (ISP/SU 2020), minor modifications
        from Cole Tamburri and Rahul Yadav (CU Boulder/NSO 2023)
    
    """
    wave_obs = np.array(wave_obs)
    if wave_idx is None:
        wave_idx = np.arange(wave_obs.size)
    else:
        wave_idx = np.atleast_1d(wave_idx)

    # Correct for limb-darkening if profile to calibrate on is not from
    # disc centre (and presumably at same mu as observations)
    spec_obs = spec_obs/limbdark_fact
    
    weights = np.ones_like(wave_obs)
    if wave_idx.size is not wave_obs.size:
        weights[wave_idx] = extra_weight

    def func_to_optimise(x):
        x0 = x[0]
        x1 = x[1]
        ospec = spec_obs * x0
        atlas = np.interp(wave_obs, wave_atlas-x1, spec_atlas)
        chi2 = np.sum( (atlas-ospec)**2 * weights)
        return chi2

    if bounds is None:
        bounds = [(spec_atlas[0]/spec_obs[0]*0.05, spec_atlas[0]/spec_obs[0]*50.), (-0.3, 0.3)]
    optim = differential_evolution(func_to_optimise, bounds)
    calibration = optim.x
    
    # ONLY if observations have no quiet sun - choose index that most likely 
    # corresponds to quiet sun, e.g. for pid_1_38
    if noqs_flag == 1:
        calibration[0] = spec_atlas[noqs_ind]/spec_obs[noqs_ind]
    
    calibrated_qs = spec_obs * calibration[0]
    
    new_dispersion_range2 = wave_obs + (calibration[1])
    



    return calibration, calibrated_qs, new_dispersion_range2

def gaussian_psf(x, fwhm):
	#x = wavelength [A]
	# fwhm in [A]
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))  # Calculate sigma from FWHM
    tr = np.exp(-(x)**2 / (2 * (sigma**2)))
    tr /= tr.sum()
    return tr
#write adjustment for point spread function from atlas
def psf_adjust(wlsel,ilamsel,fwhm,new_dispersion_range,calibrated_qs,limbdarkqs,ntw,gaussian_psf):
   
    func=interp1d(wlsel,ilamsel,kind='linear',fill_value='extrapolate')
    yatlas = func(new_dispersion_range)
    #dw=new_dispersion_range[1]-new_dispersion_range[0]
    #dw=dw
    #tw=(np.arange(ntw)-ntw//2)*dw
    
    stddev_pixels = fwhm/(new_dispersion_range[1]-new_dispersion_range[0])
    gaussian_kernel = Gaussian1DKernel(stddev=stddev_pixels)
   	
    #psf = gaussian_psf(tw, fwhm) # guassian transmission profiles
    #yconv = convolve(yatlas, psf, mode='same', method='fft')
    yconv = convolve(yatlas,gaussian_kernel)
    
    return yconv
    
# write calibration function for polynomial-fitting intensity calibration
def get_calibration_poly(wave_obs, spec_obs, wave_atlas, spec_atlas,find_nearest,
                         absline1, absline2, indmins, indmaxs,cont_vals = [396.628,396.71,396.77,396.9,396.9531,396.9675,397.0346,397.0829],limbdark_fact = 1.0, 
                         wave_idx=None, extra_weight=20., bounds=None,
                         noqs_flag = 0,noqs_ind = 20,ratioshift=0.0000016,order=2):
    
    dispersion_range = np.array(wave_obs)
    
    # define locations of reference lines
    shiftsel_obs1 = np.where((dispersion_range>dispersion_range[indmins[0]]) & 
                             (dispersion_range<dispersion_range[indmaxs[0]]))
    shiftsel_obs2 = np.where((dispersion_range>dispersion_range[indmins[1]]) & 
                             (dispersion_range<dispersion_range[indmaxs[1]]))
    
    wavesel_obs1 = np.take(dispersion_range,shiftsel_obs1)[0]
    fluxsel_obs1 = np.take(spec_obs,shiftsel_obs1)[0]
    
    wavesel_obs2 = np.take(dispersion_range,shiftsel_obs2)[0]
    fluxsel_obs2 = np.take(spec_obs,shiftsel_obs2)[0]

    pixsel_obs1 = np.arange(0, wavesel_obs1.shape[0])
    pixsel_obs2 = np.arange(0, wavesel_obs2.shape[0])
    
    xfit1, yfit1, w0_1 = fitline(pixsel_obs1,fluxsel_obs1)
    xfit2, yfit2, w0_2 = fitline(pixsel_obs2,fluxsel_obs2)
    
    # shift in lines
    dw1 = w0_1-absline1
    dw2 = w0_2-absline2
    
    dt = absline2 - absline1

    # correction in wavelengths of reference lines
    corr1 = w0_1 + indmins[0]
    corr2 = w0_2 + indmins[1]
    
    obsdiff = corr2 - corr1
    
    #change per nm
    rat = (dt/obsdiff)+ratioshift
    
    # definition of new dispersion range
    new_dispersion_range = (np.arange(0,len(dispersion_range))-corr1)*rat + \
        absline1
    
    wave_obs = np.array(new_dispersion_range)
    if wave_idx is None:
        wave_idx = np.arange(wave_obs.size)
    else:
        wave_idx = np.atleast_1d(wave_idx)

    # Correct for limb-darkening if profile to calibrate on is not from
    # disc centre (and presumably at same mu as observations)
    spec_obs2 = spec_obs/limbdark_fact
    
    #this really does not work for either the intensity scale or plate scale
    #Using just the new_dispersion_range above, and the polynomial intensity calibration
    # below

    weights = np.ones_like(wave_obs)
    if wave_idx.size is not wave_obs.size:
        weights[wave_idx] = extra_weight

    def func_to_optimise(x):
        x0 = x[0]
        x1 = x[1]
        ospec = spec_obs * x0
        atlas = np.interp(wave_obs, wave_atlas-x1, spec_atlas)
        chi2 = np.sum( (atlas-ospec)**2 * weights)
        return chi2

    if bounds is None:
        bounds = [(spec_atlas[0]/spec_obs[0]*0.05, spec_atlas[0]/spec_obs[0]*50.), (-0.3, 0.3)]
    optim = differential_evolution(func_to_optimise, bounds)
    calibration = optim.x
    
    # ONLY if observations have no quiet sun - choose index that most likely 
    # corresponds to quiet sun, e.g. for pid_1_38
    if noqs_flag == 1:
        calibration[0] = spec_atlas[noqs_ind]/spec_obs[noqs_ind]
        calibrated_qs = spec_obs * calibration[0]
    
    new_dispersion_range2 = wave_obs + (calibration[1])
    
    sample_flaretime = spec_obs
    dispersion_range_fin = new_dispersion_range2
    
    # contwind0_1 = np.mean(sample_flaretime[low0:high0])
    # contwind0_1_wave = np.mean(dispersion_range_fin[low0:high0])
    # contwind1 = np.mean(sample_flaretime[low1:high1])
    # contwind1_wave = np.mean(dispersion_range_fin[low1:high1])
    # contwind2 = np.mean(sample_flaretime[low2:high2])
    # contwind2_wave = np.mean(dispersion_range_fin[low2:high2])
    # contwind3 = np.mean(sample_flaretime[low3:high3])
    # contwind3_wave = np.mean(dispersion_range_fin[low3:high3])
    # contwind4 = np.mean(sample_flaretime[low4:high4])
    # contwind4_wave = np.mean(dispersion_range_fin[low4:high4])
    # contwind5 = np.mean(sample_flaretime[low5:high5])
    # contwind5_wave = np.mean(dispersion_range_fin[low5:high5])
    # contwind6 = np.mean(sample_flaretime[low6:high6])
    # contwind6_wave = np.mean(dispersion_range_fin[low6:high6])


    # cont_int_array = [contwind0_1,contwind1,contwind2,contwind3,
    #                   contwind4,contwind5,contwind6]
    # cont_int_wave_array = [contwind0_1_wave,contwind1_wave,
    #                        contwind2_wave,contwind3_wave,
    #                        contwind4_wave,contwind5_wave,
    #                        contwind6_wave]
    if order == 0: # if we just want to scale based on one continuum point
        obs_cont_loc = []
        fts_cont_loc = []
    
        for i in cont_vals:
            obs_cont_loc.append(find_nearest(new_dispersion_range,i)[1])
            fts_cont_loc.append(find_nearest(wave_atlas,i)[1])  
            
        flux_obs = spec_obs2

        cont_flux_obs = np.take(flux_obs,obs_cont_loc)
        cont_flux_fts = np.take(spec_atlas,fts_cont_loc)
        
        cont_mult_facts = cont_flux_fts/cont_flux_obs
        print('here')
        fit_vals = cont_mult_facts
        print(fit_vals)
    else:
        obs_cont_loc = []
        fts_cont_loc = []
    
        for i in cont_vals:
            obs_cont_loc.append(find_nearest(new_dispersion_range,i)[1])
            fts_cont_loc.append(find_nearest(wave_atlas,i)[1])
    
        flux_obs = spec_obs2
    
        cont_flux_obs = np.take(flux_obs,obs_cont_loc)
        cont_flux_fts = np.take(spec_atlas,fts_cont_loc)
        
        cont_mult_facts = cont_flux_fts/cont_flux_obs
    
        
        mult_fit = np.polyfit(new_dispersion_range[obs_cont_loc],cont_mult_facts,order)
        fit_cont_mult = np.poly1d(mult_fit)
        fit_vals = fit_cont_mult(new_dispersion_range)
        

    
    # if noqs_flag==1:
    #     return cont_mult_facts, fit_vals, new_dispersion_range,calibrated_qs
    # else:
        
    return cont_mult_facts, fit_vals, new_dispersion_range, dispersion_range_fin, rat


def plot_calibration(new_dispersion_range, visp_qs_obs, wlsel, ilamsel,
                     pid='pid_1_84',wl=854.207):
    
    """
    Plotting routine for calibration - quiet sun compared to reference atlas 
    spectrum, with calibrated units
    """
    
    fig,ax = plt.subplots()
    
    lns1 = ax.plot(new_dispersion_range,visp_qs_obs,color='red',label='DKIST')
    ax0 = ax.twinx()
    lns2 = ax0.plot(wlsel,ilamsel,label='Atlas')
    ax.legend()
    ax.grid()
    ax.set_xlabel('Wavelength [nm]',fontsize=13)
    ax.set_ylabel('Calibrated Intensity [$erg/s/cm^2/sr/\mathring A$]',fontsize=13)
    ax0.set_ylabel(r'Disk Center Intensity [$erg/s/cm^2/sr/\mathring A$]',fontsize=13)
    ax.set_title('Calibrated Quiet Sun Intensity',fontsize=20)
    lns = lns1+lns2
    ax.axvline(wl)
    
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)
    
    plt.show()
    ax0.set_ylim([min(ilamsel)-.15e6,max(ilamsel)+.15e6])
    ax.set_ylim([min(ilamsel)-.15e6,max(ilamsel)+.15e6])
    
    fig.savefig('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/'+pid+\
                '/pltprofile.png')
    
    return None

def maxintind(new_dispersion_range,image_data_arr_arr,linelow,linehigh,spacelow,spacehigh):
    """
    Stores spatial indices of max averaged intensity in emission line in
    question, to better track flare kernel when working with multiple steps
    in the same raster scan (or just to better identify the kernel center
    in a single image)
    """
    
    indices = []
    for i in range(np.shape(image_data_arr_arr)[0]):
        lineavg = list(np.mean(image_data_arr_arr[i,linelow:linehigh,spacelow:spacehigh],axis=0))
        indices.append(spacelow+lineavg.index(max(lineavg)))
    
    return indices
        
def comp_fit_results_gauss2(fits_2g,times):
    """
    

    Parameters
    ----------
    fits_2g : TYPE
        DESCRIPTION.
    times : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    fig,ax = plt.subplots()
    
def conv_to_vel(mu1,mu2,mu,lamb0 = 396.847,c=2.99e5):
    vel1 = []
    vel2 = []
    
    shift1 = []
    shift2 = []
    
    for i in range(len(mu1)):
        shift1.append(mu1[i]-lamb0)
        shift2.append(mu2[i]-lamb0)
        
    for i in range(len(mu1)):
        
        vel1.append(((((shift1[i]+lamb0)/lamb0)-1)*c)/mu)
        vel2.append(((((shift2[i]+lamb0)/lamb0)-1)*c)/mu)
        
    return vel1, vel2

def plt_allrasterprofile(maxindices,scaled_flare_time,dispersion_range,end=5):
    
    fig,ax = plt.subplots(2,2,dpi=200)

    colors = pl.cm.turbo(np.linspace(0,1,44))
    cmap = pl.cm.get_cmap('turbo', 44)
    
    for i in range(0,44,4):
        ax[0,0].plot(dispersion_range[end:],scaled_flare_time[i,end:,maxindices[i]]/1e6,color=colors[i],label='Flare-Time')
        ax[0,1].plot(dispersion_range[end:],scaled_flare_time[i+1,end:,maxindices[i+1]]/1e6,color=colors[i+1],label='Flare-Time')
        ax[1,0].plot(dispersion_range[end:],scaled_flare_time[i+2,end:,maxindices[i+2]]/1e6,color=colors[i+2],label='Flare-Time')
        ax[1,1].plot(dispersion_range[end:],scaled_flare_time[i+3,end:,maxindices[i+3]]/1e6,color=colors[i+3],label='Flare-Time')
    
        ax[0,0].grid()
        ax[0,1].grid()
        ax[1,0].grid()
        ax[1,1].grid()
        
        ax[0,0].set_ylim([0,7])
        ax[0,1].set_ylim([0,7])
        ax[1,0].set_ylim([0,7])
        ax[1,1].set_ylim([0,7])
        
        ax[0,0].set_xlim([396.7,397.1])
        ax[0,1].set_xlim([396.7,397.1])
        ax[1,0].set_xlim([396.7,397.1])
        ax[1,1].set_xlim([396.7,397.1])
        
        ax[0,0].set_title('Position 1',fontsize=10)
        ax[0,1].set_title('Position 2',fontsize=10)
        ax[1,0].set_title('Position 3',fontsize=10)
        ax[1,1].set_title('Position 4',fontsize=10)
    
        #plt.suptitle(r'DKIST/ViSP Ca II H and H$\epsilon$, 2022 August 19')
    
        ax[0,0].set_xlabel(r'Wavelength [nm]',fontsize=7)
        ax[0,1].set_xlabel(r'Wavelength [nm]',fontsize=7)
        ax[1,0].set_xlabel(r'Wavelength [nm]',fontsize=7)
        ax[1,1].set_xlabel(r'Wavelength [nm]',fontsize=7)
    
        ax[0,0].set_ylabel(r'Intensity [$10^6\;W\; cm^{-2}\; sr^{-1}\; \mathring A^{-1}$]',fontsize=7)
        ax[0,1].set_ylabel(r'Intensity [$10^6\;W\; cm^{-2}\; sr^{-1}\; \mathring A^{-1}$]',fontsize=7)
        ax[1,0].set_ylabel(r'Intensity [$10^6\;W\; cm^{-2}\; sr^{-1}\; \mathring A^{-1}$]',fontsize=7)
        ax[1,1].set_ylabel(r'Intensity [$10^6\;W\; cm^{-2}\; sr^{-1}\; \mathring A^{-1}$]',fontsize=7)
        
        ax[0,0].tick_params(axis='both', which='major', labelsize=7)
        ax[0,1].tick_params(axis='both', which='major', labelsize=7)
        ax[1,0].tick_params(axis='both', which='major', labelsize=7)
        ax[1,1].tick_params(axis='both', which='major', labelsize=7)

    plt.show()
    
    fig.tight_layout()
    
    fig.savefig('/Users/coletamburri/Desktop/allrasters_separate.png',dpi=400)
    
    fig,ax = plt.subplots(figsize=(3,3),dpi=200)
    
    cmap = pl.cm.get_cmap('turbo', 44)
        
    for i in range(0,44,4):
        ax.plot(dispersion_range[end:],scaled_flare_time[i,end:,maxindices[i]]/1e6,linewidth=1,color=colors[i],label='Flare-Time')
        ax.plot(dispersion_range[end:],scaled_flare_time[i+1,end:,maxindices[i+1]]/1e6,linewidth=1,color=colors[i+1],label='Flare-Time')
        ax.plot(dispersion_range[end:],scaled_flare_time[i+2,end:,maxindices[i+2]]/1e6,linewidth=1,color=colors[i+2],label='Flare-Time')
        ax.plot(dispersion_range[end:],scaled_flare_time[i+3,end:,maxindices[i+3]]/1e6,linewidth=1,color=colors[i+3],label='Flare-Time')
    
        ax.grid()
    
        
        ax.set_ylim([0,7])
    
        
        ax.set_xlim([396.7,397.1])

        #plt.suptitle(r'DKIST/ViSP Ca II H and H$\epsilon$, 2022 August 19')
    
        ax.set_xlabel(r'Wavelength [nm]',fontsize=7)
        ax.set_ylabel(r'Intensity [$10^6\;W\; cm^{-2}\; sr^{-1}\; \mathring A^{-1}$]',fontsize=7)
        ax.tick_params(axis='both', which='major', labelsize=7)
    
    fig.show()

    fig.savefig('/Users/coletamburri/Desktop/allrasters_together.png',dpi=400)

    return None

        
def prep_arms(base,folder_arm1,folder_arm2,file_arm1,file_arm2,startind=2548,\
              nslit = 91):
    
    # function to get coordinates (from L1 header) for two DKIST arms, could 
    # be any.  startind and nslit are the start index for a good-seeing scan
    # to do the alignment and the number of slits per scan.
    
    dir_list_arm1 = pathdef(base,folder_arm1) 
    dir_list_arm2 = pathdef(base,folder_arm2) 

    goodscan_arm1 = dir_list_arm1[startind:startind+nslit]
    goodscan_arm2 = dir_list_arm2[startind:startind+nslit]

    #function to get co-ordinates along slit
    firststep = fits.open(base+folder_arm1+goodscan_arm1[0])
    
    center = firststep[1].header['CRVAL1']
    dy = firststep[1].header['CDELT1']
    center_unit = firststep[1].header['CUNIT1']
    n1 = firststep[1].header['NAXIS1']
    
    yarr_arm1 = np.arange(center+(n1*dy/2),center-(n1*dy/2),-dy)
    
    #function to get co-ordinates along slit
    firststep = fits.open(base+folder_arm2+goodscan_arm2[0])
    
    center = firststep[1].header['CRVAL1']
    dy = firststep[1].header['CDELT1']
    center_unit = firststep[1].header['CUNIT1']
    n1 = firststep[1].header['NAXIS1']
    
    yarr_arm2 = np.arange(center+(n1*dy/2),center-(n1*dy/2),-dy)
    
    # x coordinates
    step = fits.open(base+folder_arm2+goodscan_arm2[0])
    center = firststep[1].header['CRVAL3']
    dx = firststep[1].header['CDELT3']
    center_unit = firststep[1].header['CUNIT3']
    n3 = 91
    
    xarr_arm2 = np.arange(center+(n3*dx/2),center-(n3*dx/2),-dx)
    
    #hbeta
    
    step = fits.open(base+folder_arm1+goodscan_arm1[0])
    center = firststep[1].header['CRVAL3']
    dx = firststep[1].header['CDELT3']
    center_unit = firststep[1].header['CUNIT3']
    n3 = 91
    
    xarr_arm1 = np.arange(center+(n3*dx/2),center-(n3*dx/2),-dx)
    
    return xarr_arm1, yarr_arm1, xarr_arm2, yarr_arm2, dir_list_arm1, \
        dir_list_arm2
    
            
        
    

        


    
    
    
    
    
    
    
    
        
    

    




    
    

    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    







