#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 16:10:47 2024

@author: coletamburri
A script to compare RADYN+RH to DKIST
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib as mpl
import tol_colors
from scipy.interpolate import interp1d
from scipy.signal import convolve

def normalize(arr):
    norm_arr = []
    minimum = np.nanmin(arr)
    maximum = np.nanmax(arr)
    diff_arr = maximum - minimum   
    for i in arr:
        temp = (i - minimum)/(maximum-minimum)
        norm_arr.append(temp)
    return norm_arr

def normalize_range(arr,lowlim,highlim):
    norm_arr = []
    minimum = np.nanmin(arr[lowlim:highlim])
    maximum = np.nanmax(arr[lowlim:highlim])
    diff_arr = maximum - minimum   
    for i in arr:
        temp = (i - minimum)/(maximum-minimum)
        norm_arr.append(temp)
    return norm_arr    

c=2.99e5
lamb0 = 396.847
lamb1 = 397.01
mu = 0.4760111410077789
mu2 = 0.4266927415494022

#ViSP instrument
fwhm = .005 # in nm
ntw = 45
flagh20 = 1
flagh20sum = 0

#only one of the following should be 1
flagb = 0
flagvt = 0
flagt = 0

heplowh20 = 690
hephighh20 = 850


def veltrans(x):
    return ((((x+lamb0)/lamb0)-1)*c)/mu

def veltrans2(x):
    return ((((x+lamb0)/lamb0)-1)*c)/mu2

def wltrans(x):
    return ((((x/c)+1)*lamb0)-lamb0)

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

def tocgs(self, w, s):
    clight=2.99792458e10         #speed of light [cm/s]
    joule_2_erg=1e7
    aa_to_cm=1e-8
    s *=joule_2_erg/aa_to_cm # from Watt /(cm2 ster AA) to erg/(s cm2 ster cm)
    s *=(w*aa_to_cm)**2/clight   # to erg/
    return s

def tosi(self, wav, s):
    clight=2.99792458e8      #speed of light [m/s]                                  
    aa_to_m=1e-10                                                                        
    cm_to_m=1e-2                       
    s /= cm_to_m**2 * aa_to_m # from from Watt /(s cm2 ster AA) to Watt/(s m2 ster m) 
    s *= (wav*aa_to_m)**2 / clight # to Watt/(s m2 Hz ster)
    return s

#example convolution with gaussian psf
def gaussian_psf(x, fwhm):
	#x = wavelength [nm]
	# fwhm in [nm]
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))  # Calculate sigma from FWHM
    tr = np.exp(-(x)**2 / (2 * (sigma**2)))
    tr /= tr.sum()
    return tr
#write adjustment for point spread function from atlas
def psf_adjust(wlsel,ilamsel,fwhm,new_dispersion_range,ntw,gaussian_psf):
    func=interp1d(wlsel,ilamsel,kind='linear',fill_value='extrapolate')
    yatlas = func(new_dispersion_range)
    dw=new_dispersion_range[1]-new_dispersion_range[0]
    dw=dw
    tw=(np.arange(ntw)-ntw//2)*dw
    
    for i in range(1):

    	psf = gaussian_psf(tw, fwhm) # guassian transmission profiles
    	yconv = convolve(yatlas, psf, mode='same', method='fft')
        
    return yconv

times=[500]

ncol=len(times)
map = tol_colors.tol_cmap(colormap='rainbow_discrete',lut=ncol+5)
cmap_choice = map(np.linspace(0,1,ncol))

base = '/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_output_files_npz/'
modelnameqs = 'cat_15_8_5e10_20_600_0.npz'

filename_dkist = '/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_output_files_npz/dkist_obs_file.npz'
filename_dkist_sub = '/Users/coletamburri/Desktop/Misc_DKIST/subbed.npz'

model_choiceqs = np.load(base+modelnameqs)

dkist_file = np.load(filename_dkist)
dkist_file_sub = np.load(filename_dkist_sub)

caiih_inds = np.where((model_choiceqs['wl_rh']>396.7) & (model_choiceqs['wl_rh']< 397.94))

model_choiceqs_wl = model_choiceqs['wl_rh'][caiih_inds]
model_choiceqs_int = model_choiceqs['int_rh'][caiih_inds]
model_choiceqs_wlshift = model_choiceqs_wl-lamb0
model_choiceqs_wlshift_hep = model_choiceqs_wl-lamb1

maxind = dkist_file['arr_3'][0]

dkist_wl = dkist_file_sub['arr_2']
#dkist_int = dkist_file['arr_0'][0,:,maxind] #before background subtraction
dkist_int = dkist_file_sub['arr_1'][0,:,maxind] #after background subtraction

dkist_int_trail = dkist_file_sub['arr_0'][0,:,maxind+120]
dkist_int_front = dkist_file_sub['arr_0'][0,:,maxind-100]

dkist_wl_shift = dkist_wl-lamb0
dkist_wl_shift_hep = dkist_wl-lamb1

yconvqs = psf_adjust(model_choiceqs_wl,model_choiceqs_int,fwhm,dkist_wl,ntw,gaussian_psf)


modelname1 = 'cat_15_8_5e10_20_600_500.npz'

#chosen model to compare (can/will be many)
model_choice1 = np.load(base+modelname1)
caiih_indsh20 = np.where((model_choice1['wl_rh']>396.7) & (model_choice1['wl_rh']< 397.94))

model_choice1_wl = model_choice1['wl_rh'][caiih_indsh20]
model_choice1_int = model_choice1['int_rh'][caiih_indsh20]
model_choice1_wlshift = model_choice1_wl-lamb0
model_choice1_wlshift_hep = model_choice1_wl-lamb1
    
# adjust for instrument PSF

model1_copy_int = model_choice1_int
model1_copy_wl = model_choice1_wl
yconvqs = psf_adjust(model_choiceqs_wl,model_choiceqs_int,fwhm,dkist_wl,ntw,gaussian_psf)

yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                               gaussian_psf)

model_subtract = yconv1-yconvqs


    
    
    
    