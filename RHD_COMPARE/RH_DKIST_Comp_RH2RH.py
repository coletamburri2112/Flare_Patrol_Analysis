#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 18:59:56 2026

@author: coletamburri
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

hepwl=397.01

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

base = '/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_output_files_npz/'

modelnameqs = 'cat_15_8_5e10_20_600_0.npz'
model_choiceqs = np.load(base+modelnameqs)


model_choiceqs_wl = model_choiceqs['wl_rh']
model_choiceqs_int = model_choiceqs['int_rh']
model_choiceqs_wlshift = model_choiceqs_wl-lamb0
model_choiceqs_wlshift_hep = model_choiceqs_wl-lamb1

EB3_43s = 'fchroma31_20b_5vt_43s_H20.npz'
EB3_50s = 'rhf1d_fchroma_50s_1broadc_7vt_EB3.npz'

EB4_43s = 'fchroma30_20b_5vt_43s_H20.npz'
EB4_50s = 'rhf1d_fchroma_50s_1broadc_7vt_EB4.npz'


heplowh20 = 690
hephighh20 = 850

filename_dkist = '/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_output_files_npz/dkist_obs_file.npz'
filename_dkist_sub = '/Users/coletamburri/Desktop/Misc_DKIST/subbed.npz'

dkist_file = np.load(filename_dkist)
dkist_file_sub = np.load(filename_dkist_sub)

caiih_inds = np.where((model_choiceqs['wl_rh']>396.7) & (model_choiceqs['wl_rh']< 397.94))

dkist_wl = dkist_file_sub['arr_2']
maxind = 1350
#dkist_int = dkist_file['arr_0'][0,:,maxind] #before background subtraction
dkist_int = dkist_file_sub['arr_1'][0,:,maxind] #after background subtraction

EB3_43s_load = np.load(base+EB3_43s)
EB3_50s_load = np.load(base+EB3_50s)

EB4_43s_load = np.load(base+EB4_43s)
EB4_50s_load = np.load(base+EB4_50s)

wavelength = EB3_43s_load['wl_rh']

EB3_43s_int = EB3_43s_load['int_rh']
EB3_50s_int = EB3_50s_load['int_rh']

EB4_43s_int = EB4_43s_load['int_rh']
EB4_50s_int = EB4_50s_load['int_rh']

lowind,lowval = find_nearest(wavelength,396.97)
highind,highval = find_nearest(wavelength,397.1)

fig,[ax,ax0]=plt.subplots(1,2)
ax.plot(wavelength,normalize_range((EB3_43s_int-model_choiceqs_int),lowind,highind),c='black',label='43s')
ax.plot(wavelength,normalize_range((EB3_50s_int-model_choiceqs_int),lowind,highind),c='red',label='50s')
ax.plot(dkist_wl[heplowh20:hephighh20],normalize(dkist_int[heplowh20:hephighh20]),label='DKIST/ViSP',linewidth=3,c='blue',zorder=6,alpha=1)
ax.set_xlim([396.95,397.073])
ax.set_ylim([-0.15,1.05])
ax.legend()
ax.set_title('EB3')

ax0.plot(wavelength,normalize_range((EB4_43s_int-model_choiceqs_int),lowind,highind),c='black',label='43s')
ax0.plot(wavelength,normalize_range((EB4_50s_int-model_choiceqs_int),lowind,highind),c='red',label='50s')
ax0.plot(dkist_wl[heplowh20:hephighh20],normalize(dkist_int[heplowh20:hephighh20]),label='DKIST/ViSP',linewidth=3,c='blue',zorder=6,alpha=1)
ax0.set_xlim([396.95,397.073])
ax0.set_ylim([-0.15,1.05])
ax0.set_title('EB4')
ax0.legend()
fig.show()

fig,ax=plt.subplots()
ax.plot(wavelength,EB3_43s_int,c='black')
ax.plot(wavelength,EB3_50s_int,c='red')
ax.set_xlim([396.6,397.12])
#ax.set_ylim([-0.15,1.05])
fig.show()

fig,ax=plt.subplots()
ax.plot(wavelength,EB4_43s_int,c='black')
ax.plot(wavelength,EB4_50s_int,c='red')
ax.set_xlim([396.6,397.12])
#ax.set_ylim([-0.15,1.05])
fig.show()