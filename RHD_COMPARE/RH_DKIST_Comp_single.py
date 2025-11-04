#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  4 17:57:07 2025

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
lamb2 = 486.1
mu = 0.4760111410077789
mu2 = 0.4266927415494022

#ViSP instrument
fwhm = .004 # in nm
ntw = 50
flagh20 = 1
flagh20sum = 0

#only one of the following should be 1
flagb = 0
flagvt = 0
flagt = 0

hepwl=397.01

heplowh20 = 690
hephighh20 = 850

caII_low = 630
caII_high = 670
hep_low = 730
hep_high = 900



lowvisp=148
highvisp=148+91
#lowvisp=0
#highvisp=-1

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

hbeta = 0


base = '/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_output_files_npz/'
modelname = '11Aug_Cclass_A_final_36s_mu0.5.npz'

modelnameqs = 'cat_15_8_5e10_20_600_0.npz'

dkist_file = '/Users/coletamburri/Desktop/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/ViSP_spectra_processed_11Aug24_CaII.npz'
if hbeta==1:
    dkist_file_hb = '/Users/coletamburri/Desktop/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/ViSP_spectra_processed_11Aug24_Hbeta.npz'

model_choiceqs = np.load(base+modelnameqs)

dkist_file = np.load(dkist_file)
if hbeta==1:
    dkist_file_hbeta = np.load(dkist_file_hb)

caiih_inds = np.where((model_choiceqs['wl_rh']>396.7) & (model_choiceqs['wl_rh']< 397.94))
if hbeta==1:
    hbeta_inds= np.where((model_choiceqs['wl_rh']>485.95) & (model_choiceqs['wl_rh']< 486.25))

model_choiceqs_wl = model_choiceqs['wl_rh'][caiih_inds]
model_choiceqs_int = model_choiceqs['int_rh'][caiih_inds]
model_choiceqs_wlshift = model_choiceqs_wl-lamb0
model_choiceqs_wlshift_hep = model_choiceqs_wl-lamb1

if hbeta==1:
    model_choiceqs_wl_hbeta = model_choiceqs['wl_rh'][hbeta_inds]
    model_choiceqs_int_hbeta = model_choiceqs['int_rh'][hbeta_inds]
    model_choiceqs_wlshift_hbeta = model_choiceqs_wl_hbeta-lamb2

#maxind =
dkist_wl = dkist_file['wl']
dkist_int = dkist_file['flare']
dkist_avg = np.mean(dkist_int[lowvisp:highvisp,caII_low:caII_high,:],1)
dkist_time = dkist_file['time']

if hbeta==1:
    dkist_wl_hbeta = dkist_file_hbeta['wl']
    dkist_int_hbeta = dkist_file_hbeta['flare']
    dkist_avg_hbeta = np.mean(dkist_int_hbeta[lowvisp:highvisp,422:706,:],1)
    dkist_time_hbeta = dkist_file_hbeta['time']

dkist_wl_shift = dkist_wl-lamb0
dkist_wl_shift_hep = dkist_wl-lamb1

yconvqs = psf_adjust(model_choiceqs_wl,model_choiceqs_int,fwhm,dkist_wl,ntw,gaussian_psf)

model_choice1 = np.load(base+modelname)
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

model_subtract1 = yconv1-yconvqs

#for hbeta
if hbeta==1:
    model_choice1_hb = np.load(base+modelname)
    hb_inds = np.where((model_choice1['wl_rh']>485.95) & (model_choice1['wl_rh']< 486.25))
    
    model_choice1_wl_hb = model_choice1_hb['wl_rh'][hb_inds]
    model_choice1_int_hb = model_choice1_hb['int_rh'][hb_inds]
    model_choice1_wlshift_hb = model_choice1_wl_hb-lamb2
        
    # adjust for instrument PSF
    
    model1_copy_int_hb = model_choice1_int_hb
    model1_copy_wl_hb = model_choice1_wl_hb
    yconvqs_hb = psf_adjust(model_choiceqs_wl_hbeta,model_choiceqs_int_hbeta,fwhm,dkist_wl_hbeta,ntw,gaussian_psf)
    
    yconv1_hb=psf_adjust(model1_copy_wl_hb,model1_copy_int_hb,fwhm,dkist_wl_hbeta,ntw,
                                   gaussian_psf)
    
    model_subtract1_hb = yconv1_hb-yconvqs_hb

#get pixel
npoints = 9

fig,ax=plt.subplots()
ax.pcolormesh(np.transpose(dkist_avg),cmap='inferno')
ax.invert_xaxis()
ax.invert_yaxis()
fig.show()

cc = plt.ginput(npoints,timeout=120)

# fig,[ax,ax2]=plt.subplots(1,2)
# ax.plot(dkist_wl,normalize_range(model_subtract1/1e6,caII_low,hep_high),alpha=1,c='b',linewidth=2,label='EB1, 500s')
# for i in range(len(cc)):
#     ax.plot(dkist_wl,normalize_range(dkist_int[int(cc[i][1]),:,int(cc[i][0])]/1e6,caII_low,hep_high),alpha=1,linewidth=2,label='DKIST/ViSP')
# ax.set_xlim([396.6,397.1])
# ax2.pcolormesh(dkist_avg)
# for i in range(len(cc)):
#     ax2.scatter(int(cc[i][0]),int(cc[i][1]),marker='x')
# ax.axvline(396.85)
# ax.axvline(397.01)

# fig.show()


colors = plt.cm.turbo(np.linspace(0,1,npoints))

fig,ax=plt.subplots(3,3,dpi=150)
for i in range(len(cc)):
    ax.flatten()[i].plot(dkist_wl,model_subtract1/1e6,alpha=.7,c='k',linewidth=2,linestyle='dashed')
    ax.flatten()[i].plot(dkist_wl,dkist_int[lowvisp+int(cc[i][0]),:,int(cc[i][1])]/1e6,alpha=1,linewidth=2,c=colors[i])
    ax.flatten()[i].set_xlim([396.6,397.1])
    ax.flatten()[i].axvline(396.85,c='black',linestyle='dotted',linewidth=1)
    ax.flatten()[i].axvline(397.01,c='black',linestyle='dotted',linewidth=1)
    ax.flatten()[i].set_xticks([396.6,396.8,397])
    #ax.flatten()[i].set_xlabel('Wavelength [nm]')
    #ax.flatten()[i].set_ylabel('Intensity')

fig.show()

fig,ax2=plt.subplots(dpi=200,figsize=(2,10))
ax2.pcolormesh(np.transpose(dkist_avg),cmap='Reds')
for i in range(len(cc)):
    ax2.scatter(int(cc[i][0]),int(cc[i][1]),marker='x',c=colors[i])
ax2.invert_xaxis()
ax2.invert_yaxis()

ax2.set_ylim([1500,950])
fig.show()

#now for hbeta
if hbeta==1:
    npoints = 4
    
    fig,ax=plt.subplots()
    ax.pcolormesh(np.transpose(dkist_avg_hbeta),cmap='Reds')
    ax.invert_xaxis()
    ax.invert_yaxis()
    fig.show()
    
    cc = plt.ginput(npoints,timeout=120)
    
    # fig,[ax,ax2]=plt.subplots(1,2)
    # ax.plot(dkist_wl,normalize_range(model_subtract1/1e6,caII_low,hep_high),alpha=1,c='b',linewidth=2,label='EB1, 500s')
    # for i in range(len(cc)):
    #     ax.plot(dkist_wl,normalize_range(dkist_int[int(cc[i][1]),:,int(cc[i][0])]/1e6,caII_low,hep_high),alpha=1,linewidth=2,label='DKIST/ViSP')
    # ax.set_xlim([396.6,397.1])
    # ax2.pcolormesh(dkist_avg)
    # for i in range(len(cc)):
    #     ax2.scatter(int(cc[i][0]),int(cc[i][1]),marker='x')
    # ax.axvline(396.85)
    # ax.axvline(397.01)
    
    # fig.show()
    
    
    colors = plt.cm.jet(np.linspace(0,1,npoints))
    
    fig,ax=plt.subplots(2,2,dpi=200)
    for i in range(len(cc)):
        ax.flatten()[i].plot(dkist_wl_hbeta,model_subtract1_hb/1e6,alpha=.7,c='k',linewidth=2)
        ax.flatten()[i].plot(dkist_wl_hbeta,dkist_int_hbeta[lowvisp+int(cc[i][0]),:,int(cc[i][1])]/1e6,alpha=1,linewidth=1,c=colors[i])
        #ax.flatten()[i].set_xlim([396.6,397.1])
        ax.flatten()[i].axvline(486.1)

        #ax.flatten()[i].axvline(397.01)
    fig.show()
    
    fig,ax2=plt.subplots(dpi=200,figsize=(2,10))
    ax2.pcolormesh(np.transpose(dkist_avg_hbeta),cmap='inferno')
    for i in range(len(cc)):
        ax2.scatter(int(cc[i][0]),int(cc[i][1]),marker='x',c=colors[i])
    ax2.invert_xaxis()
    ax2.invert_yaxis()
    ax2.set_xlim([100,50])
    ax2.set_ylim([1500,950])
    fig.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
