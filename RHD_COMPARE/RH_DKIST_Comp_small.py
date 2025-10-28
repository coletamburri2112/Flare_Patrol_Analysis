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

times=[500]

flag = 'f_chroma'
h20mods = [79,19,25,31]
h20vals = [3e10,1e11,3e11,1e12]

ncol=len(times)
map = tol_colors.tol_cmap(colormap='rainbow_discrete',lut=ncol+5)
cmap_choice = map(np.linspace(0,1,ncol))

base = '/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_output_files_npz/'
modelnameqs = 'cat_15_8_5e10_20_600_0.npz'
#modelnameqs = 'radyn_TC_0s'

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


modelname1_orig = 'cat_15_8_5e10_20_600_500.npz'
model_choice1 = np.load(base+modelname1_orig)
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

model_subtract1_orig = yconv1-yconvqs


modelname1 = 'EBLD_19Aug2022_7kms.npz'
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

model_subtract1 = yconv1-yconvqs



modelnameTC_orig = 'rhf1d_5_TC_90s.npz'
#chosen model to compare (can/will be many)
model_choice1 = np.load(base+modelnameTC_orig)
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

model_subtract_TC_orig = yconv1-yconvqs

modelnameTC = 'TC_19Aug2022_7kms.npz'

#chosen model to compare (can/will be many)
model_choice1 = np.load(base+modelnameTC)
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

model_subtract_TC = yconv1-yconvqs

# ax.set_xlabel(r'Wavelength ($\lambda-\lambda_c$) [nm]')
# # secaxx = ax.secondary_xaxis('top', functions=(veltrans,wltrans))
# # secaxx.set_xlabel(r'Velocity $[km\; s^{-1}]$')
# # secaxx.set_xticks([-120,-80,-40,0,40,80,120])
# # secaxx.grid('on')

fig,ax=plt.subplots(dpi=200,figsize=(5,4))
# ax.plot(dkist_wl,normalize_range(model_subtract1/1e6,heplowh20,hephighh20),alpha=1,c='b',linewidth=2,label=r'EB1, 500s, $v_{turb} = 7\;km\;s^{-1}$')
# ax.plot(dkist_wl,normalize_range(model_subtract_TC/1e6,heplowh20,hephighh20),alpha=1,c='red',linewidth=2,label=r'TC1, 90s, $v_{turb} = 7\;km\;s^{-1}$')
# ax.plot(dkist_wl,normalize_range(model_subtract1_orig/1e6,heplowh20,hephighh20),alpha=.7,c='b',linewidth=2,linestyle='dashed',label=r'EB1, 500s, $v_{turb} = 2\;km\;s^{-1}$')
# ax.plot(dkist_wl,normalize_range(model_subtract_TC_orig/1e6,heplowh20,hephighh20),alpha=.7,c='red',linewidth=2,linestyle='dashed',label=r'TC1, 90s, $v_{turb} = 2\;km\;s^{-1}$')
# ax.plot(dkist_wl,normalize_range(dkist_int/1e6,heplowh20,hephighh20),alpha=1,c='black',linewidth=2,label='DKIST/ViSP')

# ax.plot(dkist_wl,normalize_range(model_subtract1/1e6,400,heplowh20),alpha=1,c='b',linewidth=2,label=r'EB1, 500s, $v_{turb} = 7\;km\;s^{-1}$')
# ax.plot(dkist_wl,normalize_range(model_subtract_TC/1e6,400,heplowh20),alpha=1,c='red',linewidth=2,label=r'TC1, 90s, $v_{turb} = 7\;km\;s^{-1}$')
# ax.plot(dkist_wl,normalize_range(model_subtract1_orig/1e6,400,heplowh20),alpha=.7,c='b',linewidth=2,linestyle='dashed',label=r'EB1, 500s, $v_{turb} = 2\;km\;s^{-1}$')
# ax.plot(dkist_wl,normalize_range(model_subtract_TC_orig/1e6,400,heplowh20),alpha=.7,c='red',linewidth=2,linestyle='dashed',label=r'TC1, 90s, $v_{turb} = 2\;km\;s^{-1}$')
# ax.plot(dkist_wl,normalize_range(dkist_int/1e6,400,heplowh20),alpha=1,c='black',linewidth=2,label='DKIST/ViSP')


ax.plot(dkist_wl,model_subtract1/1e6,alpha=1,c='b',linewidth=2,label=r'EB1, 500s, $v_{turb} = 7\;km\;s^{-1}$')
ax.plot(dkist_wl,model_subtract_TC/1e6,alpha=1,c='red',linewidth=2,label=r'TC1, 90s, $v_{turb} = 7\;km\;s^{-1}$')
ax.plot(dkist_wl,model_subtract1_orig/1e6,alpha=.7,c='b',linewidth=2,linestyle='dashed',label=r'EB1, 500s, $v_{turb} = 2\;km\;s^{-1}$')
ax.plot(dkist_wl,model_subtract_TC_orig/1e6,alpha=.7,c='red',linewidth=2,linestyle='dashed',label=r'TC1, 90s, $v_{turb} = 2\;km\;s^{-1}$')
ax.plot(dkist_wl,dkist_int/1e6,alpha=1,c='black',linewidth=2,label='DKIST/ViSP')

# ax.plot(dkist_wl,model_subtract1/1e6,alpha=1,c='b',linewidth=3)
# ax.plot(dkist_wl,model_subtract_TC/1e6,alpha=1,c='#CC6677',linewidth=3)
# ax.plot(dkist_wl,dkist_int/1e6,alpha=1,c='black',linewidth=3)

ax.legend(fontsize=8)
ax.set_xlim([396.75,397.1])
ax.grid()
ax.axvline(396.846,linestyle='dashed',c='grey',linewidth=2)
ax.axvline(397.01,linestyle='dashed',c='grey',linewidth=2)
ax.set_ylabel(r'Intensity (Normalized to $H\epsilon$)')
ax.set_xlabel('Wavelength [nm]')
ax.set_xticks([396.8,396.9,397.0,397.1])
plt.show()
    
models_tocomp=[]
# add H20
if flag =='f_chroma':
    if flagh20sum == 1:
        for i in times:
            modelname1 = 'fchroma30_20b_5vt_'+str(int(i))+'s_H20.npz'
        
            #chosen model to compare (can/will be many)
            model_choice1 = np.load(base+modelname1)
            caiih_indsh20 = np.where((model_choice1['wl_rh']>396.7) & (model_choice1['wl_rh']< 397.94))
            
            model_choice1_wl = model_choice1['wl_rh'][caiih_indsh20]
            model_choice1_int = model_choice1['int_rh'][caiih_indsh20]
            model_choice1_wlshift = model_choice1_wl-lamb0
            model_choice1_wlshift_hep = model_choice1_wl-lamb1
            
            modelnameqs = 'fchroma30_50b_5vt_0s_H20.npz'
            
            model_choiceqs_wl = model_choiceqs['wl_rh'][caiih_inds]
            model_choiceqs_int = model_choiceqs['int_rh'][caiih_inds]
            model_choiceqs_wlshift = model_choiceqs_wl-lamb0
            model_choiceqs_wlshift_hep = model_choiceqs_wl-lamb1
                
            # adjust for instrument PSF
            
            model1_copy_int = model_choice1_int
            model1_copy_wl = model_choice1_wl
            yconvqs = psf_adjust(model_choiceqs_wl,model_choiceqs_int,fwhm,dkist_wl,ntw,gaussian_psf)
            
            yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                                           gaussian_psf)
            
            model_subtract = yconv1-yconvqs
            #model_subtract = yconv1
            
            models_tocomp.append(model_subtract)
            
    elif flagh20 == 1:
        for i in h20mods:
            
            if i == 25 or i == 31:
                modelname1 = 'fchroma'+str(int(i))+'_20b_5vt_43s_H20_2ndadjust.npz'
            else:
                modelname1 = 'fchroma'+str(int(i))+'_20b_5vt_43s_H20.npz'
        
            #chosen model to compare (can/will be many)
            model_choice1 = np.load(base+modelname1)
            caiih_indsh20 = np.where((model_choice1['wl_rh']>396.7) & (model_choice1['wl_rh']< 397.94))
            
            model_choice1_wl = model_choice1['wl_rh'][caiih_indsh20]
            model_choice1_int = model_choice1['int_rh'][caiih_indsh20]
            model_choice1_wlshift = model_choice1_wl-lamb0
            model_choice1_wlshift_hep = model_choice1_wl-lamb1
            
            modelnameqs = 'fchroma30_50b_5vt_0s_H20.npz'
            
            model_choiceqs_wl = model_choiceqs['wl_rh'][caiih_inds]
            model_choiceqs_int = model_choiceqs['int_rh'][caiih_inds]
            model_choiceqs_wlshift = model_choiceqs_wl-lamb0
            model_choiceqs_wlshift_hep = model_choiceqs_wl-lamb1
                
            # adjust for instrument PSF
            
            model1_copy_int = model_choice1_int
            model1_copy_wl = model_choice1_wl
            yconvqs = psf_adjust(model_choiceqs_wl,model_choiceqs_int,fwhm,dkist_wl,ntw,gaussian_psf)
            
            yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                                           gaussian_psf)
            
            model_subtract = yconv1-yconvqs
            #model_subtract = yconv1
            
            models_tocomp.append(model_subtract)
        
        modelname1 = 'fchroma30_50b_5vt_43s_H20.npz'
    
        #chosen model to compare (can/will be many)
        model_choice1 = np.load(base+modelname1)
        caiih_indsh20 = np.where((model_choice1['wl_rh']>396.7) & (model_choice1['wl_rh']< 397.94))
        
        model_choice1_wl = model_choice1['wl_rh'][caiih_indsh20]
        model_choice1_int = model_choice1['int_rh'][caiih_indsh20]
        model_choice1_wlshift = model_choice1_wl-lamb0
        model_choice1_wlshift_hep = model_choice1_wl-lamb1
        
        modelnameqs = 'fchroma30_50b_5vt_0s_H20.npz'
        
        model_choiceqs_wl = model_choiceqs['wl_rh'][caiih_inds]
        model_choiceqs_int = model_choiceqs['int_rh'][caiih_inds]
        model_choiceqs_wlshift = model_choiceqs_wl-lamb0
        model_choiceqs_wlshift_hep = model_choiceqs_wl-lamb1
            
        # adjust for instrument PSF
        
        model1_copy_int = model_choice1_int
        model1_copy_wl = model_choice1_wl
        yconvqs = psf_adjust(model_choiceqs_wl,model_choiceqs_int,fwhm,dkist_wl,ntw,gaussian_psf)
        
        yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                                       gaussian_psf)
        
        model_subtract = yconv1-yconvqs
        #model_subtract = yconv1
        h20mods.append(30)
        
        models_tocomp.append(model_subtract)
        
        modelname1 = 'fchroma30_20b_5vt_43s_H20_2ndadjust.npz'
    
        #chosen model to compare (can/will be many)
        model_choice1 = np.load(base+modelname1)
        caiih_indsh20 = np.where((model_choice1['wl_rh']>396.7) & (model_choice1['wl_rh']< 397.94))
        
        model_choice1_wl = model_choice1['wl_rh'][caiih_indsh20]
        model_choice1_int = model_choice1['int_rh'][caiih_indsh20]
        model_choice1_wlshift = model_choice1_wl-lamb0
        model_choice1_wlshift_hep = model_choice1_wl-lamb1
        
        modelnameqs = 'fchroma30_50b_5vt_0s_H20.npz'
        
        model_choiceqs_wl = model_choiceqs['wl_rh'][caiih_inds]
        model_choiceqs_int = model_choiceqs['int_rh'][caiih_inds]
        model_choiceqs_wlshift = model_choiceqs_wl-lamb0
        model_choiceqs_wlshift_hep = model_choiceqs_wl-lamb1
            
        # adjust for instrument PSF
        
        model1_copy_int = model_choice1_int
        model1_copy_wl = model_choice1_wl
        yconvqs = psf_adjust(model_choiceqs_wl,model_choiceqs_int,fwhm,dkist_wl,ntw,gaussian_psf)
        
        yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                                       gaussian_psf)
        
        model_subtract = yconv1-yconvqs
        #model_subtract = yconv1
        h20mods.append(30)
        
        models_tocomp.append(model_subtract)
        
        modelname1 = 'fchroma30_50b_5vt_21s_H20_2ndadjust.npz'
    
        #chosen model to compare (can/will be many)
        model_choice1 = np.load(base+modelname1)
        caiih_indsh20 = np.where((model_choice1['wl_rh']>396.7) & (model_choice1['wl_rh']< 397.94))
        
        model_choice1_wl = model_choice1['wl_rh'][caiih_indsh20]
        model_choice1_int = model_choice1['int_rh'][caiih_indsh20]
        model_choice1_wlshift = model_choice1_wl-lamb0
        model_choice1_wlshift_hep = model_choice1_wl-lamb1
        
        modelnameqs = 'fchroma30_50b_5vt_0s_H20.npz'
        
        model_choiceqs_wl = model_choiceqs['wl_rh'][caiih_inds]
        model_choiceqs_int = model_choiceqs['int_rh'][caiih_inds]
        model_choiceqs_wlshift = model_choiceqs_wl-lamb0
        model_choiceqs_wlshift_hep = model_choiceqs_wl-lamb1
            
        # adjust for instrument PSF
        
        model1_copy_int = model_choice1_int
        model1_copy_wl = model_choice1_wl
        yconvqs = psf_adjust(model_choiceqs_wl,model_choiceqs_int,fwhm,dkist_wl,ntw,gaussian_psf)
        
        yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                                       gaussian_psf)
        
        model_subtract = yconv1-yconvqs
        #model_subtract = yconv1
        h20mods.append(30)
        
        models_tocomp.append(model_subtract)
        
        modelname1 = 'fchroma30_1b_2vt_21s_H20_2ndadjust.npz'
    
        #chosen model to compare (can/will be many)
        model_choice1 = np.load(base+modelname1)
        caiih_indsh20 = np.where((model_choice1['wl_rh']>396.7) & (model_choice1['wl_rh']< 397.94))
        
        model_choice1_wl = model_choice1['wl_rh'][caiih_indsh20]
        model_choice1_int = model_choice1['int_rh'][caiih_indsh20]
        model_choice1_wlshift = model_choice1_wl-lamb0
        model_choice1_wlshift_hep = model_choice1_wl-lamb1
        
        modelnameqs = 'fchroma30_50b_5vt_0s_H20.npz'
        
        model_choiceqs_wl = model_choiceqs['wl_rh'][caiih_inds]
        model_choiceqs_int = model_choiceqs['int_rh'][caiih_inds]
        model_choiceqs_wlshift = model_choiceqs_wl-lamb0
        model_choiceqs_wlshift_hep = model_choiceqs_wl-lamb1
            
        # adjust for instrument PSF
        
        model1_copy_int = model_choice1_int
        model1_copy_wl = model_choice1_wl
        yconvqs = psf_adjust(model_choiceqs_wl,model_choiceqs_int,fwhm,dkist_wl,ntw,gaussian_psf)
        
        yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                                       gaussian_psf)
        
        model_subtract = yconv1-yconvqs
        #model_subtract = yconv1
        h20mods.append(30)
        
        models_tocomp.append(model_subtract)



fig,ax=plt.subplots(dpi=200,figsize=(5,4))
ncol2 = 20
map = tol_colors.tol_cmap(colormap='rainbow_discrete',lut=ncol2+4)
cmap_choice2 = map(np.linspace(0,1,ncol2+2))

ax.plot(dkist_wl-hepwl,normalize_range(model_subtract1/1e6,heplowh20,hephighh20),alpha=1,c='b',linewidth=2,label='EB1, 500s')

lns3 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[2][heplowh20:hephighh20]),linestyle=(0, (5, 1)),alpha=.7,label='EB2, 43s',c=cmap_choice2[12],linewidth=2)
lns4 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[3][heplowh20:hephighh20]),'-.',alpha=.7,label=r'EB3, 43s',c=cmap_choice2[6],linewidth=2)
lns6 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[5][heplowh20:hephighh20]),'-',alpha=.7,label=r'EB4, 43s',c='#CC6677',linewidth=2)
lns7 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[6][heplowh20:hephighh20]),'-',alpha=.7,label=r'EB4, 21s',c='#999933',linewidth=2)
#lns8 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[7][heplowh20:hephighh20]),'-.',alpha=.7,label=r'EB5, 21s',c='black',linewidth=2,zorder=7)
ax.plot(dkist_wl-hepwl,normalize_range(model_subtract_TC/1e6,heplowh20,hephighh20),alpha=1,c='red',linewidth=2,label='TC1, 90s')

ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(dkist_int[heplowh20:hephighh20]),label='DKIST/ViSP',linewidth=3,c='black',zorder=6,alpha=1)

# lns3 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[2][heplowh20:hephighh20]),linestyle=(0, (5, 1)),alpha=.7,label='EB2, 43s',c=cmap_choice2[12],linewidth=2)
# lns4 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[3][heplowh20:hephighh20]),'-.',alpha=.7,label=r'EB3, 43s',c=cmap_choice2[6],linewidth=2)
# lns6 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[5][heplowh20:hephighh20]),'-',alpha=.7,label=r'EB4, 43s',c='#CC6677',linewidth=2)
# lns7 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[6][heplowh20:hephighh20]),'-',alpha=.7,label=r'EB4, 21s',c='#999933',linewidth=2)
# #lns8 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[7][heplowh20:hephighh20]),'-.',alpha=.7,label=r'EB5, 21s',c='black',linewidth=2,zorder=7)
# ax.plot(dkist_wl-hepwl,normalize_range(model_subtract_TC/1e6,heplowh20,hephighh20),alpha=1,c='red',linewidth=2,label='TC1, 90s')

# ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(dkist_int[heplowh20:hephighh20]),label='DKIST/ViSP',linewidth=3,c='black',zorder=6,alpha=1)



ax.axvline(.0175,color=cmap_choice2[7],alpha=.7,linestyle='-.')
ax.axvline(.03,color=cmap_choice2[7],alpha=.7,linestyle='-.')
ax.axvline(.04,color=cmap_choice2[7],alpha=.7,linestyle='-.')
ax.axvline(-.04,color=cmap_choice2[7],alpha=.7,linestyle='-.')
ax.axvline(-.047,color=cmap_choice2[7],alpha=.7,linestyle='-.')


ax.legend(fontsize=8)
ax.set_ylabel('Normalized Intensity')
ax.grid()
ax.set_xlabel(r'Wavelength ($\lambda-\lambda_c$ [nm])')
#ax.set_xlabel(r'Wavelength')
secaxx = ax.secondary_xaxis('top', functions=(veltrans,wltrans))
secaxx.set_xlabel(r'Velocity $[km\; s^{-1}]$')
ax.set_xlim([-0.06,0.06])
ax.set_ylim([-.1,1.1])

fig.show()
