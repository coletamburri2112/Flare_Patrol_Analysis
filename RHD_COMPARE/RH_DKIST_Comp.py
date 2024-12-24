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
flagh20sum = 0 # this is the only change as of 10 Dec

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

flag = 'f-chroma' # longdur/f-chrom
if flag == 'f-chroma':
    #times = [1,5,7,17,19,21,25,30,35,36,37,38,39,40,41,42,43,44,46,47,48,49,50]
    #times = [1,5,7,17,19,21,25,30]
    #times = [25,30,35,37,40,43,44,47,50]
    #times= [37,40,43,44,47]
    times=[5,11,17,19,21,25,30,35,43]
elif flag == 'longdur':
    #times = [0,3,5,10,13,17,20,25,30,40,50,60,70,100]
    #times = [0,3,5,10,13,17,20]'
    times=[70]
    #times = [20,25,30,40,50,60,70,100]

broads = [1,20,50,70,100]
#broads=[1,50]
vturbs = [2,5,7,10]
h20mods = [79,19,25,31]
h20vals = [3e10,1e11,3e11,1e12]

if flagb ==1:
    ncol=len(broads)
elif flagvt ==1:
    ncol=len(vturbs)
elif flagt == 1:
    ncol = len(times)
elif flagh20 == 1:
    ncol = len(h20mods)+3


ncol=len(times)
map = tol_colors.tol_cmap(colormap='rainbow_discrete',lut=ncol+5)
cmap_choice = map(np.linspace(0,1,ncol))

base = '/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_output_files_npz/'


# hep_inds = np.where((rh_file3['wl_rh']>396.95) & (rh_file3['wl_rh']< 397.07))

filename_dkist = '/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_output_files_npz/dkist_obs_file.npz'
if flag == 'f-chroma':
    modelnameqs = 'fchroma30_50b_5vt_0s_H6.npz'
elif flag == 'longdur':
    modelnameqs = 'longduration_0s_H6.npz'

dkist_file = np.load(filename_dkist)



model_choiceqs = np.load(base+modelnameqs)   

caiih_inds = np.where((model_choiceqs['wl_rh']>396.7) & (model_choiceqs['wl_rh']< 397.94))

model_choiceqs_wl = model_choiceqs['wl_rh'][caiih_inds]
model_choiceqs_int = model_choiceqs['int_rh'][caiih_inds]
model_choiceqs_wlshift = model_choiceqs_wl-lamb0
model_choiceqs_wlshift_hep = model_choiceqs_wl-lamb1

maxind = dkist_file['arr_3'][0]

dkist_wl = dkist_file['arr_2']
dkist_int = dkist_file['arr_0'][0,:,maxind] #before background subtraction
dkist_int_trail = dkist_file['arr_0'][0,:,maxind+120]
dkist_int_front = dkist_file['arr_0'][0,:,maxind-100]

dkist_wl_shift = dkist_wl-lamb0
dkist_wl_shift_hep = dkist_wl-lamb1

yconvqs = psf_adjust(model_choiceqs_wl,model_choiceqs_int,fwhm,dkist_wl,ntw,gaussian_psf)

models_tocomp = []


if flag =='f-chroma':
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
        
        modelname1 = 'fchroma30_20b_5vt_43s_H20.npz'
    
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
        
        modelname1 = 'fchroma30_50b_5vt_21s_H20.npz'
    
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
        
        modelname1 = 'fchroma30_1b_2vt_21s_H20.npz'
    
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
            
    else:
        if flagb == 1:
            for i in broads:
                modelname1 = 'fchroma30_'+str(int(i))+'b_5vt_43s_H6.npz'
                #chosen model to compare (can/will be many)
                model_choice1 = np.load(base+modelname1)
                
                model_choice1_wl = model_choice1['wl_rh'][caiih_inds]
                model_choice1_int = model_choice1['int_rh'][caiih_inds]
                model_choice1_wlshift = model_choice1_wl-lamb0
                model_choice1_wlshift_hep = model_choice1_wl-lamb1
                    
                # adjust for instrument PSF
                
                model1_copy_int = model_choice1_int
                model1_copy_wl = model_choice1_wl
                
                yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                                               gaussian_psf)
                
                model_subtract = yconv1-yconvqs
                
                models_tocomp.append(model_subtract)
                
            modelname1 = 'fchroma30_1b_2vt_43s.npz'
            
            model_choice1 = np.load(base+modelname1)
            
            model_choice1_wl = model_choice1['wl_rh'][caiih_inds]
            model_choice1_int = model_choice1['int_rh'][caiih_inds]
            model_choice1_wlshift = model_choice1_wl-lamb0
            model_choice1_wlshift_hep = model_choice1_wl-lamb1
                
            # adjust for instrument PSF
            
            model1_copy_int = model_choice1_int
            model1_copy_wl = model_choice1_wl
            
            yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                                           gaussian_psf)
            
            model_subtract = yconv1-yconvqs
            
            models_tocomp.append(model_subtract)
                
        elif flagvt == 1:
            for i in vturbs:
                modelname1 = 'fchroma30_20b_'+str(int(i))+'vt_43s.npz'
                #chosen model to compare (can/will be many)
                model_choice1 = np.load(base+modelname1)
                
                model_choice1_wl = model_choice1['wl_rh'][caiih_inds]
                model_choice1_int = model_choice1['int_rh'][caiih_inds]
                model_choice1_wlshift = model_choice1_wl-lamb0
                model_choice1_wlshift_hep = model_choice1_wl-lamb1
                    
                # adjust for instrument PSF
                
                model1_copy_int = model_choice1_int
                model1_copy_wl = model_choice1_wl
                
                yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                                               gaussian_psf)
                
                model_subtract = yconv1-yconvqs
                
                models_tocomp.append(model_subtract)        
                
                
            modelname1 = 'fchroma30_1b_2vt_43s.npz'
            
            model_choice1 = np.load(base+modelname1)
            
            model_choice1_wl = model_choice1['wl_rh'][caiih_inds]
            model_choice1_int = model_choice1['int_rh'][caiih_inds]
            model_choice1_wlshift = model_choice1_wl-lamb0
            model_choice1_wlshift_hep = model_choice1_wl-lamb1
                
            # adjust for instrument PSF
            
            model1_copy_int = model_choice1_int
            model1_copy_wl = model_choice1_wl
            
            yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                                           gaussian_psf)
            
            model_subtract = yconv1-yconvqs
            
            models_tocomp.append(model_subtract)
        else: 
            for i in times:
                print('here2')
                
                if i < 35.5:
                    modelname1 = 'rhf1d_fchroma_'+str(int(i))+'s_30_20broadc_5vt.npz'
                
                    #chosen model to compare (can/will be many)
                    model_choice1 = np.load(base+modelname1)
                    
                    model_choice1_wl = model_choice1['wl_rh'][caiih_inds]
                    model_choice1_int = model_choice1['int_rh'][caiih_inds]
                    model_choice1_wlshift = model_choice1_wl-lamb0
                    model_choice1_wlshift_hep = model_choice1_wl-lamb1
                        
                    # adjust for instrument PSF
                    
                    model1_copy_int = model_choice1_int
                    model1_copy_wl = model_choice1_wl
                    
                    yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                                                   gaussian_psf)
                    
                    model_subtract = yconv1-yconvqs
                    
                    models_tocomp.append(model_subtract)
                elif i > 35.5:
                    #modelname1 = 'fchroma30_20b_5vt_'+str(int(i))+'s_H6.npz'
                    modelname2 = 'fchroma30_50b_5vt_'+str(int(i))+'s_H6.npz'
                    #modelname3 = 'fchroma30_70b_5vt_'+str(int(i))+'s_H6.npz'
                
                    #chosen model to compare (can/will be many)
                    #model_choice1 = np.load(base+modelname1)
                    model_choice2 = np.load(base+modelname2)
                    #model_choice3 = np.load(base+modelname3)
                    
                    #model_choice1_wl = model_choice1['wl_rh'][caiih_inds]
                    #model_choice1_int = model_choice1['int_rh'][caiih_inds]
                    #model_choice1_wlshift = model_choice1_wl-lamb0
                    #model_choice1_wlshift_hep = model_choice1_wl-lamb1
     
                    model_choice2_wl = model_choice2['wl_rh'][caiih_inds]
                    model_choice2_int = model_choice2['int_rh'][caiih_inds]
                    model_choice2_wlshift = model_choice2_wl-lamb0
                    model_choice2_wlshift_hep = model_choice2_wl-lamb1                   
                    # adjust for instrument PSF
                    
                    #model1_copy_int = model_choice1_int
                    #model1_copy_wl = model_choice1_wl
                    
                    #yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                                                   #gaussian_psf)
                    
                    model2_copy_int = model_choice2_int
                    model2_copy_wl = model_choice2_wl
                    
                    yconv2=psf_adjust(model2_copy_wl,model2_copy_int,fwhm,dkist_wl,ntw,
                                                   gaussian_psf)
                    
                    #model_choice3_wl = model_choice3['wl_rh'][caiih_inds]
                    #model_choice3_int = model_choice3['int_rh'][caiih_inds]
                    #model_choice3_wlshift = model_choice3_wl-lamb0
                    #model_choice3_wlshift_hep = model_choice3_wl-lamb1 
    
                    #model3_copy_int = model_choice3_int
                    #model3_copy_wl = model_choice3_wl  
    
                    #yconv3=psf_adjust(model3_copy_wl,model3_copy_int,fwhm,dkist_wl,ntw,
                    #                              gaussian_psf)                
                    # adjust for instrument PSF
                    
                   # model_subtract1 = yconv1-yconvqs
                    model_subtract2 = yconv2-yconvqs
                   # model_subtract3 = yconv3-yconvqs
                    
                    #models_tocomp.append(model_subtract)
                   # models_tocomp.append(model_subtract1)
                    models_tocomp.append(model_subtract2)
                   # models_tocomp.append(model_subtract3)

if flag == 'longdur':
    for i in times:
        modelname1 = 'longduration_'+str(int(i))+'s_H6.npz'
    
        #chosen model to compare (can/will be many)
        model_choice1 = np.load(base+modelname1)
        
        model_choice1_wl = model_choice1['wl_rh'][caiih_inds]
        model_choice1_int = model_choice1['int_rh'][caiih_inds]
        model_choice1_wlshift = model_choice1_wl-lamb0
        model_choice1_wlshift_hep = model_choice1_wl-lamb1
            
        # adjust for instrument PSF
        
        model1_copy_int = model_choice1_int
        model1_copy_wl = model_choice1_wl
        
        yconv1=psf_adjust(model1_copy_wl,model1_copy_int,fwhm,dkist_wl,ntw,
                                       gaussian_psf)
        
        model_subtract = yconv1-yconvqs
        
        models_tocomp.append(model_subtract)
        
        # added for new broadc
        modelname2 = 'longduration_'+str(int(i))+'s_H6_1broadc.npz'
    
        #chosen model to compare (can/will be many)
        model_choice2 = np.load(base+modelname2)
        
        model_choice2_wl = model_choice2['wl_rh'][caiih_inds]
        model_choice2_int = model_choice2['int_rh'][caiih_inds]
        model_choice2_wlshift = model_choice2_wl-lamb0
        model_choice2_wlshift_hep = model_choice2_wl-lamb1
            
        # adjust for instrument PSF
        
        model2_copy_int = model_choice2_int
        model2_copy_wl = model_choice2_wl
        
        yconv2=psf_adjust(model2_copy_wl,model2_copy_int,fwhm,dkist_wl,ntw,
                                       gaussian_psf)
        
        model_subtract2 = yconv2-yconvqs
        
        models_tocomp.append(model_subtract2)
        

if flagh20sum == 1:
    fig,ax=plt.subplots(dpi=200)
    
    ncol = len(times)
    map = tol_colors.tol_cmap(colormap='rainbow_discrete',lut=ncol+4)
    cmap_choice = map(np.linspace(0,1,ncol+2))
    markers = ['-.','-','--','-.','-','-']
    
    ncol2 = 20
    map = tol_colors.tol_cmap(colormap='rainbow_discrete',lut=ncol2+4)
    cmap_choice2 = map(np.linspace(0,1,ncol2+2))
    
    for i in range(len(models_tocomp)):
        #ax.plot(dkist_wl,models_tocomp[i]/1e6,alpha=1,label=r't3F10-15-8/20w-5v, '+str(int(times[i]))+'s',c=cmap_choice[-i],linewidth=3)
        ax.plot(dkist_wl,models_tocomp[i]/1e6,alpha=1,label=str(int(times[i]))+'s',c=cmap_choice[-i],linewidth=3)
        ax.legend()
        ax.set_xlim([396.75,397.1])
        ax.grid()
        ax.axvline(396.844)
        ax.axvline(397.01)
    ax.set_ylabel(r'Intensity $[10^6\;erg\;s^{-1}\;cm^{-2}\;\AA^{-1}\;sr^{-1}]$')
    ax.set_xlabel('Wavelength [nm]')
    ax.set_xticks([396.8,396.9,397.0,397.1])
        
    
else:
    #ax1 = ax.twinx()
    fig,ax=plt.subplots(dpi=200)
    ncol = 4
    map = tol_colors.tol_cmap(colormap='rainbow_discrete',lut=ncol+4)
    cmap_choice = map(np.linspace(0,1,ncol+2))
    markers = ['-.','-','--','-.','-','-','-','-','-','-']
    
    ncol2 = 20
    map = tol_colors.tol_cmap(colormap='rainbow_discrete',lut=ncol2+4)
    cmap_choice2 = map(np.linspace(0,1,ncol2+2))
    
    corrmodsbvar = [11,5,6,12,13]
    corrmodsvtvar = [8,5,9,10]
    #for i in range(10):
    if flag == 'f-chroma':
        for i in range(len(models_tocomp)-1):
            if flagb == 1:
                if i==1:
                    #lns2 = ax.plot(dkist_wl-396.846,models_tocomp[i]/1e6,markers[i],alpha=1,label=r't3F10-15-8/'+str(broads[i])+'w-5v, 43s',c=cmap_choice[-i],linewidth=3,zorder=5)
                    lns2 = ax.plot(dkist_wl-396.846,models_tocomp[i]/1e6,markers[i],alpha=1,label='Model '+str(corrmodsbvar[i]),c=cmap_choice[-i],linewidth=3,zorder=5)
                else:
                    #lns2 = ax.plot(dkist_wl-396.846,models_tocomp[i]/1e6,markers[i],alpha=1,label=r't3F10-15-8/'+str(broads[i])+'w-5v, 43s',c=cmap_choice[-i],linewidth=3)
                    lns2 = ax.plot(dkist_wl-396.846,models_tocomp[i]/1e6,markers[i],alpha=1,label='Model '+str(corrmodsbvar[i]),c=cmap_choice[-i],linewidth=3,zorder=5)
    
    
    
            elif flagvt == 1:
                lns2 = ax.plot(dkist_wl-396.846,models_tocomp[i]/1e6,markers[i],alpha=1,label='Model '+str(corrmodsvtvar[i]),c=cmap_choice[-i],linewidth=3)
                
        if flagb ==1 or flagvt ==1:
            lns2 = ax.plot(dkist_wl-396.846,models_tocomp[len(models_tocomp)-1]/1e6,'-',alpha=.5,label='Model 7',c='#999933',linewidth=3,zorder=1)
    
        #for h20 comparisons
        #lns1 = ax.plot(dkist_wl,normalize_range(models_tocomp[0]/1e6,heplowh20,hephighh20),'--',alpha=.9,label=r't3F9-15-3/20w-5v, 43s',c=cmap_choice2[0],linewidth=2)
        #lns2 = ax.plot(dkist_wl,normalize_range(models_tocomp[1]/1e6,heplowh20,hephighh20),'--',alpha=.9,label=r'tF10-15-3/20w-5v, 43s',c=cmap_choice2[4],linewidth=2)
        # lns3 = ax.plot(dkist_wl,normalize_range(models_tocomp[2]/1e6,heplowh20,hephighh20),linestyle=(0, (5, 1)),alpha=.9,label=r't3F10-15-3/20w-5v, 43s',c=cmap_choice2[12],linewidth=2)
        # lns4 = ax.plot(dkist_wl,normalize_range(models_tocomp[3]/1e6,heplowh20,hephighh20),'-.',alpha=.9,label=r'tF11-15-3/20w-5v, 43s',c=cmap_choice2[6],linewidth=2)
        # # #lns5 = ax.plot(dkist_wl,models_tocomp[4],'-',alpha=.9,label=r't3F10-15-8/50w-5v, 43s',c=cmap_choice[0],linewidth=2)
        # lns6 = ax.plot(dkist_wl,normalize_range(models_tocomp[5]/1e6,heplowh20,hephighh20),'-*',markersize=1,alpha=.8,label=r't3F10-15-8/20w-5v, 43s',c='red',linewidth=2)
        # lns7 = ax.plot(dkist_wl,normalize_range(models_tocomp[6]/1e6,heplowh20,hephighh20),'-',alpha=.9,label=r't3F10-15-8/50w-5v, 21s',c='#999933',linewidth=2,zorder=7)
        # lns8 = ax.plot(dkist_wl,normalize_range(models_tocomp[7]/1e6,heplowh20,hephighh20),'-.',alpha=.9,label=r't3F10-15-8/1w-2v, 21s',c='#AA3377',linewidth=2,zorder=8)

        lns3 = ax.plot(dkist_wl,normalize_range(models_tocomp[2]/1e6,heplowh20,hephighh20),linestyle=(0, (5, 1)),alpha=.9,label=r'Model 3, 43s',c=cmap_choice2[12],linewidth=2)
        lns4 = ax.plot(dkist_wl,normalize_range(models_tocomp[3]/1e6,heplowh20,hephighh20),'-.',alpha=.9,label=r'Model 4, 43s ',c=cmap_choice2[6],linewidth=2)
        # #lns5 = ax.plot(dkist_wl,models_tocomp[4],'-',alpha=.9,label=r't3F10-15-8/50w-5v, 43s',c=cmap_choice[0],linewidth=2)
        lns6 = ax.plot(dkist_wl,normalize_range(models_tocomp[5]/1e6,heplowh20,hephighh20),'-*',markersize=1,alpha=.8,label=r'Model 5, 43s',c='red',linewidth=2)
        lns7 = ax.plot(dkist_wl,normalize_range(models_tocomp[6]/1e6,heplowh20,hephighh20),'-',alpha=.9,label=r'Model 6, 21s',c='#999933',linewidth=2,zorder=7)
        lns8 = ax.plot(dkist_wl,normalize_range(models_tocomp[7]/1e6,heplowh20,hephighh20),'-.',alpha=.9,label=r'Model 7, 21s',c='#AA3377',linewidth=2,zorder=8)
        
        # # for no normalization with H20 comparison
        #lns1 = ax.plot(dkist_wl,normalize_range(models_tocomp[0]/1e6,heplowh20,hephighh20),'--',alpha=.9,label=r't3F9-15-3/20w-5v, 43s',c=cmap_choice2[0],linewidth=2)
        # #lns2 = ax.plot(dkist_wl,normalize_range(models_tocomp[1]/1e6,heplowh20,hephighh20),'--',alpha=.9,label=r'tF10-15-3/20w-5v, 43s',c=cmap_choice2[4],linewidth=2)
        # lns3 = ax.plot(dkist_wl,models_tocomp[2]/1e6,linestyle=(0, (5, 1)),alpha=.9,label=r't3F10-15-3/20w-5v, 43s',c=cmap_choice2[12],linewidth=2)
        # lns4 = ax.plot(dkist_wl,models_tocomp[3]/1e6,'-.',alpha=.9,label=r'tF11-15-3/20w-5v, 43s',c=cmap_choice2[6],linewidth=2)
        # # #lns5 = ax.plot(dkist_wl,models_tocomp[4],'-',alpha=.9,label=r't3F10-15-8/50w-5v, 43s',c=cmap_choice[0],linewidth=2)
        # lns6 = ax.plot(dkist_wl,models_tocomp[5]/1e6,'-*',markersize=1,alpha=.8,label=r't3F10-15-8/20w-5v, 43s',c='red',linewidth=2)
        # lns7 = ax.plot(dkist_wl,models_tocomp[6]/1e6,'-',alpha=.9,label=r't3F10-15-8/50w-5v, 21s',c='#999933',linewidth=2,zorder=7)
        # lns8 = ax.plot(dkist_wl,models_tocomp[7]/1e6,'-.',alpha=.9,label=r't3F10-15-8/1w-2v, 21s',c='#AA3377',linewidth=2,zorder=8)
        
    
     
        #ax.set_xlim([-.1,.1])
        #ax.set_xlim([396.7,397.15])
    elif flag == 'longdur':
        #for i in range(len(models_tocomp)):
        #lns2 = ax.plot(dkist_wl,models_tocomp[i],'--',alpha=1,label=r'long_duration, '+str(times[i])+'s',c=cmap_choice[i],linewidth=1)
        lns2 = ax.plot(dkist_wl,models_tocomp[0],'--',alpha=1,label=r'long_duration, 70s, broadc = 50',c=cmap_choice[0],linewidth=1)
        lns2 = ax.plot(dkist_wl,models_tocomp[1],'--',alpha=1,label=r'long_duration, 70s, broadc = 1',c=cmap_choice[1],linewidth=1)
    
    ax.set_xlim([396.75,396.95])
    
    fig.show()
    # if flagh20 ==1:
    #     lns2 = ax.plot(dkist_wl,models_tocomp[4]),'--',alpha=1,label=r'F-CHROMA25/50w-5v, 43s',c=cmap_choice[4],linewidth=1)
    ax.set_xticks([396.7,396.8,396.9,397,397.1,397.2])
    #ax.set_xticks([-.08,-0.04,0,0.04,0.08])
    
    #     ax1.plot(dkist_wl[heplowh20:hephighh20],normalize(dkist_int[heplowh20:hephighh20]),label='ViSP, ribbon center',linewidth=1,c='red',zorder=2)
    #ax.set_xticks([396.76,396.8,396.84,396.88,396.92])
    ax.plot(dkist_wl,normalize_range(dkist_int/1e6,heplowh20,hephighh20),label='ViSP, ribbon center',linewidth=3,c='black',zorder=6,alpha=1)
    ax.set_ylim([-0.5,7])
    ax.set_xlim([396.7,397.1])
    #ax.plot(dkist_wl,dkist_int/1e6,label='ViSP, ribbon center',linewidth=3,c='black',zorder=6,alpha=1)
    
    ax.legend(fontsize=8)
    #ax.set_xticks([-0.08,-.04,0,.04,.08])
    
    ax.set_ylabel(r'Intensity $[10^6\;erg\;s^{-1}\;cm^{-2}\;\AA^{-1}\;sr^{-1}]$')
    ax.set_ylabel('Intensity (Normalized to H$\epsilon$)')
    # if flagh20 == 1:
    #     ax.set_ylabel('Normalized Intensity')
    #ax1.set_ylabel('DKIST Intensity',color='red')
    #ax1.set_ylim([-.25e6,4.5e6])
    # ax.set_ylim([-.5,5.5])
    #ax.axvline(396.846,linestyle='--',c='grey')#
    #ax.axvline(397.01,linestyle='--',c='grey')
    #ax.set_xlim([396.7,397.1])
    ax.set_xticks([396.75,396.85,396.95,397.05])
    ax.set_xlabel('Wavelength [nm]')
    #ax.axvline(397.01)
    # ax.set_xlabel(r'Wavelength ($\lambda-\lambda_c$) [nm]')
    # secaxx = ax.secondary_xaxis('top', functions=(veltrans,wltrans))
    # secaxx.set_xlabel(r'Velocity $[km\; s^{-1}]$')
    # secaxx.set_xticks([-120,-80,-40,0,40,80,120])
    # secaxx.grid('on')
    ax.grid()
    hepflag=0
    hepwl = 397.01
    if hepflag ==1:      
        fig,ax=plt.subplots(dpi=200,figsize=[20,10])
        # #lns1 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[0][heplowh20:hephighh20]),'--',alpha=.9,label=r't3F9-15-3/20w-5v, 43s',c=cmap_choice2[0],linewidth=3)
        # #lns2 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[1][heplowh20:hephighh20]),'--',alpha=.9,label=r't1F10-15-3/20w-5v, 43s',c=cmap_choice2[4],linewidth=3)
        # lns3 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[2][heplowh20:hephighh20]),linestyle=(0, (5, 1)),alpha=.9,label=r't3F10-15-3/20w-5v, 43s',c=cmap_choice2[12],linewidth=3)
        # lns4 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[3][heplowh20:hephighh20]),'-.',alpha=.9,label=r'tF11-15-3/20w-5v, 43s',c=cmap_choice2[6],linewidth=3)
        # #lns5 = ax.plot(dkist_wl[heplowh20:hephighh20],normalize(models_tocomp[4][heplowh20:hephighh20]),'-',alpha=.9,label=r't3F10-15-8/50w-5v, 43s',c=cmap_choice[0],linewidth=2)
        # lns6 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[5][heplowh20:hephighh20]),'-',alpha=.9,label=r't3F10-15-8/20w-5v, 43s',c='red',linewidth=3)
        # lns7 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[6][heplowh20:hephighh20]),'-',alpha=.9,label=r't3F10-15-8/50w-5v, 21s',c='#999933',linewidth=3)
        # lns8 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[7][heplowh20:hephighh20]),'-.',alpha=.9,label=r't3F10-15-8/1w-2v, 21s',c='#AA3377',linewidth=3,zorder=7)

        # #lns1 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[0][heplowh20:hephighh20]),'--',alpha=.9,label=r't3F9-15-3/20w-5v, 43s',c=cmap_choice2[0],linewidth=3)
        # #lns2 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[1][heplowh20:hephighh20]),'--',alpha=.9,label=r't1F10-15-3/20w-5v, 43s',c=cmap_choice2[4],linewidth=3)
        # lns3 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[2][heplowh20:hephighh20]),linestyle=(0, (5, 1)),alpha=.9,label='Model 3, 43s',c=cmap_choice2[12],linewidth=3)
        # lns4 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[3][heplowh20:hephighh20]),'-.',alpha=.9,label=r'Model 4, 43s',c=cmap_choice2[6],linewidth=3)
        # #lns5 = ax.plot(dkist_wl[heplowh20:hephighh20],normalize(models_tocomp[4][heplowh20:hephighh20]),'-',alpha=.9,label=r't3F10-15-8/50w-5v, 43s',c=cmap_choice[0],linewidth=2)
        # lns6 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[5][heplowh20:hephighh20]),'-',alpha=.9,label=r'Model 5, 43s',c='red',linewidth=3)
        # lns7 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[6][heplowh20:hephighh20]),'-',alpha=.9,label=r'Model 6, 21s',c='#999933',linewidth=3)
        # lns8 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[7][heplowh20:hephighh20]),'-.',alpha=.9,label=r'Model 7, 21s',c='#AA3377',linewidth=3,zorder=7)
    
        # ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(dkist_int[heplowh20:hephighh20]),label='ViSP, ribbon center',linewidth=3,c='black',zorder=6,alpha=1)

        #lns1 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[0][heplowh20:hephighh20]),'--',alpha=.9,label=r't3F9-15-3/20w-5v, 43s',c=cmap_choice2[0],linewidth=3)
        #lns2 = ax.plot(dkist_wl[heplowh20:hephighh20]-hepwl,normalize(models_tocomp[1][heplowh20:hephighh20]),'--',alpha=.9,label=r't1F10-15-3/20w-5v, 43s',c=cmap_choice2[4],linewidth=3)
        lns3 = ax.plot(dkist_wl,normalize(models_tocomp[2]),linestyle=(0, (5, 1)),alpha=.9,label='Model 3, 43s',c=cmap_choice2[12],linewidth=3)
        lns4 = ax.plot(dkist_wl,normalize(models_tocomp[3]),'-.',alpha=.9,label=r'Model 4, 43s',c=cmap_choice2[6],linewidth=3)
        #lns5 = ax.plot(dkist_wl[heplowh20:hephighh20],normalize(models_tocomp[4][heplowh20:hephighh20]),'-',alpha=.9,label=r't3F10-15-8/50w-5v, 43s',c=cmap_choice[0],linewidth=2)
        lns6 = ax.plot(dkist_wl,normalize(models_tocomp[5]),'-',alpha=.9,label=r'Model 5, 43s',c='red',linewidth=3)
        lns7 = ax.plot(dkist_wl,normalize(models_tocomp[6]),'-',alpha=.9,label=r'Model 6, 21s',c='#999933',linewidth=3)
        lns8 = ax.plot(dkist_wl,normalize(models_tocomp[7]),'-.',alpha=.9,label=r'Model 7, 21s',c='#AA3377',linewidth=3,zorder=7)
    
        ax.plot(dkist_wl,normalize(dkist_int),label='ViSP, ribbon center',linewidth=3,c='black',zorder=6,alpha=1)
                
        
        # ax.axvline(.0175,color=cmap_choice2[7],alpha=.7,linestyle='-.')
        # ax.axvline(.03,color=cmap_choice2[7],alpha=.7,linestyle='-.')
        # ax.axvline(.04,color=cmap_choice2[7],alpha=.7,linestyle='-.')
        # ax.axvline(-.04,color=cmap_choice2[7],alpha=.7,linestyle='-.')
        # ax.axvline(-.047,color=cmap_choice2[7],alpha=.7,linestyle='-.')
        
        ax.legend(fontsize=8)
        ax.set_ylabel('Normalized Intensity')
        ax.grid()
        #ax.set_xlabel(r'Wavelength ($\lambda-\lambda_c$ [nm])')
        ax.set_xlabel(r'Wavelength')
        #secaxx = ax.secondary_xaxis('top', functions=(veltrans,wltrans))
        #secaxx.set_xlabel(r'Velocity $[km\; s^{-1}]$')
        ax.set_xlim([396.72,397.1])
        
        
    
    
    
    
    
    
    
    