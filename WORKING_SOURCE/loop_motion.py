#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 19:25:26 2025

@author: coletamburri
"""


import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import skimage
import scipy
import tol_colors as tc
import matplotlib
import utilvbi

from scipy.signal import convolve2d
from scipy.signal import convolve
from scipy.ndimage import gaussian_filter
from skimage.exposure import match_histograms

from skimage.io import imread
import matplotlib.pyplot as plt
import scipy.fftpack as fp

# allow for pop-up... need pyqt5 installed
matplotlib.use('inline')

# Function definitions for Gaussian fitting
def Gauss_func(x,A,mu,sigma,m,b):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))+ m*x + b

def double_gaussian( x, c1, mu1, sigma1, c2, mu2, sigma2 ,m,b):
    res =   (c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) )) \
          + (c2 * np.exp( - (x - mu2)**2.0 / (2.0 * sigma2**2.0) )) \
          + (m * x + b)
    return res
        
#### USER-DEFINED FUNCTIONS
#filename = '/Users/coletamburri/Desktop/August_2024_DKIST_Flares/postdestretch_dataCubeX_class_decay_full.fits' #folder to save output images, array to
filename = '/Users/coletamburri/Desktop/VBI_Destretching/AXXJL/postdestretch_histomatch_dataCubeX_class_decay_coronal_rain.fits'

# Define H-alpha file (all data) and frame to work with
fullhalpha = fits.open(filename)
#fullhalpha = fits.open(path+folder_vbi+'/'+filename)
array = fullhalpha[0].data[:,:,:]
# gbref=array[0,:,:]
# timegb,yy,xx=array.shape
   
# dataCubeTrackedhist=[]
# for idata in range(timegb):
#     image=array[idata,:,:]
#     matched = match_histograms(image, gbref)
#     dataCubeTrackedhist.append(matched)

# dcarr = np.asarray(dataCubeTrackedhist)

# high pass filter (remove low-freq components; seems good for ribbon front, actually)
img1=array[2,:,:]
# F1=fp.fft2((img1).astype(float))
# F2 = fp.fftshift(F1)

# plt.figure(figsize=(10,10))
# plt.imshow( (20*np.log10( 0.1 + F2)).astype(int), cmap=plt.cm.gray)
# plt.show()

# (w, h) = img1.shape
# half_w, half_h = int(w/2), int(h/2)

# # high pass filter
# n = 50
# F2[half_w-n:half_w+n+1,half_h-n:half_h+n+1] = 0 # select all but the first 50x50 (low) frequencies
# plt.figure(figsize=(10,10))
# plt.imshow( (20*np.log10( 0.1 + F2)).astype(int))
# plt.show()

# im1 = fp.ifft2(fp.ifftshift(F2)).real
# plt.figure(figsize=(10,10))
# plt.imshow(im1, cmap='gray')
# plt.axis('off')
# plt.show()


diff=[]

normmin = np.min(img1[2000:3000,2000:3000])
normmax = np.max(img1[2000:3000,2000:3000])

for i in range(np.shape(array)[0]):
    #frame to work with
    frame = array[i,2000:3000,2000:3000]
    framelast = array[i-1,2000:3000,2000:3000]
    
    framenorm = (frame-normmin)/(normmax-normmin)
    framelastnorm = (framelast-normmin)/(normmax-normmin)
    
    
    slicefr = framenorm[:]
    slicefr_last = framelastnorm[:]
    diff.append(np.subtract(slicefr,slicefr_last))
    
utilvbi.storeSequence(np.asarray(diff[1:]),'/Users/coletamburri/Desktop/diff.mp4', dpi=300, write=True)
utilvbi.storeSequence(array[:,2000:3000,2000:3000],'/Users/coletamburri/Desktop/no_diff.mp4', dpi=300, write=True)

