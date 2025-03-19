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
matplotlib.use('qtagg')

# Function definitions for Gaussian fitting
def Gauss_func(x,A,mu,sigma,m,b):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))+ m*x + b

def double_gaussian( x, c1, mu1, sigma1, c2, mu2, sigma2 ,m,b):
    res =   (c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) )) \
          + (c2 * np.exp( - (x - mu2)**2.0 / (2.0 * sigma2**2.0) )) \
          + (m * x + b)
    return res
        
#### USER-DEFINED FUNCTIONS
filename = '/Users/coletamburri/Desktop/August_2024_DKIST_Flares/postdestretch_dataCubeX_class_decay_full.fits' #folder to save output images, array to


# Define H-alpha file (all data) and frame to work with
fullhalpha = fits.open(filename)
#fullhalpha = fits.open(path+folder_vbi+'/'+filename)
array = fullhalpha[0].data[0:50,:,:]
gbref=array[0,:,:]
timegb,yy,xx=array.shape
   
img1 = array[0,:,:]

# #hi-pass filter?

# F1=fp.fft2((img1).astype(float))
# F2 = fp.fftshift(F1)

# # plt.figure(figsize=(10,10))
# # plt.imshow( (20*np.log10( 0.1 + F2)).astype(int), cmap=plt.cm.gray)
# # plt.show()

# (w, h) = img1.shape
# half_w, half_h = int(w/2), int(h/2)

# # high pass filter
# n = 10
# F2[half_w-n:half_w+n+1,half_h-n:half_h+n+1] = 0 # select all but the first 50x50 (low) frequencies
# # plt.figure(figsize=(10,10))
# # plt.imshow( (20*np.log10( 0.1 + F2)).astype(int))
# # plt.show()

# im1 = fp.ifft2(fp.ifftshift(F2)).real

#masking of loops:
    
ff = []

# levs3 = np.arange(1e3,9e3,1e3)
# for k in levs3:
#     img2 = np.zeros(np.shape(img1))
#     for i in range(np.shape(img1)[0]):
#         for j in range(np.shape(img1)[1]):
#             if img1[i,j]<k:
#                 img2[i,j]=-1
    
#     ff.append(np.sum(img2==-1)/(np.shape(img2)[0]*np.shape(img2)[1]))

levs4 = np.arange(1e4,9e4,.1e4)
for k in levs4:
    img2 = np.zeros(np.shape(img1))
    for i in range(np.shape(img1)[0]):
        for j in range(np.shape(img1)[1]):
            if img1[i,j]<k:
                img2[i,j]=-1
    
    ff.append(np.sum(img2==-1)/(np.shape(img2)[0]*np.shape(img2)[1]))

# levs5 = np.arange(1e5,5e5,1e5)
# for k in levs5:
#     img2 = np.zeros(np.shape(img1))
#     for i in range(np.shape(img1)[0]):
#         for j in range(np.shape(img1)[1]):
#             if img1[i,j]<k:
#                 img2[i,j]=-1
    
#     ff.append(np.sum(img2==-1)/(np.shape(img2)[0]*np.shape(img2)[1]))
    
# fig,ax=plt.subplots()
# ax.imshow(img1)
# fig.show()
