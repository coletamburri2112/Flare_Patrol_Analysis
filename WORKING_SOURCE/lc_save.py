#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 10:31:50 2025

@author: coletamburri
"""

import astropy.units as u
from sunpy.net import Fido, attrs as a
import dkist.net
import dkist
import utilvbi
import destretch_fw as destr

from astropy.io import fits
import numpy as np
from IPython.display import Video
from astropy.io import fits
import matplotlib.pyplot as plt
import os

import astropy.io.fits as pf
import utilvbi
import skimage
from skimage.exposure import match_histograms

from matplotlib.patches import Rectangle

pid='2_11'
vbiexp = 'AXXJL' # 8 August 2024 X-class

savfold='/Users/coletamburri/Desktop/VBI_Destretching/'+vbiexp+'/'
filt='Halpha'

dataCube=fits.open(savfold+'postdestretch_histomatch_dataCubeX_class_decay_looptopbright.fits')[0].data

gbref=dataCube[0,:,:] #choosing here the first image as a reference
timegb,yy,xx=dataCube.shape

# fig,ax=plt.subplots();ax.imshow(dataCube[50],cmap='gray');

x_center = 2800

y_center = 1900

half_square_size = 400

# rect = Rectangle((x_center - square_size/2, y_center - square_size/2), square_size, square_size,\
#                  linewidth=2, edgecolor='red', facecolor='none')  



# # Add the patch to the axes

# ax.add_patch(rect)

# fig.show()

lc_int = []

for i in range(np.shape(dataCube)[0]):
    lc_int.append(np.sum(dataCube[i,y_center-half_square_size:y_center+half_square_size,x_center-half_square_size:x_center+half_square_size])/(1600**2))

all_lc =[]      
for i in range(np.shape(dataCube)[0]):
    all_lc.append(np.sum(dataCube[i,:,:])/(4095**2))
         
      
dataCube=fits.open(savfold+'postdestretch_histomatch_dataCubeX_class_decay_second_set.fits')[0].data

gbref=dataCube[0,:,:] #choosing here the first image as a reference
timegb,yy,xx=dataCube.shape

for i in range(np.shape(dataCube)[0]):
    lc_int.append(np.sum(dataCube[i,y_center-half_square_size:y_center+half_square_size,x_center-half_square_size:x_center+half_square_size])/(1600**2))
           

        
# dataCube=fits.open(savfold+'postdestretch_histomatch_dataCubeX_class_decay_third_set.fits')[0].data

# gbref=dataCube[0,:,:] #choosing here the first image as a reference
# timegb,yy,xx=dataCube.shape

# for i in range(np.shape(dataCube)[0]):
#     lc_int.append(np.sum(dataCube[i,2098-600:2098+600,2645-600:2645+600])/(1200**2))
                
       
# dataCube=fits.open(savfold+'postdestretch_histomatch_dataCubeX_class_decay_fourth_set.fits')[0].data

# gbref=dataCube[0,:,:] #choosing here the first image as a reference
# timegb,yy,xx=dataCube.shape

# for i in range(np.shape(dataCube)[0]):
#     lc_int.append(np.sum(dataCube[i,2098-600:2098+600,2645-600:2645+600])/(1200**2))
                
# dataCube=fits.open(savfold+'postdestretch_histomatch_dataCubeX_class_decay_fifth_set.fits')[0].data

# gbref=dataCube[0,:,:] #choosing here the first image as a reference
# timegb,yy,xx=dataCube.shape

# for i in range(np.shape(dataCube)[0]):
#     lc_int.append(np.sum(dataCube[i,2098-600:2098+600,2645-600:2645+600])/(1200**2))
                             
              
       
              
              
              
              
              
              
               