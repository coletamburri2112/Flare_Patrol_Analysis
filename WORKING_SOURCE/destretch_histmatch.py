#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 16:54:38 2025

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

pid='2_11'
vbifold = '/Volumes/VBI_External/pid_'+pid+'/'
# vbiexp ='BXWNO' # 19 August 2022
#vbiexp = 'BDJKM' # 11 August 2024 M-class
vbiexp = 'L1-HDGZO' # 11 August 2024 C-class with new ID
#vbiexp = 'AKDKX' # 11 August 2024 C-class
#vbiexp = 'AZLRR' # 8 August 2024 X-class blue continuum
# vbiexp = 'AWYMX' # 8 August 2024 M-class

# file = '/VBI_2022_08_19T20_42_07_333_00656282_I_BXWNO_L1.fits' # 19 August 2022
#file = '/VBI_2024_08_08T20_12_32_333333_00656282_I_AXXJL_L1.fits' # 11 August 2024 M-class
#file = '/VBI_2024_08_11T20_12_34_333333_00450400_I_AZLRR_L1.fits'
file = '/VBI_2024_08_11T20_12_34_333333_00656282_I_BDJKM_L1.fits'
savfold='/Users/coletamburri/Desktop/'+vbiexp+'/'
filt='Halpha'
xtraflag = 'C_class_impulsive_phase'

hdu_list = fits.open(vbifold+vbiexp+file)
image=hdu_list[1].data[0,:,:]

tileSizeInput = [128, 64, 48, 24]  # original from FW
#tileSizeInput = [64,48,24] # test from CAT

if os.path.isdir(savfold)=='False':
    os.mkdir(savfold)

# #display
# #plt.imshow(image,cmap="gray") # one could use the origin='lower' command to flip the image
# plt.imshow(image,cmap="gray",vmax=50000,vmin=5000) 
# plt.show()

#ds = dkist.load_dataset(vbifold+vbiexp+'/')

dir_list = os.listdir(vbifold+vbiexp)
dir_list.sort()

#ds[0,:,:].plot(cmap="gray",vmax=50000,vmin=5000)
#plt.show()

dataCube=[]
#loc_files=ds.files.filenames[200:450] # I think this should be the original destretching for C-class flare
loc_files = dir_list[150:450]
lowx=0
highx=-1
lowy=0
highy=-1

for i in range(len(loc_files)):
    dataCube.append(fits.open(vbifold+vbiexp+'/'+loc_files[i])[1].data[0,lowx:highx,lowy:highy])

dataCube=np.asarray(dataCube)

#save pre-destretch data cube of the first 100 (sub-)frames
fits.writeto(savfold+'pre-destretch_dataCube'+xtraflag+'.fits',dataCube,overwrite=True)

#restore datacube
dataCube=fits.open(savfold+'pre-destretch_dataCube'+xtraflag+'.fits')[0].data

#make a movie 
utilvbi.storeSequence(dataCube,savfold+'pre-destretch_'+xtraflag+'.mp4', dpi=300, write=True)

Video(savfold+'re-destretch_'+filt+'_'+xtraflag+'.mp4', embed=True, width=600, height=600)

#restore datacube
dataCube=fits.open(savfold+'pre-destretch_dataCube'+xtraflag+'.fits')[0].data

#destretch
dataCubeTracked = destr.destretchSeq(dataCube, tileSizeInput, rMean=3, globalTrack=None) 
#save as Tracked datacube 
fits.writeto(savfold+'postdestretch_dataCube'+xtraflag+'.fits',dataCubeTracked,overwrite=True)

#restore tracked datacube
dataCubeTracked=fits.open(savfold+'postdestretch_dataCube'+xtraflag+'.fits')[0].data

utilvbi.storeSequence(dataCubeTracked,savfold+'postdestretch_'+filt+'_'+xtraflag+'.mp4', dpi=300, write=True)

#histo match
#restore tracked datacube
dataCubeTracked=fits.open(savfold+'postdestretch_dataCube'+xtraflag+'.fits')[0].data

gbref=dataCubeTracked[0,:,:] #choosing here the first image as a reference
timegb,yy,xx=dataCubeTracked.shape
   
dataCubeTrackedhist=[]
for idata in range(timegb):
 image=dataCubeTracked[idata,:,:]
 matched = match_histograms(image, gbref)
 dataCubeTrackedhist.append(matched)
   
dataCubeTrackedhist=np.asarray(dataCubeTrackedhist)
    
#save result and movie  
fits.writeto(savfold+'postdestretch_histomatch_dataCube'+xtraflag+'.fits',dataCubeTrackedhist,overwrite=True)#save movie
dataCubeTracked=fits.open(savfold+'postdestretch_dataCube'+xtraflag+'.fits')[0].data

gbref=dataCubeTracked[0,:,:] #choosing here the first image as a reference
timegb,yy,xx=dataCubeTracked.shape
   
dataCubeTrackedhist=[]
for idata in range(timegb):
 image=dataCubeTracked[idata,:,:]
 matched = match_histograms(image, gbref)
 dataCubeTrackedhist.append(matched)
 
utilvbi.storeSequence(dataCubeTrackedhist,savfold+'postdestretch_histomatch_'+filt+'_'+xtraflag+'.mp4', dpi=300, write=True)

Video(savfold+'postdestretch_histomatch_'+filt+'_'+xtraflag+'.mp4', embed=True, width=600, height=600)



