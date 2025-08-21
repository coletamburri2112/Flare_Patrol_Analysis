#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 13:01:06 2025

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

import dkistpkg_ct as DKISTanalysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sunpy.visualization.colormaps as cm
import matplotlib


  

matplotlib.rcParams.update(matplotlib.rcParamsDefault)

# color scheme for plotting
muted = DKISTanalysis.color_muted2()

# path and file ID for ViSP data
path = '/Volumes/VBI_External/pid_2_11/'
folder1 = 'AXXJL' 
#folder2 = 'BDPMQ' # for QS calibration - data from 22UT on 19 August

# list of files in directory for DKIST/ViSP
dir_list2 = DKISTanalysis.pathdef(path,folder1) #flaretime
#dir_list3 = DKISTanalysis.pathdef(path,folder2) #qs


fried=[]
aolock = []
#for i in range(len(dir_list2)):
#    i_file_raster1 = fits.open(path+folder1+'/'+dir_list2[i])  
#    fried.append(i_file_raster1[1].header['ATMOS_R0'])
#    aolock.append(i_file_raster1[1].header['AO_LOCK'])
    
dC = fits.open('/Volumes/VBI_External/postdestretch_dataCubeX_class_decay_full.fits')[0].data # Xclass
#dC = fits.open('/Users/coletamburri/Desktop/DKIST_Code/VBI_Destretching/BDJKM/postdestretch_dataCube.fits')[0].data # Mclass
#dC = fits.open('/Users/coletamburri/Desktop/DKIST_Code/VBI_Destretching/AKDKX/postdestretch_dataCubeFlareImpulsivePhase.fits')[0].data
dCslice=dC[0:300,:,:] #Xclass
#dCslice = dC #Mclass or C class (100 and 250 frames, respectively

#friedarr = np.asarray(fried)
#ndices = np.argwhere(friedarr[0:350]>0.035)
indices = np.arange(len(dCslice))

def storeSequence(data, movieName, dpi=300, write=True, inds = indices):
    fig =plt.figure(dpi=300)
    max1=np.max(data[0])
    min1=np.min(data[0])
    #norm = (data[0]-min1)/(max1-min1)
    norm=data[0]/max1
    im = plt.imshow(norm,cmap=matplotlib.colormaps['afmhot'],vmin=0,vmax=.98)
    plt.tick_params(
      axis='x',          # changes apply to the x-axis
      which='both',      # both major and minor ticks are affected
      bottom=False,      # ticks along the bottom edge are off
      top=False,         # ticks along the top edge are off
      labelbottom=False) # labels along the bottom edge are off
    
    plt.tick_params(
      axis='y',          # changes apply to the x-axis
      which='both',      # both major and minor ticks are affected
      left=False,      # ticks along the bottom edge are off
      right=False,         # ticks along the top edge are off
      labelleft=False)
    fig.tight_layout()

    def animate(n):
        if n in inds:
            print(n)
            norm=data[n]/max1
            #norm=(data[n]-min1)/(max1-min1)
            im.set_data(norm)
            return im


    ani = animation.FuncAnimation(fig, animate, frames=data.shape[0], interval=100)
    
    if write:
        writer = animation.writers['ffmpeg'](fps=40)
        ani.save(movieName, writer=writer, dpi=dpi)
    plt.close(fig)  



#dataCubeTracked=fits.open('/Volumes/VBI_External/postdestretch_dataCubeX_class_decay_full.fits')[0].data

#dataCubeTracked=fits.open('/Volumes/VBI_External/postdestretch_dataCubeX_class_decay_full.fits')[0].data

storeSequence(dCslice,'/Users/coletamburri/Desktop/movie_Xclass_faster.mp4', dpi=300, write=True,inds=indices)
