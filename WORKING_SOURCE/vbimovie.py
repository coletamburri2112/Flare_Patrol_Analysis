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

from datetime import datetime, timedelta

  

matplotlib.rcParams.update(matplotlib.rcParamsDefault)

# color scheme for plotting
muted = DKISTanalysis.color_muted2()

# path and file ID for ViSP data
path = '/Users/coletamburri/Desktop/'
folder1 = 'DXHIEL' 
#folder2 = 'BDPMQ' # for QS calibration - data from 22UT on 19 August

# list of files in directory for DKIST/ViSP
#dir_list2 = DKISTanalysis.pathdef(path,folder1) #flaretime
#dir_list3 = DKISTanalysis.pathdef(path,folder2) #qs


fried=[]
aolock = []
#for i in range(len(dir_list2)):
#    i_file_raster1 = fits.open(path+folder1+'/'+dir_list2[i])  
#    fried.append(i_file_raster1[1].header['ATMOS_R0'])
#    aolock.append(i_file_raster1[1].header['AO_LOCK'])
    
#dC = fits.open('/Volumes/VBI_External/postdestretch_dataCubeX_class_decay_full.fits')[0].data # Xclass
dC = fits.open('/Users/coletamburri/Desktop/DKIST_Code/VBI_Destretching/MBVIDS/postdestretch_dataCube_Halpha_C_class_impulsive_phase_Halpha_177_347.fits')
#dC = fits.open('/Users/coletamburri/Desktop/DXHIEL/postdestretch_dataCube_blue_cont_C_class_impulsive_phase.fits')[0].data # Xclass


#dC = fits.open('/Users/coletamburri/Desktop/DKIST_Code/VBI_Destretching/BDJKM/postdestretch_dataCube.fits')[0].data # Mclass
#dC = fits.open('/Users/coletamburri/Desktop/DKIST_Code/VBI_Destretching/AKDKX/postdestretch_dataCubeFlareImpulsivePhase.fits')[0].data
#dCslice=dC[0:300,:,:] #Xclass
dCslice = dC #Mclass or C class (100 and 250 frames, respectively

#friedarr = np.asarray(fried)
#ndices = np.argwhere(friedarr[0:350]>0.035)
indices = np.arange(len(dCslice))

# define times
t3 = np.arange(datetime(2024,8,11,22,31,26),
              datetime(2024,8,11,22,38,57), 
              timedelta(seconds=2.666)).astype(datetime)



#limits for entire flare region
ylow = 1600
yhigh = 2300
xlow = 900
xhigh = 2700

#limits for small part of region r1A ("bead d")
ylow1 = 1800
yhigh1 = 2300
xlow1 = 2000
xhigh1 = 2500

#limits for three small beads
ylow2 = 1900
yhigh2 = 2100
xlow2 = 2270
xhigh2 = 2470

#limits for regionr1B
ylow3 = 1600
yhigh3 = 2300
xlow3 = 900
xhigh3 = 1750

#limits for central outflow region
ylow4 = 1600
yhigh4 = 2300
xlow4 = 1600
xhigh4 = 2000

#limits for r2
ylow5= 200
yhigh5 = 1300
xlow5 = 1200
xhigh5 = 2500

loadvbilc = np.load('/Users/coletamburri/Desktop/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/vbi_lc.npz',allow_pickle='True')
#timesvbi=loadvbilc['times']
lcvbi=loadvbilc['lc']


def storeJPG(data,outfolder,times,lc,end=137,dpi=300,lowx=0,highx=-1,lowy=0,highy=-1,leftsubplot=0.05):
    for i in range(end):
        image = data[0].data[i,lowx:highx,lowy:highy]
        
        imagenorm = (image-image.min())/(image.max()-image.min())
        
        fig,ax=plt.subplots(dpi=dpi,figsize=(10,10))
        ax.imshow(imagenorm,cmap='sdoaia304')
        ax.set_title(str(times[i])+' UT',fontsize=12)
        
        ax_inset = ax.inset_axes([leftsubplot, 0.8, 0.4, 0.15]) 
        ax_inset.plot(lcvbi[:end],color='black')
        ax_inset.set_xticks([])
        ax_inset.set_yticks([])
        ax_inset.patch.set_alpha(0.7)
        ax_inset.axvline(i,color='red',linestyle='dashed')
        

        fig.savefig(outfolder+str(i)+'.png')
        
outfolder = '/Users/coletamburri/Desktop/fullframe/'
os.mkdir(outfolder)
storeJPG(dC,outfolder,t3,lcvbi)        

outfolder = '/Users/coletamburri/Desktop/regionr1B/'
os.mkdir(outfolder)
storeJPG(dC,outfolder,t3,lcvbi,lowx=xlow,highx=xhigh,lowy=ylow,highy=yhigh)

outfolder = '/Users/coletamburri/Desktop/regionr1A/'
os.mkdir(outfolder)
storeJPG(dC,outfolder,t3,lcvbi,lowx=xlow1,highx=xhigh1,lowy=ylow1,highy=yhigh1)

outfolder = '/Users/coletamburri/Desktop/threebeadsr1A/'
os.mkdir(outfolder)
storeJPG(dC,outfolder,t3,lcvbi,lowx=xlow2,highx=xhigh2,lowy=ylow2,highy=yhigh2)

outfolder = '/Users/coletamburri/Desktop/regionr1B/'
os.mkdir(outfolder)
storeJPG(dC,outfolder,t3,lcvbi,lowx=xlow3,highx=xhigh3,lowy=ylow3,highy=yhigh3)

outfolder = '/Users/coletamburri/Desktop/outflow/'
os.mkdir(outfolder)
storeJPG(dC,outfolder,t3,lcvbi,lowx=xlow4,highx=xhigh4,lowy=ylow4,highy=yhigh4)

outfolder = '/Users/coletamburri/Desktop/r2/'
os.mkdir(outfolder)
storeJPG(dC,outfolder,t3,lcvbi,lowx=xlow5,highx=xhigh5,lowy=ylow5,highy=yhigh5)





                  

def storeSequence(data, movieName, dpi=300, write=True, inds = indices):
    fig =plt.figure(dpi=300)
    max1=np.max(data[0])
    min1=np.min(data[0])
    norm = (data[0]-min1)/(max1-min1)
    #norm=data[0]/max1
    #im = plt.imshow(norm,cmap=matplotlib.colormaps['afmhot'],vmin=0,vmax=.98)
    im = plt.imshow(norm)
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

storeSequence(dCslice,'/Users/coletamburri/Desktop/movie_Cclass_bluecont_faster.mp4', dpi=300, write=True,inds=indices)
