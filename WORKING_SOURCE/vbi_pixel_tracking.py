#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 09:00:58 2025

@author: coletamburri
"""

import dkistpkg_ct as DKISTanalysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.visualization import ImageNormalize, SqrtStretch

import sunpy.coordinates
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
import tol_colors

from astropy.io import fits
from datetime import datetime
from datetime import time
import matplotlib.dates as mdates

muted = DKISTanalysis.color_muted2()

l=0
CUT=3.8 # cutoff for mask-making

# limits for flare region
xlow = 1700
xhigh = 2350
ylow = 1000
yhigh = 2800

#define start and end times for series
starttime = 200
endtime = 450

# use destretched (0) or pre-destretched (1)?
stretched = 0

#VBI directory
path_vbi = '/Volumes/VBI_External/pid_2_11/'
folder1_vbi = 'AKDKX'
#folder2_vbi = 'BYMOL'
dir_list2_vbi = DKISTanalysis.pathdef(path_vbi,folder1_vbi)

#destretched dataset
filename = ['/Users/coletamburri/Desktop/DKIST_Code/VBI_Destretching/AKDKX/'+
                'postdestretch_dataCubeFlareImpulsivePhase.fits'][0]


times=[]
for i in range(len(dir_list2_vbi)):
    times.append(time(int(dir_list2_vbi[i][15:17]),
                      int(dir_list2_vbi[i][18:20]),
                      int(dir_list2_vbi[i][21:23])))
    
times = times[200:450]

#another time array
string = times[0].strftime("%H:%M:%S")
strtime =[]
for i in range(len(times)):
    strtime.append(times[i].strftime("%H:%M:%S"))

#optional - processing just from the files on HD (raw, not destretched)
if stretched == 1:
    vbi_X, vbi_Y, hdul1_vbi, dat0_vbi = DKISTanalysis.vbi_process(path_vbi,
                                                                  folder1_vbi)
    #just pixels
    vbix0 = np.arange(4096)
    vbiy0 = np.arange(4096)
    vbiX0,vbiY0= np.meshgrid(vbix0,vbiy0)
else:
    vbi_DS = fits.open(filename)

vbi_DSimgs = vbi_DS[0].data

fig,ax=plt.subplots()
ax.imshow(vbi_DS[0].data[40,:,:])
fig.show()

# create cumulative mask
mask = np.zeros(np.shape(vbi_DSimgs))
timing = np.zeros(np.shape(vbi_DSimgs)[1:3])

for i in range(100):
    l+=1
    if i>0:
        mask[i,:,:]=mask[i-1,:,:]
    maskvals = np.nonzero((vbi_DSimgs[i,:,:]>CUT*np.median(vbi_DSimgs[i,:,:])))
    for j in range(np.shape(maskvals)[1]):
        #mask[i,maskvals[0][j],maskvals[1][j]] += vbi_DSimgs[i,maskvals[0][j],maskvals[1][j]]
        mask[i,maskvals[0][j],maskvals[1][j]] += 1
        #logic for timing array
        if timing[maskvals[0][j],maskvals[1][j]]==0:
            timing[maskvals[0][j],maskvals[1][j]]=l

#create timing series
for i in range(np.shape(timing)[0]):
    for j in range(np.shape(timing)[1]):
        if timing[i,j]==0:
            timing[i,j]='NaN'
            

# plot final cumulative mask
fig,ax=plt.subplots(dpi=200)
ax.imshow(mask[99,:,:],cmap='jet')
ax.set_xlim([1700,2200])
ax.set_ylim([2600,1200])
fig.show()
            
# plot timing series (evolution of ribbon mask)
fig,ax=plt.subplots(dpi=200)
pcm=ax.pcolormesh(timing,cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'))
ax.set_xlim([1600,2300])
ax.set_ylim([2700,900])
ax.set_aspect('equal')
cbar = fig.colorbar(pcm, ax=ax,ticks=[1,20,40,60,80,100])
cbar.ax.set_yticklabels([times[1],times[20],times[40],times[60],times[80],times[100]])
ax.set_xticks([])
ax.set_yticks([])
fig.show()

#indices for light curves are in last cumulative mask
i_vals = []
j_vals = []

final_mask = mask[99,:,:]

for i in range(np.shape(final_mask)[0]):
    for j in range(np.shape(final_mask)[1]):
        if final_mask[i,j] > 0.0 and j>xlow and j<xhigh and i>ylow and i<yhigh:
            i_vals.append(i)
            j_vals.append(j)
            
# make light curve
lc=[]

for i in range(250):
    lc.append(np.sum(vbi_DSimgs[i,ylow:yhigh,xlow:xhigh]))
    
#compare light curves
numcolor_timing = int(np.nanmax(timing))
maps = tol_colors.tol_cmap(colormap='rainbow_PuRd',lut=numcolor_timing)
cmap_choice2 = maps(np.linspace(0,1,numcolor_timing))

fig,ax=plt.subplots(figsize=(20,10))
# for i in np.arange(0,len(i_vals)):
for i in np.arange(1,len(i_vals)):
    if np.nanmax(vbi_DSimgs[0:119,i_vals[i],j_vals[i]])>50000:
        ax.plot(vbi_DSimgs[0:119,i_vals[i],j_vals[i]],c=cmap_choice2[int(timing[i_vals[i],j_vals[i]])-1],linewidth=0.1)

ax.plot(strtime[0:100],vbi_DSimgs[0:100,i_vals[0],j_vals[0]],c=cmap_choice2[int(timing[i_vals[i],j_vals[i]])-1],linewidth=0.1)
ax2=ax.twinx()
ax2.plot(strtime[0:100],lc[0:100],c='crimson',linewidth=5,marker='.')
ax.set_xticks(strtime[0:101:20])
ax.set_ylim([-50000,300000])



    
























