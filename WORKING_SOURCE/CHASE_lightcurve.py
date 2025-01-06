#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 07:29:29 2024

@author: coletamburri
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from scipy.integrate import cumulative_trapezoid


filesdir = '/Users/coletamburri/Desktop/Halpha_CHASE_data/'

files = os.listdir(filesdir)

files.sort()

def veltrans(x):
    return ((((x+lamb0)/lamb0)-1)*c)/mu

def veltrans2(x):
    return ((((x+lamb0)/lamb0)-1)*c)/mu2

def wltrans(x):
    return ((((x/c)+1)*lamb0)-lamb0)

def normalize(data):
    l=0
    if np.sum(data) == 0.0:
        l+=1
        normarr = dataf
    else:
        normarr=(data-np.nanmin(data))/(np.nanmax(data)-np.nanmin(data)) 
    return normarr

centwl = 656.46
lamb0=centwl
c=2.99e5
mu=1

halphadat = []
for i in range(len(files)-1):
    file = fits.open(filesdir+files[i+1])
    data = file[1].data[:,740:855,1200:1290]
    halphadat.append(data/0.000525)
    file.close()
    
centwl = file[1].header['WAVE_LEN']/10
    
wavelengths = np.arange(file[1].header['CRVAL3'],file[1].header['CRVAL3']+\
                        file[1].header['CDELT3']*118,file[1].header['CDELT3'])
    
meanhalpha = np.mean(halphadat[0],0)

# fig,ax=plt.subplots(8,3,dpi=200)

# for i in range(np.shape(halphadat)[0]):
#     meanhalpha = np.mean(halphadat[i],0)

#     ax.flatten()[i].pcolormesh(meanhalpha[740:855,1200:1290])
    
# halphadatarr = np.asarray(halphadat)

intensity = []

for i in range(23):
    #slice = halphadatarr[i]
    bit = halphadat[i]
    flux = 0
    #flux = cumtrapz(wavelengths,slice[:,710:870,1150:1330])
    for j in range(np.shape(halphadat)[2]):
        for k in range(np.shape(halphadat)[3]):
    #for j in np.arange(1000,1160,1):
    #    for k in np.arange(1000,1180,1):
            flux+= cumulative_trapezoid(bit[:,j,k],wavelengths[:-1])[-1]
    intensity.append(flux)
    
times=[]
for i in range(len(files)-1):
    times.append(files[i+1][-19:-17]+':'+files[i+1][-17:-15]+':'+files[i+1][-15:-13])

fig,ax=plt.subplots(1,2,figsize=(10,5),dpi=150)
ax.flatten()[0].plot(times,intensity,linestyle='--', marker='o')
ax.flatten()[1].plot(times,normalize(intensity),linestyle='--', marker='o')
ax.flatten()[1].set_xticks(times[0::4],times[0::4],rotation=30)
ax.flatten()[0].set_xticks(times[0::4],times[0::4],rotation=30)
ax.flatten()[1].set_ylabel('Normalized Intensity')
ax.flatten()[1].set_xlabel('Time [UT]')
ax.flatten()[0].set_xlabel('Time [UT]')
