#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 12:42:33 2025

@author: coletamburri
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

from sunpy.net import Fido, attrs as a
import dkist.net
from datetime import date, datetime, timedelta, time
import matplotlib.dates as mdates
import sunpy.visualization.colormaps as cm
import matplotlib

def perdelta(start, end, delta):
    curr = start
    while curr <= end:
        yield curr
        curr += delta
        
def normalize_3d_array(arr):
    min_val = np.nanmin(arr)
    max_val = np.nanmax(arr)

    return (arr - min_val) / (max_val - min_val)

#make times
sttime = datetime(2024,8,8,20,12,32,333333)
endtime = datetime(2024,8,8,21,5,7,0)


stack=[]
for result in perdelta(sttime , endtime, timedelta(seconds=2.666666)):
    stack.append(str(result))
    
timeshhmmss = []

# timesdt = []

# for i in stack:
#     timesdt.append(datetime.strptime(i, '%Y-%m-%dT%H:%M:%S.%f'))
    
# timesdthr = []

# for i in timesdt:
#     timesdthr.append(datetime.strftime(i, 
#                                  "%H:%M:%S"))

# times in correct format for plotting
for i in range(len(stack)):
    timeshhmmss.append(stack[i][-15:-7])
    
path = '/Volumes/VBI_External/pid_2_11/'

folder_vbi = 'AWYMX'
dir_list = os.listdir(path+folder_vbi)
dir_list.sort()
dir_list2 = []

#for i in range(len(dir_list)):
for i in range(1182):
    filename = dir_list[i]
    if filename[-5:] == '.fits' and '_I_' in filename:
        dir_list2.append(filename)
    
# load fits file
times = []
xarrs = []
yarrs = []

for i in range(50):
    hdul = fits.open(path+folder_vbi+'/'+dir_list2[i])
    times.append(hdul[1].header['DATE-BEG'])
    xloc = hdul[1].header['CRVAL1']
    yloc = hdul[1].header['CRVAL2']
    xpix = hdul[1].header['CRPIX1']
    ypix = hdul[1].header['CRPIX2']
    xdelt = hdul[1].header['CDELT1']
    ydelt = hdul[1].header['CDELT2']
    
    xarr = np.zeros(4095)
    yarr = np.zeros(4095)
    
    xarr[int(xpix)]=xloc
    yarr[int(ypix)]=yloc
    
    for j in range(int(xpix)+1,4095,1):
        xarr[j] = xarr[j-1] + xdelt
    for j in range(int(ypix)+1,4095,1):
        yarr[j] = yarr[j-1] + ydelt
        
    for j in range(int(xpix),-1,-1):
        xarr[j] = xarr[j+1] - xdelt
    for j in range(int(ypix),-1,-1):
        yarr[j] = yarr[j+1] - ydelt
        
    xarrs.append(xarr)
    yarrs.append(yarr)
    
    
    
    
#this file is indices 150 to 250
#data1 = np.load('/Users/coletamburri/Desktop/VBI_Destretching/AXXJL/AXXJLselections.npz')
#data = normalize_3d_array(data1['first50'])

#data = fits.open('/Volumes/VBI_External/postdestretch_dataCubeX_class_decay_full.fits')
data = fits.open('/Users/coletamburri/Desktop/VBI_Destretching/AWYMX/postdestretch_histomatch_dataCubeX_class_decay_blue_continuum.fits')

props = dict(edgecolor='black',facecolor='white', alpha=0.8,boxstyle='square,pad=0.4')

# fig,ax=plt.subplots(1,1,dpi=300,figsize=(10,5))
# for i in range(1):
#     #X,Y=np.meshgrid(xarrs[i],np.transpose(yarrs[i]))
#     ax.flatten()[i].pcolormesh(data[0].data[i*6,:,:],cmap='grey')
#     ax.flatten()[i].set_xticklabels([])
#     ax.flatten()[i].set_yticklabels([])
#     ax.flatten()[i].invert_yaxis()
#     # place a text box in upper left in axes coords
#     #ax.flatten()[i].text(375, -150, str(timeshhmmss[i*6])+' UT',fontsize=5,
#     #        verticalalignment='top', bbox=props)
#     ax.flatten()[i].set_aspect('equal')
#     ax.flatten()[i].tick_params(axis='both', which='minor', labelsize=5)
    

# fig.subplots_adjust(wspace=0, hspace=0)
# fig.show()

times=[]
xarrs=[]
yarrs=[]

i = 5

hdul = fits.open(path+folder_vbi+'/'+dir_list2[i])
#hdul = fits.open()
#times.append(hdul[1].header['DATE-BEG'])
xloc = hdul[1].header['CRVAL1']
yloc = hdul[1].header['CRVAL2']
xpix = hdul[1].header['CRPIX1']
ypix = hdul[1].header['CRPIX2']
xdelt = hdul[1].header['CDELT1']
ydelt = hdul[1].header['CDELT2']

xarr = np.zeros(4095)
yarr = np.zeros(4095)

xarr[int(xpix)]=xloc
yarr[int(ypix)]=yloc

for j in range(int(xpix)+1,4095,1):
    xarr[j] = xarr[j-1] + xdelt
for j in range(int(ypix)+1,4095,1):
    yarr[j] = yarr[j-1] + ydelt
    
for j in range(int(xpix),-1,-1):
    xarr[j] = xarr[j+1] - xdelt
for j in range(int(ypix),-1,-1):
    yarr[j] = yarr[j+1] - ydelt
    
xarrs.append(xarr)
yarrs.append(yarr)

X,Y = np.meshgrid(xarrs,yarrs[0][::-1])
    
    
print(stack[i])

fig,ax=plt.subplots(dpi=200,figsize=(20,20))
#ax.pcolormesh(X,Y,data[0].data[i],cmap=matplotlib.colormaps['sdoaia304'])
ax.pcolormesh(X,Y,data[0].data[i],cmap=matplotlib.colormaps['sdoaia1700'])
ax.set_aspect('equal')
#ax.invert_yaxis()
ax.invert_xaxis()
ax.tick_params(axis='both', which='minor', labelsize=2)

