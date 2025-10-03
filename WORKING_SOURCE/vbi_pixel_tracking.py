#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 09:00:58 2025

@author: coletamburri
"""

import dkistpkg_ct as DKISTanalysis
import numpy as np
import matplotlib.pyplot as plt
import tol_colors
import os
from astropy.io import fits
from datetime import time
import matplotlib.patches as patches



muted = DKISTanalysis.color_muted2()

l=0
CUT=2.5 # cutoff for mask-making
binning = 0 # to bin or not to bin?
diff = 0 # to diff or not to diff
if binning == 1:
    bin_x = 7 # this bins the VBI pixels to largest ViSP spatial scale (in scan dir), plus 1 for "safety"
    bin_y = 7
else:
    bin_x = 1 #initialize this parameter - will change if binning selected, below
    bin_y = 1 #initialize this parameter - will change if binning selected, below

# limits for flare region (defined by raw image)
xlow = 500 # 500 for include r2, 1700 for only r1
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

visp_file = '/Users/coletamburri/Desktop/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/ViSP_coalign_result_11Aug_Cclass'

vispload = np.load(visp_file)

vispX = vispload['arr_0']
vispY = vispload['arr_1']
vispavg = vispload['arr_2']




timesvbi=[]
for i in range(len(dir_list2_vbi)):
    timesvbi.append(time(int(dir_list2_vbi[i][15:17]),
                      int(dir_list2_vbi[i][18:20]),
                      int(dir_list2_vbi[i][21:23])))
    
timesvbi = timesvbi[200:450]

#another time array
string = timesvbi[0].strftime("%H:%M:%S")
strtime =[]
for i in range(len(timesvbi)):
    strtime.append(timesvbi[i].strftime("%H:%M:%S"))

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



arr = vbi_DSimgs

def rebin_image(arr, new_shape):
    """
    Rebins a 2D array (image) to a new shape by averaging pixel values.

    Parameters:
    arr (numpy.ndarray): The original 2D image array.
    new_shape (tuple): The desired new shape (rows, columns).
                       Each dimension of new_shape must be a factor of
                       the corresponding dimension in arr.shape.

    Returns:
    numpy.ndarray: The binned 2D array.
    """
    if not (arr.shape[0] % new_shape[0] == 0 and arr.shape[1] % new_shape[1] == 0):
        raise ValueError("New shape dimensions must be factors of original shape dimensions.")

    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)



if binning == 1:
    new_height = vbi_DSimgs.shape[1] // bin_y
    new_width = vbi_DSimgs.shape[2] // bin_x
    
    binned = np.zeros((250,new_height,new_width))
    
    for i in range(250):
        binned[i,:,:] = rebin_image(vbi_DSimgs[i,:,:], (new_height,new_width))
    
    arr = binned

# create cumulative mask
mask = np.zeros(np.shape(arr))
timing = np.zeros(np.shape(arr)[1:3])

for i in range(100):
    l+=1
    if i>0:
        mask[i,:,:]=mask[i-1,:,:]
    maskvals = np.nonzero((arr[i,:,:]>CUT*np.median(arr[i,:,:])))
    for j in range(np.shape(maskvals)[1]):
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

            
# plot timing series (evolution of ribbon mask)
fig,ax=plt.subplots(dpi=200)
pcm=ax.pcolormesh(timing,cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'))
ax.set_xlim([1600/bin_x,2300/bin_x])
ax.set_ylim([2700/bin_y,900/bin_y])
ax.set_aspect('equal')
cbar = fig.colorbar(pcm, ax=ax,ticks=[1,20,40,60,80,100])
cbar.ax.set_yticklabels([timesvbi[1],timesvbi[20],timesvbi[40],timesvbi[60],timesvbi[80],timesvbi[100]])
ax.set_xticks([])
ax.set_yticks([])
fig.show()

#indices for light curves are in last cumulative mask
i_vals = []
j_vals = []

final_mask = mask[99,:,:]

for i in range(np.shape(final_mask)[0]):
    for j in range(np.shape(final_mask)[1]):
        if final_mask[i,j] > 0.0 and j>xlow/bin_x and j<xhigh/bin_x and i>ylow/bin_x and i<yhigh/bin_x:
            i_vals.append(i)
            j_vals.append(j)
            
# make light curve
lc=[]

for i in range(250):
    lc.append(np.sum(arr[i,int(ylow/bin_y):int(yhigh/bin_y),int(xlow/bin_x):int(xhigh/bin_x)]))
    
#compare light curves
numcolor_timing = int(np.nanmax(timing))
maps = tol_colors.tol_cmap(colormap='rainbow_PuRd',lut=numcolor_timing)
cmap_choice2 = maps(np.linspace(0,1,numcolor_timing))

fig,ax=plt.subplots(figsize=(20,10))
# for i in np.arange(0,len(i_vals)):
for i in np.arange(1,len(i_vals)):
    if np.nanmax(arr[0:119,i_vals[i],j_vals[i]])>50000:
        ax.plot(arr[0:119,i_vals[i],j_vals[i]],c=cmap_choice2[int(timing[i_vals[i],j_vals[i]])-1],linewidth=.1*bin_x)

ax.plot(strtime[0:100],arr[0:100,i_vals[0],j_vals[0]],c=cmap_choice2[int(timing[i_vals[i],j_vals[i]])-1],linewidth=.1*bin_x)
ax2=ax.twinx()
ax2.plot(strtime[0:100],lc[0:100],c='crimson',linewidth=5,marker='.')
ax.set_xticks(strtime[0:101:20])
ax.set_ylim([-50000,300000])


# difference imaging

if diff == 1:
    diffarr = np.zeros(np.shape(arr))
    
    m=0
    folder = '/Users/coletamburri/Desktop/diffimg_pre/'
    if ~os.path.isdir(folder):
        os.mkdir(folder)
    
    for i in range(np.shape(diffarr)[0]-20):
        diffarr[m,:,:] = np.subtract(arr[i,:,:],arr[i+20,:,:])
        
        
        fig,ax=plt.subplots(dpi=200);
        ax.imshow(diffarr[m,:,:],cmap='grey')
        fig.savefig(folder+str(i)+'.png')
        
        m+=1
        
    
# instantaneous mask movie
for i in range(1,100):
    inds = np.where(timing==int(i))
    arrsamp = np.zeros(np.shape(timing))
    for j in range(np.shape(inds)[1]):
        arrsamp[inds[0][j],inds[1][j]]=timing[inds[0][j],inds[1][j]]
    arrsamp[arrsamp == 0] = np.nan
    fig,ax=plt.subplots(dpi=200)
    pcm = ax.pcolormesh(arrsamp,cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'), vmin=0, vmax=100)
    ax.set_xlim([1600/bin_x,2300/bin_x])
    ax.set_ylim([2700/bin_y,900/bin_y])
    ax.set_aspect('equal')
    cbar = fig.colorbar(pcm, ax=ax,ticks=[1,20,40,60,80,100])
    cbar.ax.set_yticklabels([timesvbi[1],timesvbi[20],timesvbi[40],timesvbi[60],timesvbi[80],timesvbi[100]])
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig('/Users/coletamburri/Desktop/maskseq/'+str(i)+'.png')
    
for i in range(1,100):
    inds = np.where(timing==int(i))
    arrsamp = np.zeros(np.shape(timing))
    for j in range(np.shape(inds)[1]):
        arrsamp[inds[0][j],inds[1][j]]=timing[inds[0][j],inds[1][j]]
    arrsamp[arrsamp == 0] = np.nan
    fig,ax=plt.subplots(dpi=200)
    ax.pcolormesh(arrsamp,cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'), vmin=0, vmax=100)
    ax.set_ylim([1890/bin_x,1730/bin_x])
    ax.set_xlim([1750/bin_y,1900/bin_y])
    ax.set_aspect('equal')
    cbar = fig.colorbar(pcm, ax=ax,ticks=[1,20,40,60,80,100])
    cbar.ax.set_yticklabels([timesvbi[1],timesvbi[20],timesvbi[40],timesvbi[60],timesvbi[80],timesvbi[100]])
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig('/Users/coletamburri/Desktop/maskseq_small/'+str(i)+'.png')
    
# cumulative mask movie
cumul_mask = np.zeros([100,np.shape(timing)[0],np.shape(timing)[1]])

for i in range(1,99):
    inds = np.where(timing==int(i))
    arrsamp = np.zeros(np.shape(timing))
    for j in range(np.shape(inds)[1]):
        arrsamp[inds[0][j],inds[1][j]]=timing[inds[0][j],inds[1][j]]
    cumul_mask[i,:,:]=cumul_mask[i-1,:,:]+arrsamp
    
cumul_mask[cumul_mask == 0] = np.nan

inst_mask = np.zeros([100,np.shape(timing)[0],np.shape(timing)[1]])

for i in range(1,99):
    inds = np.where(timing==int(i))
    arrsamp = np.zeros(np.shape(timing))
    for j in range(np.shape(inds)[1]):
        arrsamp[inds[0][j],inds[1][j]]=1.0
    inst_mask[i,:,:]=arrsamp

inst_mask[inst_mask == 0] = np.nan

for i in range(1,99):
    fig,ax=plt.subplots(dpi=200)
    pcm=ax.pcolormesh(cumul_mask[i,:,:],cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'), vmin=0, vmax=100)
    ax.set_xlim([1600/bin_x,2300/bin_x])
    ax.set_ylim([2700/bin_y,900/bin_y])
    ax.set_aspect('equal')
    cbar = fig.colorbar(pcm, ax=ax,ticks=[1,20,40,60,80,100])
    cbar.ax.set_yticklabels([timesvbi[1],timesvbi[20],timesvbi[40],timesvbi[60],timesvbi[80],timesvbi[100]])
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig('/Users/coletamburri/Desktop/cumul_maskseq/'+str(i)+'.png')

# light curve of region

lcsmall=[]

for i in range(100):
    lcsmall.append(np.nansum(arr[i,1770:1840,1790:1840])) #limits for small thing
    
fig,ax=plt.subplots();ax.plot(lcsmall)




for i in range(0,100):
    fig,[ax,ax1,ax2]=plt.subplots(1,3,dpi=200)
    ax.imshow(arr[i,1770:1840,1790:1840],cmap='hot')
    ax1.imshow(inst_mask[i,1770:1840,1790:1840],cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'), vmin=0, vmax=40)
    ax2.imshow(cumul_mask[i,1770:1840,1790:1840],cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'), vmin=0, vmax=40)
    ins = ax.inset_axes([0.5,0.7,0.4,0.2])
    ins.plot(lcsmall,c='black')
    ins.axvline(i,c='red')
    ins.set_xticks([])
    ins.set_yticks([])
    fig.savefig('/Users/coletamburri/Desktop/kernzoom/'+str(i)+'.png')

for i in range(0,40):
    fig,[(ax,ax1,ax2),(ax3,ax4,ax5)]=plt.subplots(2,3,dpi=200)
    ax.imshow(arr[i,1770:1840,1790:1840],cmap='hot')
    ax1.imshow(inst_mask[i,1770:1840,1790:1840],cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'), vmin=0, vmax=40)
    ax2.imshow(cumul_mask[i,1770:1840,1790:1840],cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'), vmin=0, vmax=40)
    ax3.imshow(arr[i,900:2700,1600:2300],cmap='hot')
    ax4.imshow(inst_mask[i,900:2700,1600:2300],cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'), vmin=0, vmax=40)
    ax5.imshow(cumul_mask[i,900:2700,1600:2300],cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'), vmin=0, vmax=40)
    rect = patches.Rectangle((1790-1600, 1770-900), 50, 70, linewidth=1, edgecolor='k', facecolor='none')

    # Add the patch to the Axes
    ax3.add_patch(rect)
    rect = patches.Rectangle((1790-1600, 1770-900), 50, 70, linewidth=1, edgecolor='k', facecolor='none')

    ax4.add_patch(rect)
    rect = patches.Rectangle((1790-1600, 1770-900), 50,70, linewidth=1, edgecolor='k', facecolor='none')

    ax5.add_patch(rect)

    ins = ax.inset_axes([0.5,0.7,0.4,0.2])
    ins.plot(lcsmall,c='black')
    ins.axvline(i,c='red')
    ins.set_xticks([])
    ins.set_yticks([])
    fig.savefig('/Users/coletamburri/Desktop/fullframe/'+str(i)+'.png',format='png',dpi=200)
    




















