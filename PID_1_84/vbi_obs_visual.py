#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 16:01:32 2023

@author: coletamburri
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib

plt.rcParams['text.usetex']=True
plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif'] = ['Tahoma']
plt.rcParams['axes.labelsize'] = 25
plt.rcParams['lines.linewidth'] = 2
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20) 

def color_muted2():
    muted =['#332288', '#88CCEE', '#44AA99','#117733','#999933','#DDCC77', '#CC6677','#882255','#AA4499','#DDDDDD']
#  0= indigo
# 1 = cyan
# 2 = teal
# 3 = green
# 4 = olive
# 5= sand
# 6 = rose
# 7 = wine
# 8 = purple
# 9=grey
    return muted

muted = color_muted2()

path = '/Volumes/VBI_External/pid_1_84/'
folders = os.listdir(path)


#for i in range(len(folders)):
#folder1 = 'AKDKX' #pid_2_11
folder1 = 'BPYLQ' #pid_1_84
print(folder1)
dir_list = os.listdir(path+folder1)

# if i>1:
#     del dir_list2
#     del dir_list
#     del hdul1
#     del fullarr
#     del xarrs
#     del yarrs
#     del timestamps

dir_list2 = []

for i in range(len(dir_list)):
    filename = dir_list[i]
    if filename[-5:] == '.fits' and '_I_' in filename:
        dir_list2.append(filename)

dir_list2.sort()
hdul1 = fits.open(path+folder1+'/'+dir_list2[0])




# for pid_2_11
if len(dir_list2)<100:
    length = len(dir_list2)
elif len(dir_list2)>99:
    length = 100
length=50

# for pid_1_84
length=100

fullarr = np.zeros([length,len(hdul1[1].data[0,0,:]),len(hdul1[1].data[0,0,:])])
xarrs = np.zeros([length,len(hdul1[1].data[0,0,:])])
yarrs = np.zeros([length,len(hdul1[1].data[0,0,:])])
timestamps = []

def normalize_3d_array(arr):
    min_val = np.nanmin(arr)
    max_val = np.nanmax(arr)

    return (arr - min_val) / (max_val - min_val)






for i in range(length):
    #hdul = fits.open(path+folder1+'/'+dir_list2[i+300]) #pid_2_11, 11Aug flare
    hdul = fits.open(path+folder1+'/'+dir_list2[i+200]) #pid_1_84, 11Aug flare

    fullarr[i,:,:] = hdul[1].data[0,:,:]

    #latitude (labeled as lon, but header is wrong)
    xcent = hdul[1].header['CRVAL2']
    xnum = hdul[1].header['NAXIS2']
    xdelt = hdul[1].header['CDELT2']

    #longitude (labeled as lat, but header is wrong)
    ycent = hdul[1].header['CRVAL1']
    ynum = hdul[1].header['NAXIS1']
    ydelt = hdul[1].header['CDELT1']

    timestamps.append(hdul[1].header['DATE-BEG'])

    xarr = np.linspace(((xcent-xdelt/2)-((xnum-1)/2)*xdelt),((xcent-xdelt/2)+((xnum-1)/2)*xdelt),xnum)
    yarr = np.linspace(((ycent-ydelt/2)-((ynum-1)/2)*ydelt),((ycent-ydelt/2)+((ynum-1)/2)*ydelt),ynum)

    xarrs[i,:] = xarr
    yarrs[i,:] = yarr
    
normalized_arr = normalize_3d_array(fullarr)
    
# fig,ax = plt.subplots(4,4,figsize=(10,70))

# for i in range(length):
#     j = i+1
#     ax1 = plt.subplot(4, 4, j)
#     plt.axis('on')
#     arr = fullarr[i,:,:].T
#     #tr = Affine2D().rotate_deg(45.)

#     ax1.imshow(arr,cmap='gray',extent=[yarr[0],yarr[-1],xarr[0],xarr[-1]],origin='upper')
#     ax1.set_xticklabels([])
#     ax1.set_yticklabels([])
#     ax1.set_aspect('equal')
#     ax1.invert_xaxis()
#     props = dict(boxstyle='square', facecolor='white', alpha=0.9)
#     #ax1.text(687,-506,timestamps[i].split()[0][-12:-4]+' UT',bbox=props,fontsize=12)
#     plt.subplots_adjust(wspace=0, hspace=0)

# ax1 = plt.subplot(25,4,1)
# ax1.set_xticks([yarr[100],yarr[2000],yarr[3900]],labels=[str(round(yarr[-100],1))+'"',str(round(yarr[-2000],1))+'"',str(round(yarr[-3900],1))+'"'],fontsize=12)
# ax1.set_yticks([xarr[100],xarr[2000],xarr[3900]],labels=[str(round(xarr[100],1))+'"',str(round(xarr[2000],1))+'"',str(round(xarr[3900],1))+'"'],fontsize=12)
# ax1.set_ylabel('HPC LAT',fontsize=15)
# ax1.set_xlabel('HPC LON',fontsize=15)
# ax1.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
# ax1.xaxis.set_label_position('top') 

# plt.savefig('/Users/coletamburri/Desktop/'+folder1+'.png')