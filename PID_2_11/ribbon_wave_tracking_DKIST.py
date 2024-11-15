#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 16:27:38 2024

@author: coletamburri
"""

#import dkistpkg_ct as DKISTanalysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.visualization import ImageNormalize, SqrtStretch
from matplotlib import animation
import matplotlib.patches as patches

import os

# import sunpy.coordinates
# import sunpy.map
# from sunpy.net import Fido
# from sunpy.net import attrs as a

from astropy.io import fits

path = '/Volumes/VBI_External/'
folder_vbi = 'AWYMX'

dir_list = os.listdir(path+folder_vbi)
dir_list.sort()
dir_list.pop(0)
dir_list.pop(0)

# img_file = fits.open(path+folder_vbi+'/'+dir_list[500])
# frame = img_file[1].data[0]
# fig,ax=plt.subplots()
# ax.pcolormesh(frame)
# #ax.invert_yaxis()
# ax.set_aspect('equal')
# rect = patches.Rectangle((1100, 2700), 650, 650, linewidth=1, edgecolor='r', facecolor='none')
# ax.add_patch(rect)


img_file = fits.open(path+folder_vbi+'/'+dir_list[0])
# frame = img_file[1].data[0][3200:3500,1600:2000]
frame = img_file[1].data[0]
fig, ax = plt.subplots(dpi=150)
cax = ax.pcolormesh(frame)
ax.invert_yaxis()
ax.set_aspect('equal')
props = dict(boxstyle='round', facecolor='white', alpha=1)
ax.text(20,20, img_file[1].header['DATE-BEG'][-15:-7],bbox=props)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])

def animate(i):
    ax.cla()
    img_file = fits.open(path+folder_vbi+'/'+dir_list[i])
    
    # frame = img_file[1].data[0][2700:3350,1100:1750]
    frame = img_file[1].data[0]
    # grad = np.gradient(frame)
    
    # arr1 = np.where(grad[0]>.42*(np.max(grad[0])))
    # arr2 = np.where(grad[1]<.62*(np.min(grad[1])))
    # arr3 = np.where(grad[1]>.45*(np.max(grad[1])))
    
    #magnitude = np.sqrt(grad[0]**2 + grad[1]**2)
    # maskmag = np.where(magnitude>4.5*np.median(magnitude))
    #ax.contour(frame,alpha=.5)
    ax.pcolormesh(frame,cmap='sdoaia1600')
    #ax.scatter(maskmag[1],maskmag[0],1,c='red')
    #ax.pcolormesh(magnitude,cmap='grey',alpha=1)
    
    #ax.set_xlim([0,650])
    #ax.set_ylim([0,650])
    ax.set_aspect('equal')
    
    # xs = list(arr1[0])+list(arr2[0])+list(arr3[0])
    # ys = list(arr1[1])+list(arr2[1])+list(arr3[1])
    ax.text(20,20, img_file[1].header['DATE-BEG'][-15:-7],bbox=props)
    #ax.scatter(ys,xs,1,alpha=1,c='red')
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    ax.invert_xaxis() # if VBI blue
    ax.invert_yaxis()


anim = animation.FuncAnimation(fig, animate,frames=200)
anim.save('bluecontinuum.gif')
plt.show()
