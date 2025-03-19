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

from sunpy.net import Fido, attrs as a
import dkist.net
from datetime import date, datetime, timedelta, time
import matplotlib.dates as mdates

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
#sttime = datetime(2024,8,8,20,12,32,333333)
3endtime = datetime(2024,8,8,21,5,7,0)


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
    
# load fits file

#this file is indices 150 to 250
data1 = np.load('/Users/coletamburri/Desktop/VBI_Destretching/AXXJL/AXXJLselections.npz')
data = normalize_3d_array(data1['brightening'])

props = dict(edgecolor='black',facecolor='white', alpha=0.8,boxstyle='square,pad=0.4')

fig,[(ax1,ax2),(ax3,ax4)]=plt.subplots(2,2,dpi=300,figsize=(3,3))
ax1.imshow(data[0,:,:],cmap='grey',)
ax1.set_xticklabels([])
ax1.set_yticklabels([])

# place a text box in upper left in axes coords
ax1.text(0.65, 0.9, str(timeshhmmss[150+0])+' UT', transform=ax1.transAxes, fontsize=4,
        verticalalignment='top', bbox=props)
ax1.set_aspect('equal')

ax2.imshow(data[24,:,:],cmap='grey')
ax2.text(0.65, 0.9, str(timeshhmmss[150+24])+' UT', transform=ax2.transAxes, fontsize=4,
        verticalalignment='top', bbox=props)
ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.set_aspect('equal')

ax3.imshow(data[49,:,:],cmap='grey')
ax3.text(0.65, 0.9, str(timeshhmmss[150+49])+' UT', transform=ax3.transAxes, fontsize=4,
        verticalalignment='top', bbox=props)
ax3.set_xticklabels([])
ax3.set_yticklabels([])
ax3.set_aspect('equal')

ax4.imshow(data[74,:,:],cmap='grey')
ax4.text(0.65, 0.9, str(timeshhmmss[150+74])+' UT', transform=ax4.transAxes, fontsize=4,
        verticalalignment='top', bbox=props)
ax4.set_xticklabels([])
ax4.set_yticklabels([])
ax4.set_aspect('equal')

ax1.axis('off')
ax2.axis('off')
ax3.axis('off')
ax4.axis('off')

# Create a Rectangle patch
rect = patches.Rectangle((2000,1400), 1200, 1200, linewidth=1, edgecolor='r', facecolor='none')

# Add the patch to the Axes
ax1.add_patch(rect)
rect = patches.Rectangle((2000,1400), 1200, 1200, linewidth=1, edgecolor='r', facecolor='none')

ax2.add_patch(rect)
rect = patches.Rectangle((2000,1400), 1200, 1200, linewidth=1, edgecolor='r', facecolor='none')

ax3.add_patch(rect)
rect = patches.Rectangle((2000,1400), 1200, 1200, linewidth=1, edgecolor='r', facecolor='none')

ax4.add_patch(rect)
fig.tight_layout()



fig.show()