#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 12:58:27 2025

@author: coletamburri
"""
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import skimage
import scipy
import tol_colors as tc
import matplotlib.patches as patches

root = '/Users/coletamburri/Desktop/'
folder1 = 'small_loop_frame4/'
folder2 = 'small_loop_frame5/'
filename = 'widths_errors.npz'

sample1 = np.load(root+folder1+filename)
sample2 = np.load(root+folder2+filename)

#if gauss1
#arr1 = widthss
#arr2 = widtherrs
#arr3 = startx
#arr4 = starty
#arr5 = endx
#arr6 = endy
#arr7 = r2s
#arr8 = amps
#arr9 = note
#arr10 = time
#arr11 = ylos
#arr12 = yhis
#arr13 = xlos
#arr14 = xhis

def min_max_normalize(data):
  """
  Normalizes a list of numbers to the range [0, 1] using min-max scaling.

  Args:
    data: A list of numbers.

  Returns:
    A new list with normalized values.
  """


  min_val = min(data)
  max_val = max(data)

  if min_val == max_val:
    return [0.0] * len(data)

  normalized_data = [(x - min_val) / (max_val - min_val) for x in data]
  return normalized_data

widths1 = sample1['arr_0']
widths2 = sample2['arr_0']
amps1 = sample1['arr_2']
amps2 = sample2['arr_2']

# ylo1 = sample1['arr_10']
# yhi1 = sample1['arr_11']
# xlo1 = sample1['arr_12']
# xhi1 = sample1['arr_13']

# ylo2 = sample2['arr_10']
# yhi2 = sample2['arr_11']
# xlo2 = sample2['arr_12']
# xhi2 = sample2['arr_13']

widtherrs1 = sample1['arr_1']
widtherrs2 = sample2['arr_1']

muted = tc.tol_cset('muted')

os.mkdir('/Users/coletamburri/Desktop/fine_stats/')

fig,ax=plt.subplots(dpi=200)
ax.errorbar(range(len(widths1)),widths1,widtherrs1,linestyle='',fmt='.',
             ecolor=muted.indigo,elinewidth=2,capsize=2,zorder=0)
ax.scatter(range(len(widths1)),widths1,c=min_max_normalize(amps1),edgecolors='black',\
           zorder=1,cmap='Blues',label='Frame 4')
ax.errorbar(range(len(widths2)),widths2,widtherrs2,linestyle='',fmt='.',
             ecolor=muted.rose,elinewidth=2,capsize=2,zorder=0)
ax.scatter(range(len(widths2)),widths2,c=min_max_normalize(amps2),marker='s',edgecolors='black',\
           zorder=1,cmap='Reds',label='Frame 5')
ax.axhline(np.nanmedian(widths1),linestyle='--',c=muted.indigo,linewidth=4,label='Frame 4 median')
ax.axhline(np.nanmedian(widths2),linestyle='--',c=muted.rose,linewidth=4,label = 'Frame 5 median') 
ax.set_ylabel('Width [km]',fontsize=12,font='Tahoma')
ax.set_xticks([])

ax.axvspan(0, 10, alpha=0.5, color='black',zorder=0)
ax.axvspan(20, 30, alpha=0.5, color='black',zorder=0)
ax.axvspan(40, 50, alpha=0.5, color='black',zorder=0)



ax.legend()
fig.show()

fig.savefig('/Users/coletamburri/Desktop/fine_stats/scatter.png')

bins=np.arange(20,140,10)

fig,ax=plt.subplots(dpi=200)
ax.hist(widths1,edgecolor='k',alpha=0.3,bins=bins,zorder=1)
ax.hist(widths2,edgecolor='k',alpha=0.3,bins=bins,zorder=0)
ax.axvline(np.nanmedian(widths1),linestyle='--',c=muted.indigo,linewidth=4,label='Frame 4 median')
ax.axvline(np.nanmedian(widths2),linestyle='--',c=muted.rose,linewidth=4,label = 'Frame 5 median') 
ax.set_xlabel('Width [km]',fontsize=12,font='Tahoma')
ax.legend()
fig.show()

fig.savefig('/Users/coletamburri/Desktop/fine_stats/histo.png')

# #show cuts
# # for npz loading
# path = '/Users/coletamburri/Desktop/VBI_Destretching/'
# folder_vbi = 'AXXJL/' # 8 August X-class flare decay phase
# filename = 'AXXJLselections.npz'
# array = np.load(path+folder_vbi+filename)['first50']

# #frame to work with
# frame = array[4,:,:]

# # Constants
# spatial_samp = 0.017 # for vbi red at 656nm
# arcsec_to_km = 727 # approximate arcsec to km conversion

# # X and Y coordinates of frame
# xarr = np.arange(np.shape(frame)[0])
# yarr = np.arange(np.shape(frame)[1])

# # X and Y coordinates, in KM
# xarr_km = xarr*spatial_samp
# yarr_km = yarr*spatial_samp

# # Meshgrid for plotting
# XKM,YKM =np.meshgrid(xarr_km,yarr_km)

# # Plot first frame
# fig,ax=plt.subplots(dpi=200,figsize=(10,10))
# ax.imshow(frame,cmap='grey')
# ax.set_aspect('equal')
# for i in range(len(ylo1)):
#     rect = patches.Rectangle((ylo1[i], xlo1[i]), yhi1[i]-ylo1[i], xhi1[i]-xlo1[i], linewidth=1, edgecolor='r', \
#                              facecolor='none')
#     ax.add_patch(rect)
# ax.set_xticks([])
# ax.set_yticks([])
# plt.show()

# fig.savefig('/Users/coletamburri/Desktop/fine_stats/frame4_context')

# #frame to work with
# frame = array[5,:,:]

# # Plot first frame
# fig,ax=plt.subplots(dpi=200,figsize=(10,10))
# ax.imshow(frame,cmap='grey')
# ax.set_aspect('equal')
# for i in range(len(ylo1)):
#     rect = patches.Rectangle((ylo2[i], xlo2[i]), yhi2[i]-ylo2[i], xhi2[i]-xlo2[i], linewidth=1, edgecolor='r', \
#                              facecolor='none')
#     ax.add_patch(rect)
# ax.set_xticks([])
# ax.set_yticks([])
# plt.show()

# fig.savefig('/Users/coletamburri/Desktop/fine_stats/frame5_context')










