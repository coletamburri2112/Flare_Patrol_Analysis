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

widtherrs1 = sample1['arr_1']
widtherrs2 = sample2['arr_1']

muted = tc.tol_cset('muted')

#os.mkdir('/Users/coletamburri/Desktop/fine_stats/')

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





