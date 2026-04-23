#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:09:55 2026

@author: coletamburri
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import cv2
from skimage.measure import profile_line
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from pathlib import Path
from datetime import datetime
import matplotlib.dates as mdates ## Import required library
import sys
from scipy.optimize import curve_fit
sys.path.append('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/WORKING_SOURCE/')
import dkistpkg_ct as DKISTanalysis
from scipy.ndimage import map_coordinates
from datetime import datetime, timedelta

npoints=50 # 30 for full ribbon
image = destretch[0].data[210,:,:]
fig,ax=plt.subplots(dpi=200)

ax.imshow(image,cmap='sdoaia304')
ax.set_xlim([1700,2400])
ax.set_ylim([2550,1000])

fig.show()

ccslot = plt.ginput(npoints,timeout = 120)

def intensity_along_slot(image, points):
    points_side1 = points[0:-1:2]
    points_side2 = points[1::2]
    
    all_coords = []
    all_coords1 = []
    int_avg = []
    for i in range(len(points_side1) - 1):
        x0, y0 = points_side1[i]
        x1, y1 = points_side1[i + 1]
        
        x0_2, y0_2 = points_side2[i]
        x1_2, y1_2 = points_side2[i + 1]

        # Number of samples based on segment length - use shorter side (inside of ribbon curvature)
        length = int(np.hypot(x1_2 - x0_2, y1_2 - y0_2))
        if length == 0:
            continue

        # Evenly spaced points along the segment
        t = np.linspace(0, 1, length)
        x = x0 + t * (x1 - x0)
        y = y0 + t * (y1 - y0)
        
        x1 = x0_2 + t * (x1_2 - x0_2)
        y1 = y0_2 + t * (y1_2 - y0_2)
        
        all_coords.append(np.vstack([x, y]).T)
        all_coords1.append(np.vstack([x1, y1]).T)

        
    # now average the intensity between every two points
    
        for k in range(len(x)):
            length2 = int(np.hypot(y[k] - y1[k], x[k] - x1[k]))
            t2 = np.linspace(0, 1, length2)
    
            xp = x1[k] + t2 * (x[k]-x1[k])
            yp = y1[k] + t2 * (y[k]-y1[k])
            
            coordsp = np.vstack([xp,yp]).T
            
            sample_coordp = np.vstack([coordsp[:, 1], coordsp[:, 0]])
            
            intensitiesj = []
            
            for j in range(np.shape(coordsp)[0]):
                intensitiesj.append(image[int(coordsp[j,1]),int(coordsp[j,0])])
                
            int_avg.append(np.nanmax(intensitiesj))


    coords = np.vstack(all_coords)
    coords1 = np.vstack(all_coords1)


    return coords, coords1, int_avg

int_avg_all = []
for m in range(0,300,1):
    image = destretch[0].data[m,:,:]
    coords,coords1,int_avg = intensity_along_slot(image,ccslot)
    int_avg_all.append(int_avg)
    
fig,ax=plt.subplots();ax.pcolormesh(np.transpose(int_avg_all),cmap='sdoaia304');ax.invert_yaxis();fig.show()
    
#debugging
image = destretch[0].data[210,:,:]
fig,ax=plt.subplots(dpi=200)

ax.imshow(image,cmap='sdoaia304')
ax.set_xlim([1700,2400])
ax.set_ylim([2550,1000])
ax.scatter(all_coords[0][:,0],all_coords[0][:,1],1);ax.scatter(all_coords1[0][:,0],all_coords1[0][:,1],1)
ax.set_xlim([1900,2100]);ax.set_ylim([2500,2350]);ax.scatter(coordsp[:, 0], coordsp[:, 1]);fig.show()

