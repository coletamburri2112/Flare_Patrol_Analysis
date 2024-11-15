#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 14:01:04 2024

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
from scipy.ndimage import gaussian_filter
import cv2
from shapely import Polygon

import os

# import sunpy.coordinates
# import sunpy.map
# from sunpy.net import Fido
# from sunpy.net import attrs as a

from astropy.io import fits

path = '/Volumes/VBI_External/'
folder_vbi = 'AXXJL'

dir_list = os.listdir(path+folder_vbi)

dir_list.sort()
dir_list.pop(0)
dir_list.pop(0)
dir_list[0]

alpha = 10
#sigma = 5
mult=1.25
length=30

polyor=19
threshperc = 0.5
levels=[threshperc]

minr = 700
maxr = 1650

mass_distances = []

mass_ivs= []

mass_dvs = []

mass_idxs = []

def area(vs):
    a = 0
    x0,y0 = vs[0]
    for [x1,y1] in vs[1:]:
        dx = x1-x0
        dy = y1-y0
        a += 0.5*(y0*dx - x0*dy)
        x0 = x1
        y0 = y1
    return a

def distance(x, y, x0, y0):
    """
    Return distance between point
    P[x0,y0] and a curve (x,y)
    """
    d_x = x - x0
    d_y = y - y0
    dis = np.sqrt( d_x**2 + d_y**2 )
    
    return dis

def min_distance(x, y, P, precision=5):
    """
    Compute minimum/a distance/s between
    a point P[x0,y0] and a curve (x,y)
    rounded at `precision`.
    
    ARGS:
        x, y      (array)
        P         (tuple)
        precision (int)
        
    Returns min indexes and distances array.
    """
    # compute distance
    d = distance(x, y, P[0], P[1])
    d = np.round(d, precision)
    # find the minima
    glob_min_idxs = np.argwhere(d==np.min(d)).ravel()
    
    idx = glob_min_idxs[0]

    if np.sqrt(x[idx]**2+y[idx]**2) > np.sqrt(P[0]**2 + P[1]**2):
        mindis = d[idx]
    else:
        mindis = -d[idx]
    return idx,mindis
k=0
for i in range(length):
    print(i)
    img_file = fits.open(path+folder_vbi+'/'+dir_list[i])
    frame = img_file[1].data[0][0:1500,:]
    img_file.close()
    framecopy = frame
    highcontrast = framecopy*alpha
    #denoised_image = gaussian_filter(highcontrast, sigma=sigma)
    
    _min, _max = np.amin(highcontrast), np.amax(highcontrast)
    # fig,ax=plt.subplots(2,1)
    # ax.flatten()[0].pcolormesh(highcontrast, vmin = _min, vmax = _max)
    # ax.flatten()[1].pcolormesh(frame, vmin = _min, vmax = _max)
    
    
    # ax.flatten()[0].invert_yaxis()
    # ax.flatten()[0].set_aspect('equal')
    
    # ax.flatten()[1].invert_yaxis()
    # ax.flatten()[1].set_aspect('equal')
    
    masklim = mult*np.median(highcontrast)
    firstmask = np.copy(highcontrast)
    firstmask[firstmask < masklim] = 0
    firstmask[firstmask > masklim] = 1
    
    fig,ax=plt.subplots(figsize=(10,5),dpi=200)
    levels = [0.5]
    CS = ax.contour(firstmask,levels=levels,linewidths=.1)
    # ax.pcolormesh(frame,alpha=.3,cmap='hot')
    
    very_denoised = gaussian_filter(highcontrast, sigma=30)
    # plt.pcolormesh(very_denoised,cmap='hot')
    
    for i in range(len(levels)):
        contour = CS.collections[i]
        vs = contour.get_paths()[0].vertices
        # Compute area enclosed by vertices.
        a = area(vs)
        #print ("r = " + str(levels[i]) + ": a =" + str(a))
        
    thresh = threshperc * np.amax(very_denoised)
    
    # Isolate pixels certainly within the mask
    xc, yc = np.where(very_denoised > thresh)
    y = np.linspace(0, 4096, 4096)
    x = np.linspace(0, 1500, 1500)
    # Fitting of fourth-order polynomial to chosen pixels and generation of
    # PIL polynomial arrays
    coeffs = np.polyfit(y[yc], x[xc], polyor)
    
    ivs = y[yc]
    
    dvs = 0
    
    for i in range(len(coeffs)):
        dvs += coeffs[i] * ivs**(polyor - i)
        
    
    biggest = 0 
    for i in range(len(CS.allsegs)):
        for j in range(len(CS.allsegs[i])):
            dat = CS.allsegs[i][j]
            polygon = Polygon(dat)
            size = polygon.area
            if size > biggest:
                biggest=size
                biggesti = i
                biggestj = j
    
    dat0 = CS.allsegs[biggesti][biggestj]
    biggestx = dat0[:,0]
    biggesty = dat0[:,1]
    
    plt.scatter(biggestx,biggesty,.1)
    
    dat0= CS.allsegs[biggesti][biggestj]
    # fig,ax=plt.subplots()
    # ax.plot(dat0[:,0],dat0[:,1])
    # ax.scatter(ivs,dvs,1,c='cyan')
    # ax.axvline(700)
    # ax.axvline(1800)
    # ax.set_xlim([700,1800])
    
    line=[]
    for i in range(len(ivs)):
        line.append((ivs[i],dvs[i]))
    
    
    poly=[]
    for i in range(len(biggestx)):
        poly.append((biggestx[i],biggesty[i]))
        
    distances = []
    idxs = []

    for i in range(len(poly)):
        point_on_poly = np.asarray(poly[i])
    
        if point_on_poly[0] > minr and point_on_poly[0] < maxr:
            min_idx, dis = min_distance(ivs,dvs, point_on_poly,precision=5)
            distances.append(dis)
            idxs.append(min_idx)
            
    # fig,ax=plt.subplots(figsize=(8,3))
    # ax.plot(ivs[idxs],distances)
    # #ax.set_ylim([0,250])
    # ax.set_xlim([minr,maxr])
    
    distances1 = [np.float32(j) for j in distances]
    ivs1 = [np.float32(j) for j in ivs]
    dvs1 = [np.float32(j) for j in dvs]
    idxs1 = [np.float32(j) for j in idxs]
        
    
    mass_distances.append(distances1)

    mass_ivs.append(ivs1)

    mass_dvs.append(dvs1)

    mass_idxs.append(idxs1)

    
# fig,ax=plt.subplots();

# for i in range(20):
#     ax.scatter(mass_ivs[i][mass_idxs[i]],mass_distances[i],0.01)

# ax.set_xlim([800,1600])
# ax.set_ylim([0,300])

fig, ax = plt.subplots(dpi=150,figsize=(15,3))
cax = ax.scatter(np.array(mass_ivs[0])[mass_idxs[0]],mass_distances[0])
ax.set_xlim([800,1600])
ax.set_ylim([0,300])

def animate(i):
    ax.cla()

    ax.plot(np.array(mass_ivs[i])[mass_idxs[i]],mass_distances[i],0.1,color='red',linestyle='solid')
    
    ax.set_xlim([800,1600])
    ax.set_ylim([0,300])

anim = animation.FuncAnimation(fig, animate,frames=length-2)
anim.save('ribbonfrontevolve_longer.gif')
fig.show()








