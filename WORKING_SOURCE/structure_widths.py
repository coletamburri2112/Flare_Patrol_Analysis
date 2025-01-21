#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 12:25:21 2025

@author: coletamburri
"""
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
from astropy.modeling import models, fitting
from astropy.io import fits
import skimage
import scipy

path = '/Users/coletamburri/Desktop/VBI_Destretching/'
folder_vbi = 'AXXJL' # 8 August X-class flare decay phase

dir_list = os.listdir(path+folder_vbi)

fullhalpha = fits.open(path+folder_vbi+'/'+dir_list[1])

first_frame = fullhalpha[0].data[0,:,:]

fig,ax=plt.subplots(dpi=300,figsize=(10,10))
ax.pcolormesh(first_frame[1000:1500,500:1000],cmap='grey')
ax.set_aspect('equal')
ax.invert_yaxis()

xarr = np.arange(np.shape(first_frame)[0])
yarr = np.arange(np.shape(first_frame)[1])

spatial_samp = 0.017 # for vbi red at 656nm
arcsec_to_km = 725
xarr_km = xarr*spatial_samp
yarr_km = yarr*spatial_samp

XKM,YKM =np.meshgrid(xarr_km,yarr_km)

fig,ax=plt.subplots()
ax.pcolormesh(first_frame[1150:1250,700:900],cmap='grey')

ax.invert_yaxis()
ax.set_aspect('equal')
plt.show()

# extract
framezoom = first_frame[1150:1250,700:900]

# now use point and click to find the line you want...

aa = plt.ginput(2,timeout = 20)

# in pixel coordinates
x0,y0 = aa[0][0],aa[0][1]
x1,y1 = aa[1][0],aa[1][1]

length = int(np.hypot(x1-x0,y1-y0))

x, y = np.linspace(y0, y1, length), np.linspace(x0, x1, length)

zi = framezoom[x.astype(int), y.astype(int)]

fig, axes = plt.subplots(nrows=2)
axes[0].imshow(framezoom)
axes[0].plot([x0, x1], [y0, y1], 'ro-')
axes[0].axis('image')
axes[1].plot(zi)
plt.show()

profile = skimage.measure.profile_line(framezoom,[aa[0][1],aa[0][0]],[aa[1][1],aa[1][0]])
xdirection = np.arange(len(profile))*spatial_samp

#define limits of gaussian fitting

fig,ax=plt.subplots()
ax.plot(profile,'-x')
plt.show()

bb = plt.ginput(2,timeout = 20)

st=int(bb[0][0])
end=int(bb[1][0])+1



ynorm = [float(i)/max(profile[st:end])-1 for i in profile[st:end]]

g_init = models.Gaussian1D(amplitude=-1., mean=0.55, stddev=.1)
fit_g = fitting.LevMarLSQFitter()
gsmaller = fit_g(g_init, xdirection[st:end], ynorm)

plt.figure(figsize=(8,5))
plt.plot(xdirection[st:end], ynorm, 'ko')
plt.plot(xdirection[st:end], gsmaller(xdirection[st:end]), label='Gaussian')
plt.xlabel('Position')
plt.ylabel('Flux')
plt.legend(loc=2)
plt.show()

print(gsmaller.stddev*arcsec_to_km)