#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 09:55:15 2024

@author: coletamburri
"""

# shift of wavelength range by inspection
end=0

# package initialize
import dkistpkg_ct as DKISTanalysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits

matplotlib.rcParams.update(matplotlib.rcParamsDefault)

# color scheme for plotting
muted = DKISTanalysis.color_muted2()

# path and file ID for ViSP data
path = '/Volumes/ViSP_External/pid_1_84/'
folder1 = 'AZNMO'
folder2 = 'ANYDJ' # for QS calibration - data from 22UT on 19 August

# list of files in directory for DKIST/ViSP
dir_list2 = DKISTanalysis.pathdef(path,folder1) #flaretime
dir_list3 = DKISTanalysis.pathdef(path,folder2) #qs

# Stonyhurst lon/lat position of the AR from JHv
lon = 58.57 #degrees
lat = -29.14 #degrees

# Stonyhurst lon/lat position of the AR from JHv for QS, ANYDJ
lon = 60 #degrees
lat = -28 #degrees

# Stonyhurst lon/lat position of the AR from JHv for QS, BGL


wl = 396.847 # central wavelength, Ca II H


fried=[]
aolock = []
for i in range(len(dir_list2)):
    i_file_raster1 = fits.open(path+folder1+'/'+dir_list2[i])  
    fried.append(i_file_raster1[1].header['ATMOS_R0'])
    aolock.append(i_file_raster1[1].header['AO_LOCK'])
    
fig,ax=plt.subplots();ax.scatter(range(len(fried)),fried,label='preflare',c='black');ax.set_xlabel('Timestep');ax.set_ylabel('$r_0$ [m]');ax.set_ylim([0,.15]);ax.legend();fig.show()
