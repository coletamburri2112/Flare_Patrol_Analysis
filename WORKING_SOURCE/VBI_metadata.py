#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu July 2 2026

@author: coletamburri
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from astropy.io import fits
import os
import pandas as pd


pid='2_11'
vbifold = '/Volumes/ViSP_External/pid_'+pid+'_VBI/'
vbiexp = 'CNFRML' # 8 August 2024 X-class 

file = '/VBI_2024_08_08T20_12_32_333333_00656282_I_CNFRML_L1.fits' # 11 August 2024 X-class
savfold='/Users/coletamburri/Desktop/DKIST_Code/VBI_Destretching/'+vbiexp+'/'
filt = 'Halpha'

xtraflag = 'X_class_8_Aug_2024_Halpha_3_273'
hdu_list = fits.open(vbifold+vbiexp+file)
image=hdu_list[1].data[0,:,:]

tileSizeInput = [128, 64, 48, 24]  # original from FW

if os.path.isdir(savfold)==False:
    os.mkdir(savfold)

dir_list = os.listdir(vbifold+vbiexp)
dir_list.sort()

loc_files=dir_list[3:273]

# load header for first file

firstfile = fits.open(vbifold+'/'+vbiexp+'/'+loc_files[0])

header = firstfile[1].header

# make time array

t_vbi = []

for i in range(len(loc_files)):

    filename = loc_files[i]
    
    year = int(filename[-52:-48])
    month = int(filename[-47:-45])
    day = int(filename[-44:-42])
    
    hh1 = int(filename[-41:-39])
    mm1 = int(filename[-38:-36])
    ss1 = int(filename[-35:-33])
    ms1 = int(filename[-32:-26])

    time = datetime(year,month,day,hh1,mm1,ss1,ms1)
    t_vbi.append(time)
    

# make reported HPC arrays

xnum = np.shape(firstfile[1].data)[1]
ynum = np.shape(firstfile[1].data)[2]

x_mid = header['CRVAL1']
y_mid = header['CRVAL2']

x_delt = 0.017#header['CDELT1']
y_delt = 0.017#header['CDELT2']

x_HPC = np.arange(x_mid-x_delt*(xnum/2),x_mid+x_delt*((xnum/2)-1),x_delt)
y_HPC = np.arange(y_mid-y_delt*(ynum/2),y_mid+y_delt*((ynum/2)-1),y_delt)

# make generalized pixel arrays

x_pix = np.arange(0,4096)
y_pix = np.arange(0,4096)

# make generalized arcsec arrays

x_arcsec = 0.017*x_pix
y_arcsec = 0.017*x_pix

# make spatial coords pandas df

space_dict = {
    'x_pix_DKIST': x_pix,
    'y_pix_DKIST': y_pix,
    'x_rel_DKIST [arcsec]': x_arcsec,
    'y_rel_DKIST [arcsec]': y_arcsec,
    'x_HPC_DKIST [arcsec]': x_HPC,
    'y_HPC_DKIST [arcsec]': y_HPC,
}

df_space = pd.DataFrame(space_dict)

# make time coords pandas df

time_dict = {
    'time_VBI': t_vbi}

df_time = pd.DataFrame(time_dict)

# save space to file
spacefilename = '/Users/coletamburri/Desktop/DKIST_space_8_8_2024_20_12_32_3333_6563_I.txt'
df_space.to_csv(spacefilename, index=False)

# save time to file
timefilename = '/Users/coletamburri/Desktop/DKIST_times_8_8_2024_20_12_32_to_20_24_30.csv'
df_time.to_csv(timefilename, index=False)

# save header to file

headerfilename = '/Users/coletamburri/Desktop/DKIST_header_8_8_2024_20_12_32_3333_6563_I.txt'
header.totextfile(headerfilename, overwrite=True)

