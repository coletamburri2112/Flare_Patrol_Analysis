#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 09:08:26 2025

@author: coletamburri
"""
import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt

path = '/Volumes/VBI_External/pid_2_11/'

VBIcode = 'AXXJL'
dir_list = os.listdir(path+VBIcode)
dir_list.sort()
dir_list2 = []

#for i in range(len(dir_list)):
for i in range(300):
    filename = dir_list[i+500]
    if filename[-5:] == '.fits' and '_I_' in filename:
        dir_list2.append(filename)

dir_list2.sort()

timestep=["" for x in range(len(dir_list2))]

fried=np.zeros(len(dir_list2))
for i in np.arange(0,len(dir_list2),1):
    hdul=fits.open(path+VBIcode+'/'+dir_list2[i])
    sample=hdul[1].data[0,:,:]
    timestep[i]=hdul[1].header['DATE-BEG']
    AO_LOCK = hdul[1].header['AO_LOCK']
    if AO_LOCK == True:
        fried[i]=hdul[1].header['ATMOS_R0']*100
        print(hdul[1].header['ATMOS_R0']*100)
    else:
        fried[i]='NaN'
    hdul.close()
    
fig,ax=plt.subplots();

ax.scatter(range(len(fried)),fried,label='8 August 2024 X-class',c='black');
ax.set_xlabel('Timestep');
ax.set_ylabel('$r_0$ [cm]');
ax.set_ylim([0,15]);
ax.legend();
fig.show()




# fig,ax=plt.subplots(dpi=300)
# ax.pcolormesh(sample,cmap='grey')
# ax.set_xticklabels([])
# ax.set_yticklabels([])
# ax.set_aspect('equal')
