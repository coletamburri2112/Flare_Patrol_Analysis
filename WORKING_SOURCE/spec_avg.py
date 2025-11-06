#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 14:11:31 2025

@author: coletamburri
"""

import numpy as np
import matplotlib.pyplot as plt
import dkistpkg_ct as DKIST_Analysis

file = '/Users/coletamburri/Desktop/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/ViSP_spectra_processed_11Aug24_CaII.npz'

def spec_avg(spec):
    nx,nw,ny = np.shape(spec)
    wl_avg = np.zeros(np.shape(spec)[1])
    for i in range(nw):
        wl_avg[i] = np.nanmean(spec[:,i,:])

    return wl_avg

caII_low = 570
caII_high = 730
hep_low = 730
hep_high = 900

hep_inner_low = 810
hep_inner_high = 825

caii_inner_low = 620
caii_inner_high = 680

hbeta_low = 500
hbeta_high = 650

data = np.load(file)
wl = data['wl']
flare = data['flare']
time = data['time']

choice = flare

fig,ax=plt.subplots()
ax.pcolormesh(np.transpose(choice[:,650,:]),cmap='magma')
ax.invert_yaxis()

aa = plt.ginput(2,timeout=120)

xlo, xhi, ylo, yhi = int(aa[0][0]), int(aa[1][0]), int(aa[0][1]),\
    int(aa[1][1])
    
scan = choice[xlo:xhi,:,ylo:yhi]
    

fig,ax=plt.subplots()
ax.pcolormesh(np.transpose(scan[:,650,:]),cmap='magma')
ax.invert_xaxis()
ax.invert_yaxis()
fig.show()


cc = plt.ginput(2,timeout = 120)

xlo2, xhi2, ylo2, yhi2 = int(cc[0][0]), int(cc[1][0]), int(cc[0][1]),\
    int(cc[1][1])
    
kern = scan[xhi2:xlo2,:,ylo2:yhi2]

avgspec = spec_avg(kern)

fig,ax=plt.subplots()
ax.plot(wl,avgspec)


    









