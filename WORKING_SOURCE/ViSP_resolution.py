#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 11:29:48 2025

@author: coletamburri
"""

import numpy as np
import matplotlib.pyplot as plt
import tol_colors as tc

from scipy.signal.windows import hann
from scipy.fft import fft
from scipy.signal import medfilt

# load spec

spectra = image_data_arr_arr

# first based on continuum, select best-seeing frame
cont = spectra[148:239,290,:]
caIIcore = spectra[148:239,650,:]

step = 0.0295 #arm1

# first compute RMS intensity contract
def RMS_contrast(spec):
    rms_arr = []
    for i in range(np.shape(spec)[0]):
        arrch = spec[i,:]
        std = np.nanstd(arrch)
        mean = np.nanmean(arrch)
        rms_arr.append(100 * std/mean)
    return rms_arr

RMS_contrast_cont = RMS_contrast(cont)
plt.plot(RMS_contrast_cont)

RMS_contrast_core = RMS_contrast(caIIcore)
plt.plot(RMS_contrast_core)


ps = []
for i in range(148,239,1):
    arrch = image_data_arr_arr[i,650,:] #pid_2_11 ca II red cont
    #arrch = image_data_arr_arrhbeta[i,150,:] # pid_2_11 hbeta red cont

    hannlen = hann(len(arrch))
    windowed = arrch*hannlen
    windowed_centered = windowed-np.nanmean(windowed)
    fft_result = np.fft.fft(windowed_centered)
    #fft_result = np.fft.fft(arrch)
    power_spectrum = np.abs(fft_result)**2
    frequencies = np.fft.fftfreq(len(windowed_centered),0.0295)
    pos_freq_ind = frequencies >=0
    positive_freq = frequencies[pos_freq_ind]
    positive_ps = power_spectrum[pos_freq_ind]
    ps.append(positive_ps)

summed = np.nansum(ps,0)/np.max(ps)/np.shape(ps)[0]

plt.figure(figsize=(10, 6))
plt.scatter(positive_freq, summed,c='black',s=3)
plt.scatter(positive_freq[39:], medfilt(summed,39)[39:],c='red',s=20)
plt.yscale('log')
plt.title('11 August 2024 - continuum PSD')
plt.xlabel(r'Frequency ($arcsec^{-1}$)')
plt.ylabel('Power')
plt.grid(True)
plt.ylim([1e-7,10])
plt.xlim([0,10])
plt.show()