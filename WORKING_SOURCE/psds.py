#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 16:01:53 2026

@author: coletamburri
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage import io, color

# 1. Load image and convert to grayscale
psds2 = []

for i in range(267):
    image = destretch[0].data[i,1700:2550,1700:2400]
    # 2. Compute 2D FFT
    fft_data = np.fft.fft2(image)
    
    # 3. Shift zero-frequency component to center
    fft_shifted = np.fft.fftshift(fft_data)
    
    # 4. Compute 2D Power Spectral Density (PSD)
    # PSD = |FFT|^2
    psd2D = np.abs(fft_shifted)**2
    
    h, w = image.shape
    y, x = np.indices((h, w))
    center = (h // 2, w // 2)
    r = np.sqrt((x - center[1])**2 + (y - center[0])**2)
    
    # 4. Bin the 2D PSD by radius
    r = r.astype(int)
    tbin = np.bincount(r.ravel(), psd2D.ravel())
    nr = np.bincount(r.ravel())
    radial_profile = tbin / nr
    
    psds2.append(radial_profile)

fig,ax=plt.subplots();
ax.pcolormesh(np.arange(267),nr,np.transpose(np.log10(psds2)),cmap='seismic');
ax2=ax.twinx();
ax2.plot(friedvbi[47:347]);
ax.set_yscale('log');
ax.set_ylim([25,2000]);
fig.show()


