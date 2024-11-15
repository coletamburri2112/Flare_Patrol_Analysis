#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 12:01:35 2024

@author: coletamburri
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from specutils import Spectrum1D
from astropy import units as u

qs_wlun = u.nm
flare_wlun = u.nm
qs_intun = u.Unit('W cm-2 sr-1 AA-1') 
flare_intun = u.Unit('W cm-2 sr-1 AA-1') 



qswl = np.array(dispersion_range) * qs_wlun
qs_int = calibrated_qs * qs_intun
flare_wl = np.array(dispersion_range) *flare_wlun
flare_int = scaled_flare_time[0,:,maxindices[0]] *flare_intun

spectrumqs = Spectrum1D(spectral_axis=qs_wlun, flux=qs_intun)
spectrumflare = Spectrum1D(spectral_axis=flare_wlun, flux=flare_intun)

mask_region = (qs_intun > 397.02 * u.nm) & (qs_intun < 397.03 * u.nm)

masked_spectrum_qs = spectrum.mask_region(mask_region)



fig,ax=plt.subplots();
ax.plot(dispersion_range,scaled_flare_time[0,:,maxindices[0]],color='purple',label='DKIST/ViSP Calib.');
ax.set_xlim([396.95,397.1]);
ax.set_ylim([0e6,1.75e6]);
ax.plot(dispersion_range,bkgd_subtract_flaretime[0,:,maxindices[0]],color='red',label='Flare-Time'); 
ax.plot(dispersion_range,calibrated_qs,color='grey',label='QS');ax.axvline(792);ax.axvline(809);
ax.legend()
ax.axvline(397.04)
ax.axvline(397.055)
fig.show()