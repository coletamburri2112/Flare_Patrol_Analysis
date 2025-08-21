#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  9 09:00:29 2025

@author: coletamburri
"""

import numpy as np
import matplotlib.pyplot as plt
import tol_colors as tc

#file = '/Users/coletamburri/Desktop/August_2024_DKIST_Flares/ViSPselection8AugXclass_hbeta.npz'
#file = '/Users/coletamburri/Desktop/ViSPselectionHbeta11aug.npz'
#file = '/Users/coletamburri/Desktop/Misc_DKIST/11August2024_Cclass_imp_CaII.npz'
#file = '/Users/coletamburri/Desktop/ViSPselection8AugXclass.npz'
file = '/Users/coletamburri/Desktop/ViSPselection11August24Mclass.npz'

data = np.load(file)

raw = data['raw']
wl = data['wl']
# scaled = data['scaled']
flare = data['flare']
# time = data['time']

spectra = flare

caII_low = 570
caII_high = 730
hep_low = 730
hep_high = 900

hep_inner_low = 810
hep_inner_high = 825

caii_inner_low = 620
caii_inner_high = 680

hbeta_low = 400
hbeta_high = 500

n_points = 20

caII_avg = np.mean(spectra[:,caII_low:caII_high,:],1)
hep_avg = np.mean(spectra[:,hep_low:hep_high,:],1)
hep_avg_inner = np.mean(spectra[:,hep_inner_low:hep_inner_high,:],1)

caii_avg_inner = np.mean(spectra[:,caii_inner_low:caii_inner_high,:],1)

hbeta_avg = np.mean(spectra[:,hbeta_low:hbeta_high,:],1)

both_avg = np.mean(spectra[:,caII_low:hep_high,:],1)

all_avg = np.mean(spectra,1)

choice = caii_avg_inner
fig,ax=plt.subplots()
ax.pcolormesh(np.transpose(choice),cmap='magma')
ax.invert_yaxis()
fig.show()

aa = plt.ginput(2,timeout=120)

xlo, xhi, ylo, yhi = int(aa[0][0]), int(aa[1][0]), int(aa[0][1]),\
    int(aa[1][1])
    

fig,ax=plt.subplots()
ax.pcolormesh(np.transpose(choice[xlo:xhi,ylo:yhi]),cmap='magma')
ax.invert_yaxis()
fig.show()

colors = plt.cm.jet(np.linspace(0,1,n_points))

cc = plt.ginput(n_points,timeout = 120)

#fig,ax=plt.subplots(len(cc),1)
fig,ax=plt.subplots(int(len(cc)/4),4)
for i in range(len(cc)):
    xsel,ysel = cc[i][0],cc[i][1]
    ax.flatten()[i].plot(wl,spectra[int(xsel)+xlo,:,int(ysel)+ylo],color=colors[i])
    ax.flatten()[i].axvline(396.85)
    ax.flatten()[i].axvline(397.01)
fig.show()
    

# Extract coordinates
xlo, xhi, ylo, yhi = int(aa[0][0]), int(aa[1][0]), int(aa[0][1]),\
    int(aa[1][1])
    
# Extract zoomed-in frame from coordinates
framezoom = choice[xlo:xhi,ylo:yhi]

# Plot zoomed-in
fig,ax=plt.subplots(dpi=300)
ax.pcolormesh(np.transpose(choice[xlo:xhi,ylo:yhi]),cmap='magma')
for i in range(len(cc)):
    xsel,ysel = cc[i][0],cc[i][1]
    ax.plot(xsel,ysel,'x',color=colors[i])
    
#ax.set_xlim([xlo,xhi])
#ax.set_ylim([ylo,yhi])

ax.invert_yaxis()
#ax.set_aspect('equal')
plt.show()
    
    

