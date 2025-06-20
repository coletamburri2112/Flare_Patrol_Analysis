#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 19:27:59 2024

@author: coletamburri
"""
import numpy as np
import matplotlib.pyplot as plt

n_points = 14
shiftx = 0
shifty = 0
num = 30

# for upper portion
x_inds = np.linspace(195+shiftx,214+shiftx,n_points)
y_inds = np.linspace(1270+shifty,1414+shifty,n_points)

# for middle portion
#x_inds = np.linspace(217,237,n_points)
#y_inds = np.full((1,n_points),1100+shifty)[0]

# for bottom portion
#x_inds = np.linspace(200+shiftx,220+shiftx,n_points)
#y_inds = np.linspace(520+shifty,550+shifty,n_points)

colors = plt.cm.turbo(np.linspace(0,1,n_points))

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(10,5),dpi=150)

for i in range(len(x_inds)):
    ax1.plot(new_dispersion_range,scaled_flare_time[int(x_inds[i]),:,int(y_inds[i])],c=colors[i],
            linewidth=3)
    
ax1.set_xlim([396.7,397.1])
ax1.grid()
ax1.set_xlabel('Wavelength [nm]')
ax1.set_ylabel('Intensity [$erg\;s^{-1}\;sr^{-1}\;\AA^{-1}\;cm^{-2}$]')
ax1.set_xticks([396.8,396.9,397.0,397.1])
    

ax2.pcolormesh(np.transpose(caii_avg),cmap='magma',alpha=.6)
ax2.set_xlim([183,241])
ax2.set_ylim([400,1600])
ax2.scatter(x_inds,y_inds,15,c=colors,edgecolors='black')
ax2.set_aspect(.15)

fig.show()
fig.savefig('/Users/coletamburri/Desktop/'+str(num)+'.png')