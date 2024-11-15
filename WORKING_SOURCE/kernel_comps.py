#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 20:55:13 2024

@author: coletamburri
"""

fig,ax=plt.subplots(3,1)

ax.flatten()[0].pcolormesh(obs_avg[146:237,:],cmap='hot')
ax.flatten()[0].contour(obs_avg_hep[146:237,:],cmap='ocean')
ax.flatten()[0].set_title('CaII H with Hep contour')

ax.flatten()[1].pcolormesh(obs_avg_hep[146:237,:],cmap='hot')
ax.flatten()[1].contour(obs_avg[146:237,:],cmap='ocean')
ax.flatten()[1].set_title('Hep int, Ca II contour')

ax.flatten()[2].contour(obs_avg_hep[146:237,:],linestyles='dashed',cmap='cool')
ax.flatten()[2].contour(obs_avg[146:237,:],cmap='copper')
ax.flatten()[2].set_title('Contours, Ca II (cool), Hep (copper)')

fig.show()

fig,ax=plt.subplots();
ax.plot(new_dispersion_range,scaled_flare_time[208,:,1352],label='ViSP calibrated',color='#44AA99');
ax.plot(new_dispersion_range,bkgd_subtract_flaretime[208,:,1352],label='non-flare',color='#0077BB');
ax.plot(new_dispersion_range,nonflare_average_avg,label='pre-flare',color='#CC3311');
ax.set_title('Temporal 208, Spatial 1352')
ax.grid()
ax.axvline(new_dispersion_range[caII_low],color='blue',linestyle='dotted')
ax.axvline(new_dispersion_range[caII_high],color='red',linestyle='dotted')
ax.axvline(new_dispersion_range[hep_low],color='blue',linestyle='dashdot')
ax.axvline(new_dispersion_range[-1],color='red',linestyle='dashdot')

fig.show()