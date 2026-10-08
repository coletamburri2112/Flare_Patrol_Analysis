#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 06:18:46 2024

@author: coletamburri
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits

import KMeansFunctions as KM
from Clusterer import SpectralClusterer

from matplotlib.collections import LineCollection

## DEFINE NUMBER OF CLUSTERS (determined by inspection for now, but this should be updated)
k_values = np.arange(2,30)

## Want to output metrics?
metricstest = 0

## Define mask cutoff (median intensity)
cutoff0=1

## KEYWORDS FOR RUN
full_scan_qs = 0 # if =1, subtract the pre-flare sun pixel-by-pixel
                 # if 0, subtract averaged "non-flare" from all pixels
adjust='no' # if 'no', sorts clusters by weighted mean (or other choice); otherwise
            # 'byhand' means that some cluster order is switched for clarity
nsteps = 148 # number of slit steps per ViSP scan
start = 0 # where does the interesting bit of the ViSP data start in loaded datacube?
n_init = 100 # number of times to initialize the k-means clustering; 
            # 10 by default (though 1 is probably ok if using k-means++ as initializer)
manyscan = 1 # if =0, only one scan of the ViSP; if =1, many
nframes = 10 # number of scans (if manyscan)
c = 299792458 # speed of light in m/s

## SPECTRAL AXIS PIXELS FOR EACH LINE
linelow = 380
linehigh = 500

## LINE CENTER IN NANOMETERS
cent = 854.21

## DEFINE SAVED FILE
filename = '/Users/coletamburri/Desktop/window_1746_arm1_Ca_II_(854.21_nm)_aligned_roi.fits'
caii_file = fits.open(filename)

# datacube shape - [map repeat, Stokes param, wavelength, ROI x (slit direction), ROI y (scan direction)]
flare_arr0 = caii_file[0].data[:,0,:,:,:]

#placeholder, will need to actually extract wavelengths
wave = np.arange(913)

## DEFINE SPATIAL LIMITS TO INCLUDE IN MASKING
startspace = 0
endspace = -1

## SELECT WHICH WAVELENGTHS FOR SPECTRAL AXIS
selwls = wave[linelow:linehigh]

# Clustering object
clusterer = SpectralClusterer(flare_arr0)
clusterer.prep_data(cutoff=cutoff0)
clusterer.clustering(k_values=k_values, metricstest=metricstest, n_init=n_init)




## CLUSTERING; DEPENDS ON WHICH METHOD
if clusterer == 'scikit':
    frame_line, mask0, km0, normprofiles_line, labels0, cc, x_mask0, y_mask0 = \
        kmeans_scikit(start,0,nsteps,startspace,
                    endspace,obs_avg_line,flare_arr2,normalize,n_clusters0,\
                        cutoff0,normflag=1)    
## REDEFINE STORED ViSP ARRAY
arr_normprofs0 = normprofiles_line


dists=[]

## SORTING METHODS; DEPENDS ON WHICH CLUSTERING METHOD B/C OUTPUT IS DIFFERENT FORMAT
if clusterer == 'scikit':
    for i in range(len(cc)):
        dists.append(find_30p_height(cc[i],find_nearest))
    relint=[]
    for i in range(len(cc)):
        relint.append(find_relint(cc[i],find_nearest))
    wm=[]
    for i in range(len(cc)):
        wm.append(find_weightmean(cc[i],find_nearest))

    inds = np.arange(len(cc))

## MAKE ORDERING INTO DATAFRAME
df = pd.DataFrame({'x':inds,'y':wm,'z':dists}) # by blue wing to core - 480 to 600

## SORT VALUES BY THE SORTING METHOD USED ABOVE
df.sort_values(by=['y'])

## SILLY PYTHON TYPE STUFF
sortedinds0 = df.sort_values(['y'])['x']
sortedwls = df.sort_values(['y'])['y']
sortedinds=np.asarray(sortedinds0).copy()

## ADJUST THE SORTING?
if adjust == 'byhand':
    #first = sortedinds[0]
    last = sortedinds[-1]
    #bluest = sortedinds[2]
    redest = sortedinds[-2]
    
    #sortedinds[0]=bluest
    #sortedinds[2]=first
    sortedinds[-2]=last
    sortedinds[-1]=redest

## MORE PYTHON TYPE STUFF
sortedwls = np.asarray(sortedwls)

distlocs = []
if clusterer == 'scikit':
    for i in range(len(labels0)):
        distlocs.append(np.where(sortedinds==labels0[i])[0][0])        
        
colors = plt.cm.turbo(np.linspace(0,1,n_clusters0))

if clusterer == 'scikit':
    groupsarr = np.asarray(labels0)
    
# remove pixels from kmeans? There should be logic here in the read_in case too...
if clusterer == 'scikit':
    for i in range(len(cc)):
        if cc[i][0] > cc[i][int(len(cc[i])/2)]:
            print(i)
            x_mask0[groupsarr==i]=0
            y_mask0[groupsarr==i]=0
            
maskind = {'x': x_mask0, 'y': y_mask0,'dist': distlocs}
df_mask = pd.DataFrame(maskind)



fig,ax=plt.subplots(figsize=(8,3),dpi=200)

ax.pcolormesh(np.transpose(frame_line),cmap='grey',alpha=1)
ax.scatter(x_mask0,y_mask0,.03,color=colors[distlocs],alpha=.6,marker='s')




#fig,ax=plt.subplots(5,int(n_clusters0/5),figsize=(10,6),dpi=200) #if hep and caii
fig,ax=plt.subplots(4,5,figsize=(10,6),dpi=200) #if hep and caii
if clusterer == 'scikit':
    arr_normprofs0 = normprofiles_line
    axes = ax.ravel()
    group_to_ind = {g: i for i, g in enumerate(sortedinds)}
    x = wave[linelow:linehigh]
    lines_per_ax = [[] for _ in axes]
    for curve, group in zip(normprofiles_line, labels0):
        ind = group_to_ind[group]
        pts = np.column_stack([x,curve])
        lines_per_ax[ind].append(pts)
    b = 0
    for a, lines in zip(axes, lines_per_ax):

        lc = LineCollection(
            lines,
            colors='black',
            linewidths = 0.5,
            alpha=0.01)    
        a.add_collection(lc)
        a.axvline(cent, linewidth=0.5, c='black')
        a.autoscale
        b+=1
     
    for i in range(n_clusters0):

        axes[i].plot(wave[linelow:linehigh],cc[sortedinds[i]],marker='*',color=colors[i],markersize=.1)
        axes[group].axvline(cent,linewidth=0.6,c='black')
        if i > 27:
            axes[i].tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=True,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=True,
            labelsize=6)
            axes[i].set_xlabel('Wavelength [nm]',fontsize=6)

        else: 
            axes[i].tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False)                
        axes[i].tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft=False)# labels along the bottom edge are off
        axes[i].text(0.95, 0.95, str(i+1), transform=axes[i].transAxes, \
             ha='right', va='top', fontsize=6, fontstyle='italic')
        axes[i].set_ylim([-0.2,1.2])
        axes[i].set_xlim([wave[linelow],wave[linehigh]])
        
        obs_wl = selwls[int(sortedwls[i])]
        
        velocity = find_velocity(cent,obs_wl)
        #axes[i].set_title(str(round(velocity/1e3,1))+r' km s$^{-1}$',fontsize=6,y=.895)


    fig.subplots_adjust(hspace=0,wspace=0)
    
    fig.show()
 


