#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:09:55 2026

@author: coletamburri
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import cv2
from skimage.measure import profile_line
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from pathlib import Path
from datetime import datetime
import matplotlib.dates as mdates ## Import required library
import sys
from scipy.optimize import curve_fit
sys.path.append('/Users/coletamburri/Desktop/DKIST_Code/Flare_Patrol_Analysis/WORKING_SOURCE/')
import dkistpkg_ct as DKISTanalysis
from scipy.ndimage import map_coordinates
from datetime import datetime, timedelta
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d



# Gaussian + linear background
def gaussian_linear(x, A, mu, sigma, B, C):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2)) + B*x + C

# get times and Fried parameter for entire series
vbi_filenames=[]

path = '/Volumes/ViSP_External/pid_2_11_VBI/MBVIDS'
# Specify the directory
folder = Path(path)

# List all files and folders
for item in folder.iterdir():
    if item.is_file():
        vbi_filenames.append(item.name)

vbi_filenames.sort()
del vbi_filenames[0:3]
vbi_filenames

timesvbi=[]
friedvbi=[]

for i in range(len(vbi_filenames)-1):
    timesvbi.append(fits.open(path+'/'+vbi_filenames[i])[1].header['DATE-BEG'])
    if fits.open(path+'/'+vbi_filenames[i])[1].header['AO_LOCK']==True:
        friedvbi.append(fits.open(path+'/'+vbi_filenames[i])[1].header['ATMOS_R0']*100)    
    else:
        friedvbi.append(float(0))
        
#convert to datetime
datetimevbi=[]
for i in range(len(timesvbi)):
    date_str = timesvbi[i]
    # Format: Year-Month-Day Hour:Minute:Second
    datetimevbi.append(datetime.strptime(date_str, "%Y-%m-%dT%H:%M:%S.%f"))
    
#strip day
onlytimevbi = []
for i in range(len(timesvbi)):
    date_str = timesvbi[i]
    # Format: Year-Month-Day Hour:Minute:Second
    onlytimevbi.append(date_str[-15:-7])
    
onlytimevbi_samp = onlytimevbi[47:347]
    
file = '/Users/coletamburri/Desktop/DKIST_Code/VBI_Destretching/MBVIDS/postdestretch_dataCube_Halpha_C_class_impulsive_phase_Halpha_47_347.fits'

destretch = fits.open(file)

def intensity_along_slot(image, points):
    points_side1 = points[0:-1:2]
    points_side2 = points[1::2]
    
    all_coords = []
    all_coords1 = []
    int_avg = []
    xmaxsall = []
    ymaxsall = []
    xmidsall = []
    ymidsall = []
    for i in range(len(points_side1) - 1):
        x0, y0 = points_side1[i]
        x1, y1 = points_side1[i + 1]
        
        x0_2, y0_2 = points_side2[i]
        x1_2, y1_2 = points_side2[i + 1]

        # Number of samples based on segment length - use longer side (inside of ribbon curvature)
        length = int(np.hypot(x0_2 - x1_2, y0_2 - y1_2))
        if length == 0:
            continue

        # Evenly spaced points along the segment
        t = np.linspace(0, 1, length)
        x = x0 + t * (x1 - x0)
        y = y0 + t * (y1 - y0)
        
        x1 = x0_2 + t * (x1_2 - x0_2)
        y1 = y0_2 + t * (y1_2 - y0_2)
        
        all_coords.append(np.vstack([x, y]).T)
        all_coords1.append(np.vstack([x1, y1]).T)

        
    # now average the intensity between every two points
    
        xmaxs = []
        ymaxs = []
        xmids = []
        ymids = []
        for k in range(len(x)):
            length2 = int(np.hypot(y[k] - y1[k], x[k] - x1[k]))
            t2 = np.linspace(0, 1, length2)
            xmid = np.nanmean([x[k],x1[k]])
            ymid = np.nanmean([y[k],y1[k]])
    
            xp = x1[k] + t2 * (x[k]-x1[k])
            yp = y1[k] + t2 * (y[k]-y1[k])
            
            coordsp = np.vstack([xp,yp]).T
            
            sample_coordp = np.vstack([coordsp[:, 1], coordsp[:, 0]])
            
            intensitiesj = []
            
            for j in range(np.shape(coordsp)[0]):
                intensitiesj.append(image[int(coordsp[j,1]),int(coordsp[j,0])])
                
            int_avg.append(np.nanmax(intensitiesj)) #max shows patterns better
            
            maxind = np.where(intensitiesj==np.nanmax(intensitiesj))
            
            xmaxcoord = int(coordsp[maxind[0][0],0])
            ymaxcoord = int(coordsp[maxind[0][0],1])
            
            xmaxs.append(xmaxcoord)
            ymaxs.append(ymaxcoord)
            
            xmids.append(xmid)
            ymids.append(ymid)

        
        xmidsall.append(xmids)
        ymidsall.append(ymids)
        
        xmaxsall.append(xmaxs)
        ymaxsall.append(ymaxs)
        
        flat_xmaxs = [item for sublist in xmaxsall for item in sublist]
        flat_ymaxs = [item for sublist in ymaxsall for item in sublist]

        flat_xmids = [item for sublist in xmidsall for item in sublist]
        flat_ymids = [item for sublist in ymidsall for item in sublist]
        

    return flat_xmaxs,flat_ymaxs, flat_xmids,flat_ymids,int_avg

def calculate_distance_along_curve(x, y):
    # Calculate the difference between consecutive points
    dx = np.diff(x)
    dy = np.diff(y)
    
    # Calculate Euclidean step distance: sqrt(dx^2 + dy^2)
    step_distances = np.sqrt(dx**2 + dy**2)
    
    # Cumulative sum to get distance from start to each point
    # Prepend 0 to represent the start point (distance = 0)
    cumulative_distance = np.concatenate(([0], np.cumsum(step_distances)))
    
    
    return cumulative_distance

npoints=40# 10 for full ribbon
image = destretch[0].data[210,:,:]
fig,ax=plt.subplots(dpi=200)

ax.imshow(image,cmap='sdoaia304')
ax.set_xlim([1700,2400])
ax.set_ylim([2550,1000])

fig.show()

ccslot = plt.ginput(npoints,timeout = 120)

npoints2=10# 30 for full ribbon
image = destretch[0].data[210,:,:]
fig,ax=plt.subplots(dpi=200)

ax.imshow(image,cmap='sdoaia304')
ax.set_xlim([1700,2400])
ax.set_ylim([2550,1000])

fig.show()

ccslot2 = plt.ginput(npoints2,timeout = 120)

int_avg_all = []
xmaxsallall = []
ymaxsallall = []
xmidsallall = []
ymidsallall = []
pixdistancesall = []

int_avg_all2 = []
xmaxsallall2 = []
ymaxsallall2 = []
xmidsallall2 = []
ymidsallall2 = []
pixdistancesall2 = []

for i in range(0,300,1):
    image = destretch[0].data[i,:,:]
    xmaxsall,ymaxsall,xmidsall,ymidsall,int_avg = intensity_along_slot(image,ccslot)
    xmaxsall2,ymaxsall2,xmidsall2,ymidsall2,int_avg2 = intensity_along_slot(image,ccslot2)

    int_avg_all.append(int_avg)
    xmaxsallall.append(xmaxsall)
    ymaxsallall.append(ymaxsall)
    xmidsallall.append(xmidsall)
    ymidsallall.append(ymidsall)
    pixdistancesall.append(calculate_distance_along_curve(xmidsall,ymidsall))
    
    int_avg_all2.append(int_avg2)
    xmaxsallall2.append(xmaxsall2)
    ymaxsallall2.append(ymaxsall2)
    xmidsallall2.append(xmidsall2)
    ymidsallall2.append(ymidsall2)
    pixdistancesall2.append(calculate_distance_along_curve(xmidsall2,ymidsall2))
    
pixdistanceinterps = []
intavginterps = []

for i in range(0,300,1):
    int_avg = int_avg_all[i]
    pixdist = pixdistancesall[i]

    newlow = 0
    newhigh = np.round(np.max(pixdist))
    
    newrange = np.arange(newlow,newhigh,1)
    
    f_nearest = interp1d(pixdist, int_avg, kind='nearest', fill_value="extrapolate")

    y_new = f_nearest(newrange)
    
    pixdistanceinterps.append(newrange)
    
    intavginterps.append(y_new)

l=len(newrange)
n=len(newrange)

#compute FFT for the first

dx = l / n      # Spatial sampling interval
space_x = np.linspace(0, l, n, endpoint=False)

all_freq = []
psds = []

for i in range(len(intavginterps)):
    chintensities = intavginterps[i]
    # Compute the 1D FFT
    fft_result = np.fft.fft(chintensities/np.max(chintensities))
    #fft_result = np.fft.fft(chintensities-np.nanmean(chintensities))
    # Get the frequencies for the result
    freqs = np.fft.fftfreq(n,d=dx)
    
    all_freq.append(freqs)
    
    psds.append((np.abs(fft_result)**2)/(1/dx)/n) #power spectral density (n is number of samples, dx is sampling frequency)
    
all_psdarr = np.array(psds)
all_freqarr = np.array(all_freq)

# again, for other half of ribbon
pixdistanceinterps2 = []
intavginterps2 = []
for i in range(0,300,1):
    int_avg = int_avg_all2[i]
    pixdist = pixdistancesall2[i]

    newlow = 0
    newhigh = np.round(np.max(pixdist))
    
    newrange2 = np.arange(newlow,newhigh,1)
    
    f_nearest = interp1d(pixdist, int_avg, kind='nearest', fill_value="extrapolate")

    y_new = f_nearest(newrange2)
    
    pixdistanceinterps2.append(newrange)
    
    intavginterps2.append(y_new)

l2=len(newrange2)
n2=len(newrange2)

dx2 = l2 / n2      # Spatial sampling interval
space_x2 = np.linspace(0, l2, n2, endpoint=False)

all_freq2 = []
psds2 = []

for i in range(len(intavginterps2)):
    chintensities = intavginterps2[i]
    # Compute the 1D FFT
    fft_result = np.fft.fft(chintensities/np.max(chintensities))
    #fft_result = np.fft.fft(chintensities-np.nanmean(chintensities))
    # Get the frequencies for the result
    freqs = np.fft.fftfreq(n2,d=dx)
    
    all_freq2.append(freqs)
    
    psds2.append((np.abs(fft_result)**2)/(1/dx2)/n2) #power spectral density (n is number of samples, dx is sampling frequency)
    
all_psdarr2 = np.array(psds2)
all_freqarr2 = np.array(all_freq2)



pix_to_arcsec = 0.017
arcsec_to_km=727
km_to_Mm = 0.001
    
#loadvbilc = np.load('/Users/coletamburri/Desktop/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/vbi_lc_extended.npz',allow_pickle='True')
loadvbilc = np.load('/Users/coletamburri/Desktop/vbi_lc_extended.npz',allow_pickle='True')

#timesvbi=loadvbilc['times']
lcvbi=loadvbilc['lc']
t3 = np.arange(datetime(2024,8,11,22,31,26),
              datetime(2024,8,11,22,38,57), 
              timedelta(seconds=2.666)).astype(datetime)
choicemap = 'coolwarm'
fig,ax=plt.subplots(4,1,dpi=100,figsize=(10,15), gridspec_kw={'hspace': 0});
ax.flatten()[0].pcolormesh(np.arange(300),np.arange(l)*pix_to_arcsec*arcsec_to_km,np.transpose(intavginterps),cmap='afmhot');
ax.flatten()[0].invert_yaxis();
ax.flatten()[1].pcolormesh(np.arange(300),arcsec_to_km*pix_to_arcsec*1/all_freqarr[0,1:n//2],gaussian_filter1d(gaussian_filter1d(np.transpose(np.power(all_psdarr[:,1:n//2],.1)),axis=1,sigma=1),axis=0,sigma=.8),cmap=choicemap,vmin=0.5,vmax=1)

ax.flatten()[1].set_ylim([90,2000])
ax.flatten()[1].set_yscale('log')
ax.flatten()[1].set_xlabel('Time [UT]',fontsize=8)

ax2 = ax.flatten()[0].twinx()
ax3 = ax.flatten()[1].twinx()

ax2.plot(friedvbi[47:347],c='white',linewidth=1)
ax3.plot(friedvbi[47:347],c='black',linewidth=2)

ax4 = ax.flatten()[0].twinx()
ax4.plot(np.arange(300),lcvbi[:300],color='yellow',linewidth=2);
ax4.set_yticks([])

ax5 = ax.flatten()[1].twinx()
ax5.set_yticks([])


ax2.set_ylim([1,19])
ax3.set_ylim([1,19])


ax2.set_ylabel(r'$r_0$ [cm]',fontsize=8)
ax3.set_ylabel(r'$r_0$ [cm]',fontsize=8)

ax.flatten()[0].set_xlim([0,260])
ax.flatten()[1].set_xlim([0,260])
ax.flatten()[2].set_xlim([0,260])
ax.flatten()[3].set_xlim([0,260])


ax.flatten()[0].tick_params(axis='both', labelsize=8)
ax.flatten()[1].tick_params(axis='both', labelsize=8)
ax.flatten()[2].tick_params(axis='both', labelsize=8)
ax.flatten()[3].tick_params(axis='both', labelsize=8)
ax2.set_xticks([])

ax2.tick_params(axis='both', labelsize=8)
ax3.tick_params(axis='both', labelsize=8)
ax4.tick_params(axis='both', labelsize=8)
ax5.tick_params(axis='both', labelsize=8)

ax.flatten()[2].pcolormesh(np.arange(300),np.arange(l2)*pix_to_arcsec*arcsec_to_km,np.transpose(intavginterps2),cmap='afmhot');
ax.flatten()[2].invert_yaxis();
ax.flatten()[3].pcolormesh(np.arange(300),arcsec_to_km*pix_to_arcsec*1/all_freqarr2[0,1:n2//2],gaussian_filter1d(gaussian_filter1d(np.transpose(np.power(all_psdarr2[:,1:n2//2],.1)),axis=1,sigma=1),axis=0,sigma=.75),cmap=choicemap,vmin=0.5,vmax=1)

ax.flatten()[3].set_ylim([90,2000])
ax.flatten()[3].set_yscale('log')
ax.flatten()[3].set_xlabel('Time [UT]',fontsize=8)

ax7 = ax.flatten()[3].twinx()

ax7.plot(friedvbi[47:347],c='black',linewidth=2)


ax7.set_ylim([1,19])


ax7.set_ylabel(r'$r_0$ [cm]',fontsize=8)

# vertical lines at peak of HXR... 22:33:13 UT
ax.flatten()[0].axvline(40+130,color='white',linestyle='dashed',linewidth=1)
ax.flatten()[1].axvline(40+130,color='black',linestyle='dashed',linewidth=2)
#ax.flatten()[2].axvline(40,color='white',linestyle='dashed',linewidth=1)
ax.flatten()[3].axvline(40+130,color='black',linestyle='dashed',linewidth=2)

# ... and peak of SXR - 22:35:45 UT (these found via analysis of light curve plots)
ax.flatten()[0].axvline(97+130,color='white',linestyle='dotted',linewidth=1)
ax.flatten()[1].axvline(97+130,color='black',linestyle='dotted',linewidth=2)
#ax.flatten()[2].axvline(97,color='white',linestyle='dotted',linewidth=1)
ax.flatten()[3].axvline(97+130,color='black',linestyle='dotted',linewidth=2)


ax.flatten()[1].set_xticks([])
ax.flatten()[2].set_xticks([])
ax.flatten()[3].set_xticks([])

ax.flatten()[3].set_xticks([0,40,80,120,160,200,240],[onlytimevbi[47],onlytimevbi[47+40],onlytimevbi[47+80],onlytimevbi[47+120],onlytimevbi[47+160],onlytimevbi[47+200],onlytimevbi[47+240]])


ax.flatten()[0].tick_params(axis='both', labelsize=8)
ax.flatten()[1].tick_params(axis='both', labelsize=8)
ax.flatten()[2].tick_params(axis='both', labelsize=8)
ax.flatten()[3].tick_params(axis='both', labelsize=8)

ax2.tick_params(axis='both', labelsize=8)
ax3.tick_params(axis='both', labelsize=8)
ax4.tick_params(axis='both', labelsize=8)
ax5.tick_params(axis='both', labelsize=8)
ax7.tick_params(axis='both', labelsize=8)

ax.flatten()[1].set_ylabel('Spatial scale [km]',fontsize=8)
#ax.flatten()[1].set_ylabel(r'Frequency [pix$^{-1}$]',fontsize=8)

ax.flatten()[1].set_xlabel('Time [UT]',fontsize=8)

ax.flatten()[3].set_ylabel('Spatial scale [km]',fontsize=8)

ax.flatten()[3].set_xlabel('Time [UT]',fontsize=8)


ax.flatten()[2].set_ylabel('Distance along trace [km]',fontsize=8)
ax.flatten()[2].set_xlabel('Time [UT]',fontsize=8)


ax.flatten()[0].set_ylabel('Distance along trace [km]',fontsize=8)
ax.flatten()[0].set_xlabel('Time [UT]',fontsize=8)
ax.flatten()[0].set_yticks([0,5000,10000])
ax.flatten()[2].set_yticks([500,1500,2500])

for i in range(len(friedvbi[47:347])):
    if friedvbi[47+i]<3:
        ax.flatten()[1].axvline(i, color='lavender', alpha=1,linewidth=3)
        ax.flatten()[3].axvline(i, color='lavender', alpha=1,linewidth=3)
fig.show()


xs_full = []
ys_full = []

for i in range(npoints):
    xs_full.append(ccslot[i][0])
    ys_full.append(ccslot[i][1])



xs_part = []
ys_part = []

for i in range(npoints2):
    xs_part.append(ccslot2[i][0])
    ys_part.append(ccslot2[i][1])



xarcsec = np.arange(np.shape(destretch[0].data[0,:,:][0])[0])*0.017
yarcsec = np.arange(np.shape(destretch[0].data[0,:,:][1])[0])*0.017

fig,ax=plt.subplots(dpi=200,figsize=(3,10))

ax.pcolormesh(destretch[0].data[100+130,:,:],cmap='sdoaia304',alpha=0.9)
ax.plot(xs_full[0::2],ys_full[0::2],c='green',linestyle='-',marker='o',markersize=1,linewidth=.75)
ax.plot(xs_part[0::2],ys_part[0::2],c='violet',linestyle='-',marker='o',markersize=1,linewidth=.75)
ax.plot(xs_full[1::2],ys_full[1::2],c='green',linestyle='-',marker='o',markersize=1,linewidth=.75)
ax.plot(xs_part[1::2],ys_part[1::2],c='violet',linestyle='-',marker='o',markersize=1,linewidth=.75)

ax.plot(xs_part[0:2],ys_part[0:2],c='violet',linestyle='-',marker='o',markersize=1,linewidth=.75)
ax.plot(xs_part[-2:],ys_part[-2:],c='violet',linestyle='-',marker='o',markersize=1,linewidth=.75)

ax.plot(xs_full[0:2],ys_full[0:2],c='green',linestyle='-',marker='o',markersize=1,linewidth=.75)
ax.plot(xs_full[-2:],ys_full[-2:],c='green',linestyle='-',marker='o',markersize=1,linewidth=.75)


ax.set_xlim([1700,2400])
ax.set_ylim([2550,1000])
ax.set_xticks([1800,2000,2200],[int(xarcsec[1800]),int(xarcsec[2000]),int(xarcsec[2200])])
ax.set_yticks([2400,2000,1600,1200],[int(yarcsec[2400]),int(yarcsec[2000]),int(xarcsec[1600]),int(xarcsec[1200])])


ax.tick_params(axis='x',labelsize=5.5)
ax.tick_params(axis='y',labelsize=5.5)
ax.set_aspect('equal')

ax.xaxis.set_minor_locator(MultipleLocator(50)) 
ax.yaxis.set_minor_locator(MultipleLocator(50)) 
ax.set_xlabel('VBI-X [arcsec]',fontsize=6)
ax.set_ylabel('VBI-Y [arcsec]',fontsize=6)

for i in range(2,len(xs_full),2):
    ax.plot(xs_full[i:i+2],ys_full[i:i+2],c='green',linestyle='--',marker='o',markersize=1,linewidth=.75)

for i in range(2,len(xs_part),2):
    ax.plot(xs_part[i:i+2],ys_part[i:i+2],c='violet',linestyle='--',marker='o',markersize=1,linewidth=.75)

ax.grid(alpha=0.2)

fig.show()