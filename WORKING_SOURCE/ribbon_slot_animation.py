#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 06:50:07 2026

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
import os
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d



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
    hdul = fits.open(path+'/'+vbi_filenames[i])
    timesvbi.append(hdul[1].header['DATE-BEG'])
    if hdul[1].header['AO_LOCK']==True:
        friedvbi.append(fits.open(path+'/'+vbi_filenames[i])[1].header['ATMOS_R0']*100)   
        hdul.close()
    else:
        friedvbi.append(float(0))
        hdul.close()
        
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

# if derived
ccslot = np.load('/Users/coletamburri/Desktop/ccslot.npz')['ccslot']

ccslot2 = np.load('/Users/coletamburri/Desktop/ccslot2.npz')['ccslot2']


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

outfolder = '/Users/coletamburri/Desktop/powerspec_movie/'
outfolder2 = '/Users/coletamburri/Desktop/ribbontrace_movie/'

os.mkdir(outfolder)
#os.mkdir(outfolder2)

xarcsec = np.arange(np.shape(destretch[0].data[0,:,:][0])[0])*0.017
yarcsec = np.arange(np.shape(destretch[0].data[0,:,:][1])[0])*0.017

#for SV filter
kern1 = 8
kern2 = 8

poly1 = 3
poly2 = 3
for i in range(0,260,1):
    
    xmaxsch = xmaxsallall[i]
    ymaxsch = ymaxsallall[i]
    
    xmaxsch2 = xmaxsallall2[i]
    ymaxsch2 = ymaxsallall2[i]
    choicemap = 'coolwarm'

    fig2,ax2=plt.subplots(dpi=200,num=0, clear=True)

    ax2.pcolormesh(destretch[0].data[i,:,:],cmap='sdoaia304',alpha=0.9)
    ax2.plot(xmaxsch,ymaxsch,color='green',marker='.',linewidth=.5,markersize=1)
    ax2.plot(xmaxsch2,ymaxsch2,color='violet',marker='.',linewidth=.5,markersize=1)
    ax2.invert_yaxis()
    ax2.grid(alpha=0.2)
    
    ax2.set_xlim([1700,2400])
    ax2.set_ylim([2550,1000])
    ax2.set_xticks([1800,2000,2200],[int(xarcsec[1800]),int(xarcsec[2000]),int(xarcsec[2200])])
    ax2.set_yticks([2400,2000,1600,1200],[int(yarcsec[2400]),int(yarcsec[2000]),int(xarcsec[1600]),int(xarcsec[1200])])

    
    ax2.xaxis.set_minor_locator(MultipleLocator(50)) 
    ax2.yaxis.set_minor_locator(MultipleLocator(50)) 
    
    ax2.set_xlabel('VBI-X [arcsec]',fontsize=5)
    ax2.set_ylabel('VBI-Y [arcsec]',fontsize=5)
    
    ax2.set_title(onlytimevbi[47+i]+ ' UT',fontsize=6)
    ax2.tick_params(axis='x',labelsize=4.5)
    ax2.tick_params(axis='y',labelsize=4.5)
    ax2.set_aspect('equal')
    
    fig,ax=plt.subplots(3,2,dpi=200,num=1, clear=True)
    
    fig.suptitle(onlytimevbi[47+i]+ ' UT',fontsize=6, y=0.92)


    ax.flatten()[0].plot(arcsec_to_km*pix_to_arcsec*1/all_freqarr2[0,1:n2//2],\
                         np.power(all_psdarr2[i,1:n2//2],.25),c='black')
    # ax.flatten()[0].plot(arcsec_to_km*pix_to_arcsec*1/all_freqarr2[0,1:n2//2],\
    #                      savgol_filter(np.power(all_psdarr2[i,1:n2//2],.25),kern1,poly1),c='black')
        

    ax.flatten()[0].set_xscale('log')
    ax.flatten()[0].set_xlim([0,3000])
    #ax.flatten()[0].set_xlim([0,1500])
    ax.flatten()[0].grid(alpha=0.2)
    ax.flatten()[0].set_xlabel('Spatial scale [km]',fontsize=5)
    ax.flatten()[0].set_ylabel(r'$\sqrt[4]{Power}$',fontsize=5)
    ax.flatten()[0].axvline(100,linestyle='-.',linewidth=.5,color='grey')
    ax.flatten()[0].axvline(450,linestyle='-.',linewidth=.5,color='grey')
    
    ax.flatten()[0].axvspan(100,450,color='grey',alpha=0.2)
        
    ax.flatten()[5].pcolormesh(np.arange(300),\
                               arcsec_to_km*pix_to_arcsec*1/all_freqarr[0,1:n//2],\
                                   gaussian_filter1d(gaussian_filter1d(np.transpose(np.power(all_psdarr[:,1:n//2],.1)),axis=1,sigma=1),axis=0,sigma=.8),\
                                       cmap=choicemap,vmin=0.5,vmax=1)
        
    #ax.flatten()[5].invert_yaxis()

    ax.flatten()[5].set_yscale('log')
    
    ax.flatten()[5].axvline(i,linestyle='dashed',color='black',linewidth=.5)
    
    ax.flatten()[4].pcolormesh(np.arange(300),\
                               arcsec_to_km*pix_to_arcsec*1/all_freqarr2[0,1:n2//2],\
                                   gaussian_filter1d(gaussian_filter1d(np.transpose(np.power(all_psdarr2[:,1:n2//2],.1)),axis=1,sigma=1),axis=0,sigma=.8),\
                                       cmap=choicemap,vmin=0.5,vmax=1)
        
    #ax.flatten()[4].invert_yaxis()
        
    ax.flatten()[4].set_yscale('log')
    
    ax.flatten()[4].axvline(i,linestyle='dashed',color='black',linewidth=.5)
    
    ax.flatten()[3].pcolormesh(np.arange(300),\
                               np.arange(l)*pix_to_arcsec*arcsec_to_km,\
                                   np.transpose(intavginterps),cmap='afmhot');
    ax.flatten()[3].invert_yaxis()
    
    ax.flatten()[3].axvline(i,linestyle='dashed',color='darkorange',linewidth=.5)
    
    ax.flatten()[2].pcolormesh(np.arange(300),\
                               np.arange(l2)*pix_to_arcsec*arcsec_to_km,\
                                   np.transpose(intavginterps2),cmap='afmhot');
    ax.flatten()[2].invert_yaxis()
    
    ax.flatten()[2].axvline(i,linestyle='dashed',color='darkorange',linewidth=.5)
    
    ax.flatten()[1].plot(arcsec_to_km*pix_to_arcsec*1/all_freqarr[0,1:n//2],\
                         np.power(all_psdarr[i,1:n//2],.25),c='black')    
    # ax.flatten()[1].plot(arcsec_to_km*pix_to_arcsec*1/all_freqarr[0,1:n//2],\
    #                      savgol_filter(np.power(all_psdarr[i,1:n//2],.25),kern2,poly2),c='black')
    
    ax.flatten()[1].set_xscale('log')
    
    ax.flatten()[5].set_xlim([0,260])
    ax.flatten()[3].set_xlim([0,260])
    
    ax.flatten()[4].set_xlim([0,260])
    ax.flatten()[2].set_xlim([0,260])

    ax.flatten()[1].set_xlim([0,3000])
    

    ax.flatten()[1].grid(alpha=0.2)

    
    ax.flatten()[1].set_xlabel('Spatial scale [km]',fontsize=5)
    ax.flatten()[1].set_ylabel(r'$\sqrt[4]{Power}$',fontsize=5)
    
    ax.flatten()[3].set_xlabel('Time [UT]',fontsize=5)
    ax.flatten()[3].set_ylabel('Distance along ribbon [km]',fontsize=5)
    
    ax.flatten()[5].set_xlabel('Time [UT]',fontsize=5)
    ax.flatten()[5].set_ylabel(r'Spatial scale [km]',fontsize=5)
    
    ax.flatten()[2].set_xlabel('Time [UT]',fontsize=5)
    ax.flatten()[2].set_ylabel('Distance along ribbon [km]',fontsize=5)
    
    ax.flatten()[4].set_xlabel('Time [UT]',fontsize=5)
    ax.flatten()[4].set_ylabel(r'Spatial scale [km]',fontsize=5)
    
    ax.flatten()[0].tick_params(axis='x',labelsize=4.5)
    ax.flatten()[0].tick_params(axis='y',labelsize=4.5)
    
    ax.flatten()[1].tick_params(axis='x',labelsize=4.5)
    ax.flatten()[1].tick_params(axis='y',labelsize=4.5)

    ax.flatten()[3].tick_params(axis='x',labelsize=4.5)
    ax.flatten()[3].tick_params(axis='y',labelsize=4.5)

    ax.flatten()[5].tick_params(axis='x',labelsize=4.5)
    ax.flatten()[5].tick_params(axis='y',labelsize=4.5)
    
    ax.flatten()[0].set_ylim([-0,1.5])
    ax.flatten()[1].set_ylim([-0,1.5])
    
    
    ax.flatten()[2].tick_params(axis='x',labelsize=4.5)
    ax.flatten()[2].tick_params(axis='y',labelsize=4.5)

    ax.flatten()[4].tick_params(axis='x',labelsize=4.5)
    ax.flatten()[4].tick_params(axis='y',labelsize=4.5)

    ax.flatten()[3].set_xticks([0,80,160,240],[onlytimevbi[47],onlytimevbi[47+80],onlytimevbi[47+160],onlytimevbi[47+240]])
    ax.flatten()[5].set_xticks([0,80,160,240],[onlytimevbi[47],onlytimevbi[47+80],onlytimevbi[47+160],onlytimevbi[47+240]])

    ax.flatten()[2].set_xticks([0,80,160,240],[onlytimevbi[47],onlytimevbi[47+80],onlytimevbi[47+160],onlytimevbi[47+240]])
    ax.flatten()[4].set_xticks([0,80,160,240],[onlytimevbi[47],onlytimevbi[47+80],onlytimevbi[47+160],onlytimevbi[47+240]])
    
    ax.flatten()[1].axvline(100,linestyle='-.',linewidth=.5,color='grey')
    ax.flatten()[1].axvline(450,linestyle='-.',linewidth=.5,color='grey')
    
    ax.flatten()[1].axvspan(100,450,color='grey',alpha=0.2)
    ax.flatten()[4].set_ylim([90,2000])
    ax.flatten()[5].set_ylim([90,2000])

    
    fig.subplots_adjust(wspace=0.3)
    fig.subplots_adjust(hspace=0.5)

    fig2.savefig(outfolder2+str(i)+'.png')
    fig.savefig(outfolder+str(i)+'.png')

    plt.close('all')
    #ax.flatten()[4].plot(np.arange(l)*pix_to_arcsec*arcsec_to_km,intavginterps[i][:])


