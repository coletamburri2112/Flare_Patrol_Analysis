#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 11:04:40 2026

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


def intensity_along_polyline(image, points):

    points = np.array(points)
    
    all_coords = []

    for i in range(len(points) - 1):
        x0, y0 = points[i]
        x1, y1 = points[i + 1]

        # Number of samples based on segment length
        length = int(np.hypot(x1 - x0, y1 - y0))
        if length == 0:
            continue

        # Evenly spaced points along the segment
        t = np.linspace(0, 1, length)
        x = x0 + t * (x1 - x0)
        y = y0 + t * (y1 - y0)

        all_coords.append(np.vstack([x, y]).T)

    coords = np.vstack(all_coords)

    # map_coordinates expects (row, col) = (y, x)
    sample_coords = np.vstack([coords[:, 1], coords[:, 0]])

    ## if want to do via spline interpolant
    #intensities = map_coordinates(image, sample_coords, order=3, mode='nearest')
    
    # try wihtout interpolant
    intensities=[]
    for i in range(np.shape(sample_coords)[1]):
        intensities.append(image[int(sample_coords[0,i]),int(sample_coords[1,i])])

    return coords, intensities

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

ch=230

image = destretch[0].data[ch,:,:]

npoints=30 # 30 for full ribbon

fig,ax=plt.subplots(dpi=200)
ax.imshow(image,cmap='sdoaia304')
ax.set_xlim([1700,2400])
ax.set_ylim([2550,1000])

fig.show()

# if original
ccfull = plt.ginput(npoints,timeout = 120)

xs_full = []
ys_full = []

for i in range(npoints):
    xs_full.append(ccfull[i][0])
    ys_full.append(ccfull[i][1])
    
npoints=10 # 10 for part ribbon

fig,ax=plt.subplots(dpi=200)
ax.imshow(image,cmap='sdoaia304')
ax.set_xlim([1700,2400])
ax.set_ylim([2550,1000])

fig.show()

# if original
ccpart = plt.ginput(npoints,timeout = 120)

xs_part = []
ys_part = []

for i in range(npoints):
    xs_part.append(ccpart[i][0])
    ys_part.append(ccpart[i][1])



# # if derived
# ccfull = np.load('/Users/coletamburri/Desktop/fullribboncoords_morepreflare.npz')['points']

# ccpart = np.load('/Users/coletamburri/Desktop/lowerribboncoords_morepreflare.npz')['points']

# xs_full = []
# ys_full = []

# for i in range(len(ccfull)):
#     xs_full.append(ccfull[i][0])
#     ys_full.append(ccfull[i][1])
    
# xs_part = []
# ys_part = []

# for i in range(len(ccpart)):
#     xs_part.append(ccpart[i][0])
#     ys_part.append(ccpart[i][1])


intensities_all = []
coords_all=[]

for i in range(0,300,1):
    image = destretch[0].data[i,:,:]
    coords, intensities = intensity_along_polyline(image, ccfull)
    intensities_all.append(intensities)
    coords_all.append(coords)
    
l=len(intensities)
n=len(intensities)
    
    
intensities_part = []
coords_part=[]

for i in range(0,300,1):
    image = destretch[0].data[i,:,:]
    coords, intensities = intensity_along_polyline(image, ccpart)
    intensities_part.append(intensities)
    coords_part.append(coords)
    
l2=len(intensities)
n2=len(intensities)

#compute FFT for the first

dx = l / n      # Spatial sampling interval
space_x = np.linspace(0, l, n, endpoint=False)

all_freq = []
psds = []

for i in range(len(intensities_all)):
    chintensities = intensities_all[i]
    # Compute the 1D FFT
    fft_result = np.fft.fft(chintensities/np.max(chintensities))
    #fft_result = np.fft.fft(chintensities-np.nanmean(chintensities))
    # Get the frequencies for the result
    freqs = np.fft.fftfreq(n,d=dx)
    
    all_freq.append(freqs)
    
    psds.append((np.abs(fft_result)**2)/(1/dx)/n) #power spectral density (n is number of samples, dx is sampling frequency)
    
all_psdarr = np.array(psds)
all_freqarr = np.array(all_freq)

pix_to_arcsec = 0.017
arcsec_to_km=727
km_to_Mm = 0.001

# FFT/PSD for the smaller sample

dx2 = l2 / n2      # Spatial sampling interval
space_x2 = np.linspace(0, l2, n2, endpoint=False)

all_freq_part = []
psds_part = []

for i in range(len(intensities_all)):
    chintensities = intensities_part[i]
    # Compute the 1D FFT
    fft_result = np.fft.fft(chintensities/np.max(chintensities))
    #fft_result = np.fft.fft(chintensities-np.nanmean(chintensities))
    # Get the frequencies for the result
    freqs = np.fft.fftfreq(n2,d=dx)
    
    all_freq_part.append(freqs)
    
    psds_part.append((np.abs(fft_result)**2)/(1/dx2)/n2) #power spectral density (n is number of samples, dx is sampling frequency)
    
part_psdarr = np.array(psds_part)
part_freqarr = np.array(all_freq_part)

# fig,ax=plt.subplots(1,3,dpi=200)
# ax.flatten()[0].imshow(destretch[0].data[20,:,:],cmap='sdoaia304')
# ax.flatten()[0].plot(xs,ys,c='grey',linestyle='solid',markersize=4,linewidth=1)
# ax.flatten()[0].set_xlim([1700,2400])
# ax.flatten()[0].set_ylim([2550,1000])

# ax.flatten()[1].imshow(destretch[0].data[38,:,:],cmap='sdoaia304')
# ax.flatten()[1].plot(xs,ys,c='grey',linestyle='solid',markersize=4,linewidth=1)
# ax.flatten()[1].set_xlim([1700,2400])
# ax.flatten()[1].set_ylim([2550,1000])

# ax.flatten()[2].imshow(destretch[0].data[100,:,:],cmap='sdoaia304')
# ax.flatten()[2].plot(xs,ys,c='grey',linestyle='solid',markersize=4,linewidth=1)
# ax.flatten()[2].set_xlim([1700,2400])
# ax.flatten()[2].set_ylim([2550,1000])

# ax.flatten()[0].set_xticks([])
# ax.flatten()[1].set_xticks([])
# ax.flatten()[2].set_xticks([])

# ax.flatten()[0].set_yticks([])
# ax.flatten()[1].set_yticks([])
# ax.flatten()[2].set_yticks([])


# fig.show()

# plotting for both traces of ribbon - full ribbon and partial ribbon

X=np.arange(np.shape(destretch[0].data[0,:,:])[0])*0.017
Y=np.arange(np.shape(destretch[0].data[0,:,:])[0])*0.017

# fig,ax=plt.subplots(1,2,dpi=200,figsize=(3,20))

# ax.flatten()[0].pcolormesh(destretch[0].data[38+130,:,:],cmap='gray')
# ax.flatten()[0].plot(xs_full,ys_full,c='green',linestyle='solid',markersize=4,linewidth=2)
# ax.flatten()[0].plot(xs_part,ys_part,c='violet',linestyle='solid',markersize=4,linewidth=2)
# ax.flatten()[0].set_xlim([1700,2400])
# ax.flatten()[0].set_ylim([2550,1000])

# xarcsec = np.arange(np.shape(destretch[0].data[0,:,:][0])[0])*0.017
# yarcsec = np.arange(np.shape(destretch[0].data[0,:,:][1])[0])*0.017

# ax.flatten()[1].pcolormesh(destretch[0].data[100+130,:,:],cmap='gray')
# ax.flatten()[1].scatter(xs_full,ys_full,2,c='green',linestyle='-o')#,markersize=4,linewidth=2)
# ax.flatten()[1].scatter(xs_part,ys_part,2,c='violet',linestyle='-o')#,markersize=4,linewidth=2)
# ax.flatten()[1].set_xlim([1700,2400])
# ax.flatten()[1].set_ylim([2550,1000])
# ax.flatten()[1].set_xticks([1800,2000,2200],[xarcsec[1800],xarcsec[2000],xarcsec[2200]])
# ax.flatten()[1].set_yticks([2400,2000,1600,1200],[yarcsec[2400],yarcsec[2000],xarcsec[1600],xarcsec[1200]])


# ax.flatten()[1].tick_params(axis='x',labelsize=5)
# fig.show()

fig,ax=plt.subplots(dpi=200,figsize=(3,10))

xarcsec = np.arange(np.shape(destretch[0].data[0,:,:][0])[0])*0.017
yarcsec = np.arange(np.shape(destretch[0].data[0,:,:][1])[0])*0.017

ax.pcolormesh(destretch[0].data[100+130,:,:],cmap='sdoaia304',alpha=0.9)
ax.plot(xs_full,ys_full,c='green',linestyle='-',marker='o',markersize=1,linewidth=.75)
ax.plot(xs_part,ys_part,c='violet',linestyle='-',marker='o',markersize=1,linewidth=.75)
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


ax.grid(alpha=0.2)

fig.show()

# plotting time-distance plots (incl. fft)

# load light curve

#loadvbilc = np.load('/Users/coletamburri/Desktop/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/vbi_lc_extended.npz',allow_pickle='True')
loadvbilc = np.load('/Users/coletamburri/Desktop/vbi_lc_extended.npz',allow_pickle='True')

#timesvbi=loadvbilc['times']
lcvbi=loadvbilc['lc']
t3 = np.arange(datetime(2024,8,11,22,31,26),
              datetime(2024,8,11,22,38,57), 
              timedelta(seconds=2.666)).astype(datetime)
    
fig,ax=plt.subplots(4,1,dpi=100,figsize=(10,15), gridspec_kw={'hspace': 0});
ax.flatten()[0].pcolormesh(np.arange(300),np.arange(l)*pix_to_arcsec*arcsec_to_km,np.transpose(intensities_all),cmap='afmhot');
ax.flatten()[0].invert_yaxis()

ax.flatten()[0].set_ylabel('Distance along trace [km]',fontsize=8)
ax.flatten()[0].set_xlabel('Time [UT]',fontsize=8)
#ax.flatten()[1].pcolormesh(np.arange(170),all_freqarr[0,1:n//2],np.transpose(np.power(all_psdarr[:,1:n//2],.25)),cmap='seismic');

ax.flatten()[1].pcolormesh(np.arange(300),arcsec_to_km*pix_to_arcsec*1/all_freqarr[0,1:n//2],np.transpose(np.power(all_psdarr[:,1:n//2],.25)),cmap='seismic');
#ax.flatten()[1].axhspan(200,400, facecolor='grey', alpha=0.2)
#ax.flatten()[1].axhline(200,color='red')
#ax.flatten()[1].axhline(400,color='red')

ax.flatten()[2].pcolormesh(np.arange(300),np.arange(l2)*pix_to_arcsec*arcsec_to_km,np.transpose(intensities_part),cmap='afmhot');
ax.flatten()[2].invert_yaxis()

ax.flatten()[2].set_ylabel('Distance along trace [km]',fontsize=8)
ax.flatten()[2].set_xlabel('Time [UT]',fontsize=8)

ax.flatten()[3].pcolormesh(np.arange(300),arcsec_to_km*pix_to_arcsec*1/part_freqarr[0,1:n2//2],np.transpose(np.power(part_psdarr[:,1:n2//2],.25)),cmap='seismic');
#ax.flatten()[3].pcolormesh(np.arange(170),part_freqarr[0,1:n2//2],np.transpose(np.power(part_psdarr[:,1:n2//2],.25)),cmap='seismic');


#ax.flatten()[1].invert_yaxis()
#ax.flatten()[1].axhline(50,color='grey',linestyle='dashdot') # loop size
#ax.flatten()[1].axhline(300,color='white',linestyle='dashdot') # bead size
ax.flatten()[1].set_ylim([25,2000])

#ax.flatten()[1].set_yscale('log')
ax.flatten()[1].set_ylabel('Spatial scale [km]',fontsize=8)
#ax.flatten()[1].set_ylabel(r'Frequency [pix$^{-1}$]',fontsize=8)

ax.flatten()[1].set_xlabel('Time [UT]',fontsize=8)

ax.flatten()[3].set_ylim([25,2000])

# if frequency
# ax.flatten()[1].set_ylim([0,0.1])
# ax.flatten()[1].set_ylim([0,0.1])
# ax.flatten()[3].set_ylim([0,0.1])
# ax.flatten()[3].set_ylim([0,0.1])

#ax.flatten()[3].set_yscale('log')
ax.flatten()[3].set_ylabel('Spatial scale [km]',fontsize=8)

#ax.flatten()[3].set_ylabel(r'Frequency [pix$^{-1}$]',fontsize=8)
ax.flatten()[3].set_xlabel('Time [UT]',fontsize=8)

ax2 = ax.flatten()[0].twinx()
ax3 = ax.flatten()[1].twinx()

ax6 = ax.flatten()[3].twinx()

ax.flatten()[0].set_xticks([])
ax.flatten()[1].set_xticks([])
ax.flatten()[2].set_xticks([])

# vertical lines at peak of HXR... 22:33:13 UT
ax.flatten()[0].axvline(40+130,color='white',linestyle='dashed',linewidth=1)
ax.flatten()[1].axvline(40+130,color='white',linestyle='dashed',linewidth=1)
#ax.flatten()[2].axvline(40,color='white',linestyle='dashed',linewidth=1)
ax.flatten()[3].axvline(40+130,color='white',linestyle='dashed',linewidth=1)

# ... and peak of SXR - 22:35:45 UT (these found via analysis of light curve plots)
ax.flatten()[0].axvline(97+130,color='white',linestyle='dotted',linewidth=1)
ax.flatten()[1].axvline(97+130,color='white',linestyle='dotted',linewidth=1)
#ax.flatten()[2].axvline(97,color='white',linestyle='dotted',linewidth=1)
ax.flatten()[3].axvline(97+130,color='white',linestyle='dotted',linewidth=1)




ax2.plot(friedvbi[47:347],c='white',linewidth=1)
ax3.plot(friedvbi[47:347],c='white',linewidth=1)
ax6.plot(friedvbi[47:347],c='white',linewidth=1)



ax4 = ax.flatten()[0].twinx()
ax4.plot(np.arange(300),lcvbi[:300],color='black',linewidth=2);
ax4.set_yticks([])

ax5 = ax.flatten()[1].twinx()
ax5.set_yticks([])


ax2.set_ylim([1,19])
ax3.set_ylim([1,19])
ax6.set_ylim([1,19])

ax2.set_ylabel(r'$r_0$ [cm]',fontsize=8)
ax3.set_ylabel(r'$r_0$ [cm]',fontsize=8)
ax6.set_ylabel(r'$r_0$ [cm]',fontsize=8)



ax.flatten()[3].set_xticks([0,40,80,120,160,200,240,280],[onlytimevbi[47],onlytimevbi[47+40],onlytimevbi[47+80],onlytimevbi[47+120],onlytimevbi[47+160],onlytimevbi[47+200],onlytimevbi[47+240],onlytimevbi[47+280]])


ax.flatten()[0].set_xlim([0,267])
ax.flatten()[1].set_xlim([0,267])
ax.flatten()[2].set_xlim([0,267])
ax.flatten()[3].set_xlim([0,267])


ax.flatten()[0].tick_params(axis='both', labelsize=8)
ax.flatten()[1].tick_params(axis='both', labelsize=8)
ax.flatten()[2].tick_params(axis='both', labelsize=8)
ax.flatten()[3].tick_params(axis='both', labelsize=8)

ax2.tick_params(axis='both', labelsize=8)
ax3.tick_params(axis='both', labelsize=8)
ax4.tick_params(axis='both', labelsize=8)
ax5.tick_params(axis='both', labelsize=8)
ax6.tick_params(axis='both', labelsize=8)


fig.show()

# from scipy.signal import savgol_filter



# fig,ax=plt.subplots();
# ntimes = 6

# colors = plt.cm.cool(np.linspace(0,1,ntimes))
# l2=0
# for i in range(0,120,20):
#     ax.plot(arcsec_to_km*pix_to_arcsec*1/all_freqarr[0,1:n//2],savgol_filter(all_psdarr[i,1:n//2],window_length=20,polyorder=3),label=onlytimevbi[177+i],color=colors[l2])
#     l2+=1
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlim([0,2000])
# ax.axvline(200,color='red',linestyle='dashed')
# ax.axvline(400,color='red',linestyle='dashed')
# ax.axvspan(200,400, facecolor='red', alpha=0.2)
# ax.axvline(25,color='black',linestyle='dashed')
# ax.axvline(65,color='black',linestyle='dashed')
# ax.axvspan(25,65, facecolor='grey', alpha=0.2)

# ax.legend()
# fig.show()

# # do it again but add some number to each x coordinate in cc

# ccarrf = np.array(ccfull)
# for i in range(len(ccfull)):
#     ccarrf[i,0]+=0
    
# ccarrp = np.array(ccpart)
# for i in range(len(ccpart)):
#     ccarrp[i,0]+=0


# xs_full = []
# ys_full = []

# for i in range(len(ccfull)):
#     xs_full.append(ccarrf[i][0])
#     ys_full.append(ccarrf[i][1])
    
# xs_part = []
# ys_part = []

# for i in range(len(ccpart)):
#     xs_part.append(ccarrp[i][0])
#     ys_part.append(ccarrp[i][1])

# intensities_all = []
# coords_all=[]

# for i in range(0,170,1):
#     image = destretch[0].data[i,:,:]
#     coords, intensities = intensity_along_polyline(image, ccarrf)
#     intensities_all.append(intensities)
#     coords_all.append(coords)
    
# l=len(intensities)
# n=len(intensities)

# intensities_part = []
# coords_part=[]

# for i in range(0,170,1):
#     image = destretch[0].data[i,:,:]
#     coords, intensities = intensity_along_polyline(image, ccarrp)
#     intensities_part.append(intensities)
#     coords_part.append(coords)
    
# l2=len(intensities)
# n2=len(intensities)
# #overplot fried parameter on on this map

# #compute FFT


# dx = l / n      # Spatial sampling interval (meters)
# space_x = np.linspace(0, l, n, endpoint=False)

# all_freq = []
# all_fft = []

# for i in range(len(intensities_all)):
#     chintensities = intensities_all[i]
#     # Compute the 1D FFT
#     fft_result = np.fft.fft(chintensities/np.max(chintensities))
    
#     # Get the frequencies for the result
#     freqs = np.fft.fftfreq(n,d=dx)
    
#     all_freq.append(freqs)
    
#     all_fft.append(fft_result)
    
# all_fftarr = np.array(all_fft)
# all_freqarr = np.array(all_freq)

# pix_to_arcsec = 0.017
# arcsec_to_km=727
# km_to_Mm = 0.001

# fig,ax=plt.subplots(dpi=100,figsize=(2,5))

# ax.pcolormesh(destretch[0].data[38,:,:],cmap='gray')
# ax.plot(xs_full,ys_full,c='magenta',linestyle='solid',markersize=4,linewidth=2)
# ax.plot(xs_part,ys_part,c='darkred',linestyle='solid',markersize=4,linewidth=2)
# ax.set_xlim([1700,2400])
# ax.set_ylim([2550,1000])

# ax.pcolormesh(destretch[0].data[100,:,:],cmap='gray')
# ax.plot(xs_full,ys_full,c='magenta',linestyle='solid',markersize=4,linewidth=2)
# ax.plot(xs_part,ys_part,c='darkred',linestyle='solid',markersize=4,linewidth=2)
# ax.set_xlim([1700,2400])
# ax.set_ylim([2550,1000])

# ax.set_xticks([])
# ax.set_yticks([])

# fig.show()

    
# fig,ax=plt.subplots(3,1,dpi=100,figsize=(10,15));
# ax.flatten()[0].pcolormesh(np.arange(170),np.arange(l)*pix_to_arcsec*arcsec_to_km,np.transpose(intensities_all),cmap='afmhot');
# ax.flatten()[0].invert_yaxis()

# ax.flatten()[0].set_ylabel('Distance along trace [km]',fontsize=8)
# ax.flatten()[0].set_xlabel('Time [UT]',fontsize=8)
# ax.flatten()[1].pcolormesh(np.arange(170),arcsec_to_km*pix_to_arcsec*1/all_freqarr[0,1:n//2],np.transpose(np.power(all_psdarr[:,1:n//2],.25)),cmap='seismic');
# #ax.flatten()[1].axhspan(200,400, facecolor='grey', alpha=0.2)
# ax.flatten()[1].axhline(200,color='red')
# ax.flatten()[1].axhline(400,color='red')

# ax.flatten()[2].pcolormesh(np.arange(170),np.arange(l2)*pix_to_arcsec*arcsec_to_km,np.transpose(intensities_part),cmap='afmhot');
# ax.flatten()[2].invert_yaxis()

# ax.flatten()[2].set_ylabel('Distance along trace [km]',fontsize=8)
# ax.flatten()[2].set_xlabel('Time [UT]',fontsize=8)


# #ax.flatten()[1].invert_yaxis()
# ax.flatten()[1].axhline(50,color='grey',linestyle='dashdot') # loop size
# #ax.flatten()[1].axhline(300,color='white',linestyle='dashdot') # bead size
# ax.flatten()[1].set_ylim([25,1500])

# ax.flatten()[1].set_yscale('log')
# ax.flatten()[1].set_ylabel('Spatial scale [km]',fontsize=8)
# ax.flatten()[1].set_xlabel('Time [UT]',fontsize=8)

# ax2 = ax.flatten()[0].twinx()
# ax3 = ax.flatten()[1].twinx()

# ax2.plot(friedvbi[177:347],c='white',linewidth=2)
# ax3.plot(friedvbi[177:347],c='white',linewidth=2)

# ax4 = ax.flatten()[0].twinx()
# ax4.plot(np.arange(137),lcvbi[:137],color='black',linewidth=2);
# ax4.set_yticks([])

# ax5 = ax.flatten()[1].twinx()
# ax5.set_yticks([])


# ax2.set_ylim([0,20])
# ax3.set_ylim([0,20])

# ax2.set_ylabel(r'$r_0$ [cm]',fontsize=8)
# ax3.set_ylabel(r'$r_0$ [cm]',fontsize=8)

# ax.flatten()[0].set_xticks([0,20,40,60,80,100,120],[onlytimevbi[177],onlytimevbi[177+20],onlytimevbi[177+40],onlytimevbi[177+60],onlytimevbi[177+80],onlytimevbi[177+100],onlytimevbi[177+120]])
# ax.flatten()[1].set_xticks([0,20,40,60,80,100,120],[onlytimevbi[177],onlytimevbi[177+20],onlytimevbi[177+40],onlytimevbi[177+60],onlytimevbi[177+80],onlytimevbi[177+100],onlytimevbi[177+120]])
# ax.flatten()[2].set_xticks([0,20,40,60,80,100,120],[onlytimevbi[177],onlytimevbi[177+20],onlytimevbi[177+40],onlytimevbi[177+60],onlytimevbi[177+80],onlytimevbi[177+100],onlytimevbi[177+120]])


# ax.flatten()[0].set_xlim([0,138])
# ax.flatten()[1].set_xlim([0,138])
# ax.flatten()[2].set_xlim([0,138])

# ax.flatten()[0].tick_params(axis='both', labelsize=8)
# ax.flatten()[1].tick_params(axis='both', labelsize=8)
# ax.flatten()[2].tick_params(axis='both', labelsize=8)
# ax2.tick_params(axis='both', labelsize=8)
# ax3.tick_params(axis='both', labelsize=8)
# ax4.tick_params(axis='both', labelsize=8)
# ax5.tick_params(axis='both', labelsize=8)

# fig.show()
# # from scipy.signal import savgol_filter



# # fig,ax=plt.subplots();
# # for i in range(0,120,20):
# #     ax.plot(arcsec_to_km*pix_to_arcsec*1/all_freqarr[0,1:n//2],savgol_filter(np.abs(all_fftarr)[i,1:n//2],window_length=20,polyorder=3),label=str(i))
# # ax.set_xscale('log')
# # ax.set_yscale('log')
# # ax.set_xlim([0,2000])
# # ax.axvline(250)
# # ax.legend()
# # fig.show()































