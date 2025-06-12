#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 12:42:33 2025

@author: coletamburri
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

from sunpy.net import Fido, attrs as a
import dkist.net
from datetime import date, datetime, timedelta, time
import matplotlib.dates as mdates
import sunpy.visualization.colormaps as cm
import matplotlib
import tol_colors as tc

def perdelta(start, end, delta):
    curr = start
    while curr <= end:
        yield curr
        curr += delta
        
def normalize_3d_array(arr):
    min_val = np.nanmin(arr)
    max_val = np.nanmax(arr)

    return (arr - min_val) / (max_val - min_val)

#make times
sttime = datetime(2024,8,8,20,12,32,333333)
endtime = datetime(2024,8,8,21,5,7,0)


stack=[]
for result in perdelta(sttime , endtime, timedelta(seconds=2.666666)):
    stack.append(str(result))
    
timeshhmmss = []

# timesdt = []

# for i in stack:
#     timesdt.append(datetime.strptime(i, '%Y-%m-%dT%H:%M:%S.%f'))
    
# timesdthr = []

# for i in timesdt:
#     timesdthr.append(datetime.strftime(i, 
#                                  "%H:%M:%S"))

# times in correct format for plotting
for i in range(len(stack)):
    timeshhmmss.append(stack[i][-15:-7])
    
path = '/Volumes/VBI_External/pid_2_11/'

folder_vbi = 'AXXJL'
#'AWYMX'
dir_list = os.listdir(path+folder_vbi)
dir_list.sort()
dir_list2 = []

#for i in range(len(dir_list)):
for i in range(1182):
    filename = dir_list[i]
    if filename[-5:] == '.fits' and '_I_' in filename:
        dir_list2.append(filename)
    
# load fits file
times = []
xarrs = []
yarrs = []

for i in range(50):
    hdul = fits.open(path+folder_vbi+'/'+dir_list2[i+150])
    times.append(hdul[1].header['DATE-BEG'])
    xloc = hdul[1].header['CRVAL1']
    yloc = hdul[1].header['CRVAL2']
    xpix = hdul[1].header['CRPIX1']
    ypix = hdul[1].header['CRPIX2']
    xdelt = hdul[1].header['CDELT1']
    ydelt = hdul[1].header['CDELT2']
    
    xarr = np.zeros(4095)
    yarr = np.zeros(4095)
    
    xarr[int(xpix)]=xloc
    yarr[int(ypix)]=yloc
    
    for j in range(int(xpix)+1,4095,1):
        xarr[j] = xarr[j-1] + xdelt
    for j in range(int(ypix)+1,4095,1):
        yarr[j] = yarr[j-1] + ydelt
        
    for j in range(int(xpix),-1,-1):
        xarr[j] = xarr[j+1] - xdelt
    for j in range(int(ypix),-1,-1):
        yarr[j] = yarr[j+1] - ydelt
        
    xarrs.append(xarr)
    yarrs.append(yarr)
    
    
    
    
#this file is indices 150 to 250
data1 = np.load('/Users/coletamburri/Desktop/VBI_Destretching/AXXJL/AXXJLselections.npz')
#data = normalize_3d_array(data1['first50'])

data = fits.open('/Volumes/VBI_External/postdestretch_dataCubeX_class_decay_full.fits')[0].data
#data = fits.open('/Users/coletamburri/Desktop/VBI_Destretching/AWYMX/postdestretch_histomatch_dataCubeX_class_decay_blue_continuum.fits')
#data = fits.open('/Users/coletamburri/Desktop/VBI_Destretching/AKDKX/postdestretch_histomatch_dataCubeFlareImpulsivePhase.fits')

props = dict(edgecolor='black',facecolor='white', alpha=0.8,boxstyle='square,pad=0.4')
data1179 = fits.open('/Volumes/VBI_External/pid_2_11/AXXJL/VBI_2024_08_08T21_04_56_333333_00656282_I_AXXJL_L1.fits')[1].data

data750 = fits.open('/Volumes/VBI_External/pid_2_11/AXXJL/VBI_2024_08_08T20_45_52_333333_00656282_I_AXXJL_L1.fits')[1].data
data260 = fits.open('/Volumes/VBI_External/pid_2_11/AXXJL/VBI_2024_08_08T20_24_05_666666_00656282_I_AXXJL_L1.fits')[1].data
# fig,ax=plt.subplots(1,1,dpi=300,figsize=(10,5))
# for i in range(1):
#     #X,Y=np.meshgrid(xarrs[i],np.transpose(yarrs[i]))
#     ax.flatten()[i].pcolormesh(data[0].data[i*6,:,:],cmap='grey')
#     ax.flatten()[i].set_xticklabels([])
#     ax.flatten()[i].set_yticklabels([])
#     ax.flatten()[i].invert_yaxis()
#     # place a text box in upper left in axes coords
#     #ax.flatten()[i].text(375, -150, str(timeshhmmss[i*6])+' UT',fontsize=5,
#     #        verticalalignment='top', bbox=props)
#     ax.flatten()[i].set_aspect('equal')
#     ax.flatten()[i].tick_params(axis='both', which='minor', labelsize=5)
    

# fig.subplots_adjust(wspace=0, hspace=0)
# fig.show()

times=[]
xarrs=[]
yarrs=[]

i =0

hdul = fits.open(path+folder_vbi+'/'+dir_list2[i])
#hdul = fits.open()
#times.append(hdul[1].header['DATE-BEG'])
xloc = hdul[1].header['CRVAL1']
yloc = hdul[1].header['CRVAL2']
xpix = hdul[1].header['CRPIX1']
ypix = hdul[1].header['CRPIX2']
xdelt = hdul[1].header['CDELT1']
ydelt = hdul[1].header['CDELT2']

xarr = np.zeros(4095)
yarr = np.zeros(4095)

#for vbi-red
xarr[int(xpix)]=xloc
yarr[int(ypix)]=yloc
# for vbi-blue
#xarr[int(xpix)] = 385
#yarr[int(ypix)] = -170

for j in range(int(xpix)+1,4095,1):
    xarr[j] = xarr[j-1] + xdelt
for j in range(int(ypix)+1,4095,1):
    yarr[j] = yarr[j-1] + ydelt
    
for j in range(int(xpix),-1,-1):
    xarr[j] = xarr[j+1] - xdelt
for j in range(int(ypix),-1,-1):
    yarr[j] = yarr[j+1] - ydelt
    

xarrs.append(xarr)
yarrs.append(yarr)

X,Y = np.meshgrid(xarrs[0],yarrs[0][::-1])
#X,Y = np.meshgrid(xarrs[0],yarrs[0][::-1])

#normalized = data[0].data[i]/data[0].data[i].max()   

#data=data260[0][:-1,:-1]
#normalized = (data-data.min()) /(data.max() -data.min())

normalized = (data[i]-data[i].min()) /(data[i].max() -data[i].min())
#gamma = 0.7  # Adjust this value to control brightness (gamma < 1 brightens)
#corrected_data = np.power(normalized, gamma)

print(stack[i])

fig,ax=plt.subplots(dpi=200,figsize=(20,20))
#im = ax.pcolormesh(X,Y,np.log10(normalized),cmap=matplotlib.colormaps['afmhot'],vmin=np.log10(.15),vmax=np.log10(0.96)) #4
ax.pcolormesh(X,Y,np.log10(normalized),cmap=matplotlib.colormaps['afmhot'],vmin=np.log10(.2),vmax=np.log10(0.92)) # 1 and 3
#im = ax.pcolormesh(X,Y,np.log10(normalized),cmap=matplotlib.colormaps['afmhot'],vmin=np.log10(.25),vmax=np.log10(0.92)) #2
#im = ax.imshow(normalized,cmap=matplotlib.colormaps['afmhot'],vmin=0.07,vmax=0.7)
ax.set_aspect('equal')
#ax.invert_yaxis()
#ax.invert_xaxis()
ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
ax.tick_params(axis='both', which='minor', labelsize=2)
fig.show()

#masking

lower_intensity = 0.4  # Example lower threshold
upper_intensity = 0.285 # Example upper threshold

mask = np.copy(normalized)
mask[mask<lower_intensity] = 0

maskneg = np.copy(normalized)
maskneg[maskneg>upper_intensity]=0

fig,ax=plt.subplots();ax.imshow(maskneg);fig.show()

lower_masked = mask[2500:4000,1200:2800]
fingers = mask[200:800,800:1800]
fingersnonmasked = normalized[200:800,800:1800]
upper_masked = mask[300:900,3000:4000]
lower_masked_small = mask[3200:3500,2000:2300]
lower_nonmasked = normalized[2500:4000,1200:2800]
upper_nonmasked = normalized[300:900,3000:4000]
arcade = normalized[700:1500,2000:3000]
qs = normalized[20:200,0:2500]
arcade_masked = maskneg[800:1800,2200:3500]

ncuts = 3 # for each power spec

# first plot for the ribbons - positive mask
fig,ax=plt.subplots(dpi=300)
ax.imshow(normalized)
fig.show()
aa = plt.ginput(ncuts*2,timeout = 120)

# next plot for the arcade - negative mask
fig,ax=plt.subplots(dpi=300)
ax.imshow(normalized)
fig.show()
bb = plt.ginput(ncuts*2,timeout = 120)

# next plot for the quiet sun - no mask
fig,ax=plt.subplots(dpi=300)
ax.imshow(normalized)
fig.show()
cc = plt.ginput(ncuts*2,timeout = 120)

import skimage
#for ribbon
n=10
fig1,ax1=plt.subplots();

for j in range(0,2*ncuts,2):
    print(j)
    x0, y0, x1, y1 = aa[j][0], aa[j][1], aa[j+1][0], aa[j+1][1]

    length = int(np.hypot(x1-x0,y1-y0))
    
    # Define the x and y in # pixels along the cut
    x, y = np.linspace(y0, y1, length), np.linspace(x0, x1, length)
    
    # Find the intensities nearest to the (x,y) coordinates above
    # This essentially finds the intensity profile along the cut.
    zi = normalized[x.astype(int), y.astype(int)]
    
    # Use skimage to find the intensity profile along the line.
    # skimage.measure.profile_line returns a 1D array of intensity values 
    # in a directly line from (x0,y0) to (x1,y1), of length equal to the 
    # ceil of the computed length of the line (in units of pixels)
    profile = skimage.measure.profile_line(normalized,[y1,x0],[y0,x1])
    
    # Convert the length of the skimage output to arcsec
    xdirection = np.arange(len(profile))
    
    # Plot intensity profile in separate window
    #fig,ax=plt.subplots(dpi=300)
    a#x.plot(profile,'-x')
    #fig.show()
    
    ps = np.abs(np.fft.fft(profile))**2

    step = 0.017
    freqs = np.fft.fftfreq(profile.size, step)
    idx = np.argsort(freqs)

    ps_smooth = np.convolve(ps, np.ones(n)/n, mode='full')
    ax1.loglog((1/freqs[idx])*727, ps_smooth[idx])
ax1.set_xlim([0,3000])
ax1.axvline(0.034*727)
ax1.axvline(45)
fig1.show()
   
fig1,ax1=plt.subplots();

for j in range(0,2*ncuts,2):
    print(j)
    x0, y0, x1, y1 = bb[j][0], bb[j][1], bb[j+1][0], bb[j+1][1]

    length = int(np.hypot(x1-x0,y1-y0))
    
    # Define the x and y in # pixels along the cut
    x, y = np.linspace(y0, y1, length), np.linspace(x0, x1, length)
    
    # Find the intensities nearest to the (x,y) coordinates above
    # This essentially finds the intensity profile along the cut.
    zi = normalized[x.astype(int), y.astype(int)]
    
    # Use skimage to find the intensity profile along the line.
    # skimage.measure.profile_line returns a 1D array of intensity values 
    # in a directly line from (x0,y0) to (x1,y1), of length equal to the 
    # ceil of the computed length of the line (in units of pixels)
    profile = skimage.measure.profile_line(normalized,[y1,x0],[y0,x1])
    
    # Convert the length of the skimage output to arcsec
    xdirection = np.arange(len(profile))
    
    # Plot intensity profile in separate window
    #fig,ax=plt.subplots(dpi=300)
    a#x.plot(profile,'-x')
    #fig.show()
    
    ps = np.abs(np.fft.fft(profile))**2

    step = 0.017
    freqs = np.fft.fftfreq(profile.size, step)
    idx = np.argsort(freqs)
    

    ps_smooth = np.convolve(ps, np.ones(n)/n, mode='full')
    ax1.loglog((1/freqs[idx])*727, ps_smooth[idx])
ax1.set_xlim([0,3000])
ax1.axvline(0.034*727)
ax1.axvline(45)
fig1.show()
        

fig1,ax1=plt.subplots();
for j in range(0,2*ncuts,2):
    print(j)
    x0, y0, x1, y1 = cc[j][0], cc[j][1], cc[j+1][0], cc[j+1][1]

    length = int(np.hypot(x1-x0,y1-y0))
    
    # Define the x and y in # pixels along the cut
    x, y = np.linspace(y0, y1, length), np.linspace(x0, x1, length)
    
    # Find the intensities nearest to the (x,y) coordinates above
    # This essentially finds the intensity profile along the cut.
    zi = normalized[x.astype(int), y.astype(int)]
    
    # Use skimage to find the intensity profile along the line.
    # skimage.measure.profile_line returns a 1D array of intensity values 
    # in a directly line from (x0,y0) to (x1,y1), of length equal to the 
    # ceil of the computed length of the line (in units of pixels)
    profile = skimage.measure.profile_line(normalized,[y1,x0],[y0,x1])
    
    # Convert the length of the skimage output to arcsec
    xdirection = np.arange(len(profile))
    
    # Plot intensity profile in separate window
    #fig,ax=plt.subplots(dpi=300)
    a#x.plot(profile,'-x')
    #fig.show()
    
    ps = np.abs(np.fft.fft(profile))**2

    step = 0.017
    freqs = np.fft.fftfreq(profile.size, step)
    idx = np.argsort(freqs)

    ps_smooth = np.convolve(ps, np.ones(n)/n, mode='full')
    ax1.loglog((1/freqs[idx])*727, ps_smooth[idx])
ax1.set_xlim([0,3000])
ax1.axvline(0.034*727)
ax1.axvline(45)
fig1.show()


# def power_spectrum_1d(image):
#     """
#     Calculates the 1D power spectrum of a 2D image.

#     Args:
#         image (numpy.ndarray): A 2D numpy array representing the image.

#     Returns:
#         tuple: A tuple containing two numpy arrays:
#             - k_radii: Radial spatial frequencies.
#             - power_spectrum_1d: 1D power spectrum.
#     """
#     # 1. 2D Fast Fourier Transform (FFT)
#     fft_2d = np.fft.fft2(image)
#     fft_2d_shifted = np.fft.fftshift(fft_2d)

#     # 2. Power Spectrum
#     power_spectrum_2d = np.abs(fft_2d_shifted)**2
    
#     power_spectrum_2d /= power_spectrum_2d.size

#     #power_spectrum_2d /= 500
#     # 3. Radial Averaging
#     ny, nx = image.shape
#     cy, cx = ny // 2, nx // 2
    
#     y_indices, x_indices = np.indices((ny, nx))
#     radii = np.sqrt((x_indices - cx)**2 + (y_indices - cy)**2).astype(int)

#     # Create bins for radii
#     max_radius = np.max(radii)
#     k_radii = np.arange(0, max_radius + 1)
#     power_spectrum_1d = np.zeros(len(k_radii))

#     # Average power for each radius
#     for i, k in enumerate(k_radii):
#         mask = radii == k
#         if np.any(mask):
#             power_spectrum_1d[i] = np.mean(power_spectrum_2d[mask])

#     return k_radii, power_spectrum_1d

# k_rad_mask_up, power_spec_mask_up = power_spectrum_1d(upper_masked)
# k_rad_mask_fingers, power_spec_mask_fingers = power_spectrum_1d(fingers)
# k_rad_mask_fingers_nonmask, power_spec_mask_fingers_nonmask = power_spectrum_1d(fingersnonmasked)

# k_rad_unmask_up, power_spec_unmask_up = power_spectrum_1d(upper_nonmasked)

# k_rad_mask_down, power_spec_mask_down = power_spectrum_1d(lower_masked)
# k_rad_mask_down_small, power_spec_mask_down_small = power_spectrum_1d(lower_masked_small)

# k_rad_unmask_down, power_spec_unmask_down = power_spectrum_1d(lower_nonmasked)

# k_rad_full, power_spec_full = power_spectrum_1d(normalized)

# k_rad_qs, power_spec_qs = power_spectrum_1d(qs)

# arcade_k_rad, arcade_power_spec = power_spectrum_1d(arcade)
# arcade_k_rad_mask, arcade_power_spec_mask = power_spectrum_1d(arcade_masked)

# muted = tc.tol_cset('muted')

# fig,ax=plt.subplots()
# ax.loglog(k_rad_mask_down, power_spec_mask_down, '-o', label='Masked,lower',color=muted[0])
# ax.loglog(k_rad_mask_fingers, power_spec_mask_fingers, '-o', label='Fingers,masked',color=muted[7])

# ax.loglog(k_rad_mask_up, power_spec_mask_up, '-o', label='High right mask',color=muted[1])
# ax.loglog(k_rad_mask_down_small, power_spec_mask_down_small, '-o', label='Masked,lower,subfield',color=muted[2])
# #ax.loglog(arcade_k_rad,arcade_power_spec,'-o',label='Arcade',color=muted[3])
# ax.loglog(arcade_k_rad_mask,arcade_power_spec_mask,'-o',label='Arcade, masked',color=muted[6])

# #ax.loglog(k_rad_unmask_up, power_spec_unmask_up, '-o', label='NonMasked,higher')
# #ax.loglog(k_rad_unmask_down, power_spec_unmask_down, '-o', label='NonMasked,lower')
# ax.loglog(k_rad_full, power_spec_full, '-o', label='Full image',color='black')
# ax.loglog(k_rad_qs, power_spec_qs, '-o', label='qs',color=muted[5])
# #ax.axvline(10,label='10km')
# #ax.axvline(100,label='100km')

# ax.set_xlim([0,150])
# ax.set_ylim([0.01,10000])
# ax.legend()
# fig.show()

# # # subtractqs spec?
# # fig,ax=plt.subplots()
# # ax.loglog(k_rad_mask_down[0:30], power_spec_mask_down[0:30]-power_spec_qs[0:30], '-o', label='Masked,lower',color=muted[0])
# # ax.loglog(k_rad_mask_up[0:30], power_spec_mask_up[0:30]-power_spec_qs[0:30], '-o', label='Masked,higher',color=muted[1])
# # ax.loglog(k_rad_mask_down_small[0:30], power_spec_mask_down_small[0:30]-power_spec_qs[0:30], '-o', label='Masked,lower,small',color=muted[2])
# # ax.loglog(arcade_k_rad[0:30],arcade_power_spec[0:30]-power_spec_qs[0:30],'-o',label='Arcade',color=muted[3])
# # #ax.loglog(k_rad_unmask_up, power_spec_unmask_up, '-o', label='NonMasked,higher')
# # #ax.loglog(k_rad_unmask_down, power_spec_unmask_down, '-o', label='NonMasked,lower')
# # ax.loglog(k_rad_full[0:30], power_spec_full[0:30]-power_spec_qs[0:30], '-o', label='Full image',color='black')
# # ax.loglog(k_rad_qs[0:30], power_spec_qs[0:30]-power_spec_qs[0:30], '-o', label='qs',color=muted[5])
# # #ax.axvline(10,label='10km')
# # #ax.axvline(100,label='100km')

# # #ax.set_xlim([0,5000])
# # #ax.set_ylim([0,1e11])
# # ax.legend()
# # fig.show()

# fig,ax=plt.subplots();
# ax.imshow(normalized)
# rect = patches.Rectangle((1200, 2500), 1600, 1500, linewidth=1, edgecolor=muted[0], facecolor='none')
# ax.add_patch(rect)
# rect = patches.Rectangle((3000, 300), 1000, 600, linewidth=1, edgecolor=muted[1], facecolor='none')
# ax.add_patch(rect)
# rect = patches.Rectangle((2000, 3200), 300, 300, linewidth=1, edgecolor=muted[2], facecolor='none')
# ax.add_patch(rect)
# rect = patches.Rectangle((2000, 700), 1000, 500, linewidth=1, edgecolor=muted[3], facecolor='none')
# ax.add_patch(rect)
# rect = patches.Rectangle((2000, 20), 500, 180, linewidth=1, edgecolor=muted[5], facecolor='none')
# ax.add_patch(rect)
# rect = patches.Rectangle((800, 200), 1000, 600, linewidth=1, edgecolor=muted[7], facecolor='none')
# ax.add_patch(rect)
# fig.show()

# # analysis of just the fingers and the arcade connecting:
    
# fig,ax=plt.subplots()
# ax.plot((1/k_rad_mask_fingers), power_spec_mask_fingers, '-o', label='Fingers,masked',color=muted[7])

# ax.plot((1/arcade_k_rad_mask),arcade_power_spec_mask,'-o',label='Arcade, masked',color=muted[6])

# ax.plot((1/k_rad_qs), power_spec_qs, '-o', label='qs',color=muted[5])
# #ax.axvline(10,label='10km')
# #ax.axvline(100,label='100km')

# ax.set_ylim([-0.2,10])
# ax.set_xlim([0.0,0.2])
# # ax.axvline(12)
# # ax.axvline(45)
# # ax.axvline(55)
# # ax.axvline(80)
# ax.legend()
# fig.show()

# fig,ax=plt.subplots()
# ax.loglog(k_rad_mask_fingers_nonmask, power_spec_mask_fingers_nonmask, '-o', label='Fingers,nonmasked',color=muted[7])

# ax.loglog(arcade_k_rad,arcade_power_spec,'-o',label='Arcade, nonmasked',color=muted[6])

# ax.loglog(k_rad_qs, power_spec_qs, '-o', label='qs',color=muted[5])
# #ax.axvline(10,label='10km')
# #ax.axvline(100,label='100km')

# ax.set_xlim([0,150])
# ax.set_ylim([0.0001,10000])
# ax.legend()
# fig.show()


# fig,ax=plt.subplots();
# ax.imshow(normalized)
# rect = patches.Rectangle((800, 200), 1000, 600, linewidth=1, edgecolor=muted[7], facecolor='none')
# ax.add_patch(rect)
# rect = patches.Rectangle((0, 20), 2500, 180, linewidth=1, edgecolor=muted[5], facecolor='none')
# ax.add_patch(rect)
# rect = patches.Rectangle((2200, 800), 1300, 1000, linewidth=1, edgecolor=muted[6], facecolor='none')
# ax.add_patch(rect)
# fig.show()



