#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 08:57:28 2026

@author: coletamburri
"""

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

dkist_coord_file = '/Users/coletamburri/Desktop/DKIST_Flares/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/ViSPcoords_newcalib.npz'
dkist_coords = np.load(dkist_coord_file)

xarr_CaII = dkist_coords['xarr_caII']
yarr_CaII = dkist_coords['yarr_caII']

xarr_hbeta = dkist_coords['xarr_hbeta']
yarr_hbeta = dkist_coords['yarr_hbeta']

visp_file = '/Users/coletamburri/Desktop/DKIST_Flares/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/ViSP_coalign_result_11Aug_Cclass'

X=np.load(visp_file)['arr_0']
Y=np.load(visp_file)['arr_1']
V=np.load(visp_file)['arr_2']

ylow= int(Y[382][0])
yhigh = int(Y[1565][0])

xlow = int(X[0][-1])
xhigh = int(X[0][0])

img4define_beads = destretch[0].data[171,:,:]
img_sel = img4define_beads[ylow:yhigh,xlow:xhigh]

npoints=8

fig,ax=plt.subplots()
ax.imshow(img_sel)

fig.show()

aa = plt.ginput(npoints,timeout=120)

lcs=[]

for i in range(0,len(aa),2):
    point1 = aa[i]
    point2 = aa[i+1]
    
    point1x = int(point1[0])
    point1y = int(point1[1])

    point2x = int(point2[0])
    point2y = int(point2[1])

    size = (point2x-point1x)*(point2y-point1y)

    lc = []

    for j in range(300):
        img_sel1 = destretch[0].data[j,ylow:yhigh,xlow:xhigh]
        lcpoint = np.nansum(img_sel1[point1y:point2y,point1x:point2x])/size
        lc.append(lcpoint/1e3)
    lcs.append(lc)

colors=['darkviolet','mediumseagreen','darkorange','crimson']
labels=['bead a','bead b','bead c','bead d']

desired_y = [-234,-238,-242,-246]
desired_ystr = ['—234','—238','—242','—246']

desired_x = [760,765]
desired_xstr = ['760','765']

# need small adjustments to account for different coordinates in destretched
# data from the 47 to 347 file (here) and the frames 177 to 347 
# (on which the OG coalignment between VBI and ViSP was done)

# Values also reflected in ribbon_slot.py
XX = 27
YY = 12

vbiy=[]
for i in range(len(desired_y)):
    indch = DKISTanalysis.find_nearest(yarr_CaII,desired_y[i])[1]
    vbiy.append(int(Y[indch][0])-ylow+YY)
    
vbix=[]
for i in range(len(desired_x)):
    indch = DKISTanalysis.find_nearest(xarr_CaII,desired_x[i])[1]
    vbix.append(int(X[0][indch])-xlow+XX)

fig,ax=plt.subplots(1,2,dpi=100,figsize=(14,4));

axes = ax.ravel()
im=axes[0].imshow(img_sel/1e6,cmap='grey')
l=0
for k in range(4):
    axes[1].plot(onlytimevbi_samp,lcs[k],label=labels[k],color=colors[k])
    maxind = np.where(lcs[k]==np.nanmax(lcs[k]))[0][0]
    #ax.axvline(onlytimevbi_samp[maxind],linestyle='--',color=colors[k])
    
    point1 = aa[l]
    point2 = aa[l+1]
    l+=2
    
    point1x = int(point1[0])
    point1y = int(point1[1])

    point2x = int(point2[0])
    point2y = int(point2[1])
    
    rect = patches.Rectangle(
        (point1x, point1y),                # (x, y) coordinates
        point2x-point1x, point2y-point1y,                   # width, height
        linewidth=2,              # border thickness
        edgecolor=colors[k],            # border color
        facecolor='none'          # 'none' makes it an open box (hollow)
    )
    
    # Add the box to the Axes
    axes[0].add_patch(rect)
        
axes[0].set_xticks(vbix,desired_xstr) 
axes[0].set_yticks(vbiy,desired_ystr)

axes[0].set_xlim([30,450])
axes[0].set_ylim([2000,1000])



axes[0].xaxis.set_minor_locator(AutoMinorLocator(6)) 
axes[0].yaxis.set_minor_locator(AutoMinorLocator(6)) 
axes[0].set_xlabel('DKIST HPC-x [arcsec]')
axes[0].set_ylabel('DKIST HPC-y [arcsec]')
fig.colorbar(im,label='Intensity [arb. units]')
axes[1].legend()
axes[1].set_xlim([0,260])
axes[1].set_xticks([0,80,160,240],[onlytimevbi_samp[0],onlytimevbi_samp[80],onlytimevbi_samp[160],onlytimevbi_samp[240]])
axes[1].xaxis.set_minor_locator(AutoMinorLocator(6)) 
axes[1].yaxis.set_minor_locator(AutoMinorLocator(6)) 
axes[1].set_xlabel('Time [UT]')
axes[1].set_ylabel(r'Integrated intensity [arb. units pix$^{-2}$]')
axes[1].grid()

slices=[20,23,27,32,35,38,41,45,50,60,105]

for i in [0,3,6,10]:
    axes[1].axvline(130+slices[i],linestyle='--',color='grey',linewidth=1)
    
fig.show()










