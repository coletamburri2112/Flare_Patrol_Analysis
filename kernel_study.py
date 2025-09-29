#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 13:42:35 2025

@author: coletamburri
"""


import dkistpkg_ct as DKISTanalysis
import numpy as np
import matplotlib.pyplot as plt
import tol_colors
import os
from astropy.io import fits
from datetime import time
import matplotlib.patches as patches

def rebin_image(arr, new_shape):
    """
    Rebins a 2D array (image) to a new shape by averaging pixel values.

    Parameters:
    arr (numpy.ndarray): The original 2D image array.
    new_shape (tuple): The desired new shape (rows, columns).
                       Each dimension of new_shape must be a factor of
                       the corresponding dimension in arr.shape.

    Returns:
    numpy.ndarray: The binned 2D array.
    """
    if not (arr.shape[0] % new_shape[0] == 0 and arr.shape[1] % new_shape[1] == 0):
        raise ValueError("New shape dimensions must be factors of original shape dimensions.")

    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

muted = DKISTanalysis.color_muted2()

l=0
CUT=2.5 # cutoff for mask-making
binning = 0 # to bin or not to bin?
diff = 0 # to diff or not to diff
if binning == 1:
    bin_x = 7 # this bins the VBI pixels to largest ViSP spatial scale (in scan dir), plus 1 for "safety"
    bin_y = 7
else:
    bin_x = 1 #initialize this parameter - will change if binning selected, below
    bin_y = 1 #initialize this parameter - will change if binning selected, below

# limits for flare region (defined by raw image)
xlow = 1600
xhigh = 2300
ylow = 900
yhigh = 2700

caII_low = 570
caII_high = 730
hep_low = 730
hep_high = 900

#define start and end times for series
starttime = 200
endtime = 450

# use destretched (0) or pre-destretched (1)?
stretched = 0



#VBI directory
path_vbi = '/Volumes/VBI_External/pid_2_11/'
folder1_vbi = 'AKDKX'
#folder2_vbi = 'BYMOL'
dir_list2_vbi = DKISTanalysis.pathdef(path_vbi,folder1_vbi)

#destretched dataset
filename = ['/Users/coletamburri/Desktop/DKIST_Code/VBI_Destretching/AKDKX/'+
                'postdestretch_dataCubeFlareImpulsivePhase.fits'][0]

visp_file = '/Users/coletamburri/Desktop/co_align_res_11Aug.npz'

vispload = np.load(visp_file)



vispX = vispload['arr_0']
vispY = vispload['arr_1']
vispavg = vispload['arr_2']

file = '/Users/coletamburri/Desktop/Misc_DKIST/11August2024_Cclass_imp_CaII.npz'

#load visp spec
vispspec = np.load(file) #arrays allspc,flare,wl,range
timesvisp = vispspec['time']
specvisp0 = vispspec['flare']
wlvisp = vispspec['wl']

#specific to this save
start = 57
nstep = 91
fullframenum = 1 #1 is the good seeing frame in visp data

# First need to define which frame in the VBI sequence 
# is most appropriate for the kernel you're searching for in the ViSP scan!
# For example, fr=21 is the time corresponding to the first kernel studied, middle of scan
# For the good seeing frame, all steps are between 22:33:19 and 22:33:24...
# these correspond roughly to indices 19 through 22 in the VBI.
fr=21

specvisp = specvisp0[start+fullframenum*nstep:start+(fullframenum+1)*nstep,:,:]
timeframe = timesvisp[start+fullframenum*nstep:start+(fullframenum+1)*nstep]


caiiavg = np.mean(specvisp[:,caII_low:caII_high,:],1)
hepavg = np.mean(specvisp[:,hep_low:hep_high,:],1)



timesvbi=[]
for i in range(len(dir_list2_vbi)):
    timesvbi.append(time(int(dir_list2_vbi[i][15:17]),
                      int(dir_list2_vbi[i][18:20]),
                      int(dir_list2_vbi[i][21:23])))
    
timesvbi = timesvbi[200:450]

#another time array
string = timesvbi[0].strftime("%H:%M:%S")
strtime =[]
for i in range(len(timesvbi)):
    strtime.append(timesvbi[i].strftime("%H:%M:%S"))

#optional - processing just from the files on HD (raw, not destretched)
if stretched == 1:
    vbi_X, vbi_Y, hdul1_vbi, dat0_vbi = DKISTanalysis.vbi_process(path_vbi,
                                                                  folder1_vbi)
    #just pixels
    vbix0 = np.arange(4096)
    vbiy0 = np.arange(4096)
    vbiX0,vbiY0= np.meshgrid(vbix0,vbiy0)
else:
    vbi_DS = fits.open(filename)

vbi_DSimgs = vbi_DS[0].data

arr = vbi_DSimgs

if binning == 1:
    new_height = vbi_DSimgs.shape[1] // bin_y
    new_width = vbi_DSimgs.shape[2] // bin_x
    
    binned = np.zeros((250,new_height,new_width))
    
    for i in range(250):
        binned[i,:,:] = rebin_image(vbi_DSimgs[i,:,:], (new_height,new_width))
    
    arr = binned

# create cumulative mask
mask = np.zeros(np.shape(arr))
timing = np.zeros(np.shape(arr)[1:3])

for i in range(100):
    l+=1
    if i>0:
        mask[i,:,:]=mask[i-1,:,:]
    maskvals = np.nonzero((arr[i,:,:]>CUT*np.median(arr[i,:,:])))
    for j in range(np.shape(maskvals)[1]):
        mask[i,maskvals[0][j],maskvals[1][j]] += 1
        #logic for timing array
        if timing[maskvals[0][j],maskvals[1][j]]==0:
            timing[maskvals[0][j],maskvals[1][j]]=l

#create timing series
for i in range(np.shape(timing)[0]):
    for j in range(np.shape(timing)[1]):
        if timing[i,j]==0:
            timing[i,j]='NaN'
            
#indices for light curves are in last cumulative mask
i_vals = []
j_vals = []

final_mask = mask[99,:,:]

for i in range(np.shape(final_mask)[0]):
    for j in range(np.shape(final_mask)[1]):
        if final_mask[i,j] > 0.0 and j>xlow/bin_x and j<xhigh/bin_x and i>ylow/bin_x and i<yhigh/bin_x:
            i_vals.append(i)
            j_vals.append(j)
            
# make light curve
lc=[]

for i in range(250):
    lc.append(np.sum(arr[i,int(ylow/bin_y):int(yhigh/bin_y),int(xlow/bin_x):int(xhigh/bin_x)]))

# cumulative mask movie
cumul_mask = np.zeros([100,np.shape(timing)[0],np.shape(timing)[1]])

for i in range(1,99):
    inds = np.where(timing==int(i))
    arrsamp = np.zeros(np.shape(timing))
    for j in range(np.shape(inds)[1]):
        arrsamp[inds[0][j],inds[1][j]]=timing[inds[0][j],inds[1][j]]
    cumul_mask[i,:,:]=cumul_mask[i-1,:,:]+arrsamp
    
cumul_mask[cumul_mask == 0] = np.nan

inst_mask = np.zeros([100,np.shape(timing)[0],np.shape(timing)[1]])

for i in range(1,99):
    inds = np.where(timing==int(i))
    arrsamp = np.zeros(np.shape(timing))
    for j in range(np.shape(inds)[1]):
        arrsamp[inds[0][j],inds[1][j]]=1.0
    inst_mask[i,:,:]=arrsamp

inst_mask[inst_mask == 0] = np.nan
if diff == 1:
    diffarr = np.zeros(np.shape(arr))
    
    m=0
    folder = '/Users/coletamburri/Desktop/diffimg_pre/'
    if ~os.path.isdir(folder):
        os.mkdir(folder)
    
    for i in range(np.shape(diffarr)[0]-20):
        diffarr[m,:,:] = np.subtract(arr[i,:,:],arr[i+20,:,:])
        
        
        # fig,ax=plt.subplots(dpi=200);
        # ax.imshow(diffarr[m,:,:],cmap='grey')
        # fig.savefig(folder+str(i)+'.png')
        
        m+=1
        


#select the thing you want to look at
fig,[ax,ax1]=plt.subplots(1,2);
ax.imshow(inst_mask[fr,ylow:yhigh,xlow:xhigh])
ax1.imshow(cumul_mask[fr,ylow:yhigh,xlow:xhigh])
fig.show()
cc=plt.ginput(2,timeout=120)

upperleftx = int(cc[0][0])
upperlefty = int(cc[0][1])
lowerrightx = int(cc[1][0])
lowerrighty = int(cc[1][1])

lcsmall=[]
for i in range(100):
    lcsmall.append(np.nansum(arr[i,ylow+upperlefty:ylow+lowerrighty,xlow+upperleftx:xlow+lowerrightx])) #limits for small thing
    
fig,ax=plt.subplots();ax.plot(lcsmall);fig.show()



fig,[(ax,ax1,ax2),(ax3,ax4,ax5)]=plt.subplots(2,3,dpi=200)
ax.imshow(arr[fr,ylow+upperlefty:ylow+lowerrighty,xlow+upperleftx:xlow+lowerrightx],cmap='hot')
ax1.imshow(inst_mask[fr,ylow+upperlefty:ylow+lowerrighty,xlow+upperleftx:xlow+lowerrightx],cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'), vmin=0, vmax=40)
ax2.imshow(cumul_mask[fr,ylow+upperlefty:ylow+lowerrighty,xlow+upperleftx:xlow+lowerrightx],cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'), vmin=0, vmax=40)
ax3.imshow(arr[fr,ylow:yhigh,xlow:xhigh],cmap='hot')
ax4.imshow(inst_mask[fr,ylow:yhigh,xlow:xhigh],cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'), vmin=0, vmax=40)
ax5.imshow(cumul_mask[fr,ylow:yhigh,xlow:xhigh],cmap=tol_colors.tol_cmap(colormap='rainbow_PuRd'), vmin=0, vmax=40)
rect = patches.Rectangle((upperleftx, upperlefty), lowerrightx-upperleftx, lowerrighty-upperlefty, linewidth=1, edgecolor='k', facecolor='none')

# Add the patch to the Axes
ax3.add_patch(rect)
rect = patches.Rectangle((upperleftx, upperlefty), lowerrightx-upperleftx, lowerrighty-upperlefty, linewidth=1, edgecolor='k', facecolor='none')
ax4.add_patch(rect)
rect = patches.Rectangle((upperleftx, upperlefty), lowerrightx-upperleftx, lowerrighty-upperlefty, linewidth=1, edgecolor='k', facecolor='none')
ax5.add_patch(rect)
ins = ax.inset_axes([0.5,0.7,0.4,0.2])
ins.plot(lcsmall,c='black')
ins.axvline(fr,c='red')
ins.set_xticks([])
ins.set_yticks([])
fig.show()


#find points
fig,[ax0,ax1]=plt.subplots(1,2)
ax0.imshow(arr[fr,ylow+upperlefty:ylow+lowerrighty,xlow+upperleftx:xlow+lowerrightx],cmap='hot')
ax1.pcolormesh(vispX,vispY,np.transpose(vispavg),cmap='hot',vmin=0.1,vmax=1)
ax1.set_xlim([xlow+upperleftx,xlow+lowerrightx]);
ax1.set_ylim([ylow+upperlefty,ylow+lowerrighty]);
ax1.invert_yaxis();
plt.show()

aa = plt.ginput(6,timeout=120)

def find_nearest_numpy(array, value):
    """
    Finds the element in a NumPy array closest to the given value.
    """
    idx = (np.abs(array - value)).argmin()
    return idx
n_points = len(aa)
colors = plt.cm.jet(np.linspace(0,1,n_points))

fig,[ax0,ax1]=plt.subplots(1,2)
ax0.imshow(arr[fr,ylow+upperlefty:ylow+lowerrighty,xlow+upperleftx:xlow+lowerrightx],cmap='hot')
ax1.pcolormesh(vispX,vispY,np.transpose(vispavg),cmap='hot',vmin=0.1,vmax=1.2)
ax1.set_xlim([xlow+upperleftx,xlow+lowerrightx]);
ax1.set_ylim([ylow+upperlefty,ylow+lowerrighty]);
ax1.invert_yaxis();
for i in range(len(aa)):
    xsel,ysel = aa[i][0],aa[i][1]
    ax1.plot(xsel,ysel,'x',color=colors[i])
plt.show()

vispx_1 = vispX[0]
vispy_1 = vispY[:,0]

xsspec = []
ysspec = []

for i in range(len(aa)):
    xsel,ysel = aa[i][0],aa[i][1]
    idx_x = find_nearest_numpy(vispx_1,xsel)
    idx_y = find_nearest_numpy(vispy_1,ysel)
    
    xsspec.append(idx_x)
    ysspec.append(idx_y)

fig,ax=plt.subplots(2,3)
for i in range(len(aa)):
    xsel,ysel = aa[i][0],aa[i][1]
    idx_x = find_nearest_numpy(vispx_1,xsel)
    idx_y = find_nearest_numpy(vispy_1,ysel)
    ax.flatten()[i].plot(wlvisp,specvisp[idx_x,:,idx_y],color=colors[i])
    ax.flatten()[i].axvline(396.85)
    ax.flatten()[i].axvline(397.01)
    ax.flatten()[i].set_xlim([396.7,397.07])
fig.show()









