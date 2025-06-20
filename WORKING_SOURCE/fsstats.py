#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 12:58:27 2025

@author: coletamburri
"""
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import skimage
import scipy
import tol_colors as tc
import matplotlib.patches as patches
from scipy.stats import skewtest
from scipy.stats import wilcoxon

valnum = '4'

#root = '/Users/coletamburri/Desktop/August_2024_DKIST_Flares/FS_Runs/'
root = '/Users/coletamburri/Downloads/for_Cole_try2/'
#root = '/Users/coletamburri/Desktop/'
#folder0 = 'small_loop_frame0_validate_SSD/'
folder0 = 'small_loop_frame0_validate'+str(valnum)+'/'
folder1 = 'small_loop_frame1_validate'+str(valnum)+'/'
folder2 = 'small_loop_frame2_validate'+str(valnum)+'/'
folder3 = 'small_loop_frame3_validate'+str(valnum)+'/'
folder4 = 'small_loop_frame4_validate'+str(valnum)+'/'

# folder0 = 'small_loop_frame0_pre_destretch_separate_all/'
# folder1 = 'small_loop_frame1_pre_destretch_separate_all/'
# folder2 = 'small_loop_frame2_pre_destretch_separate_all/'
# folder3 = 'small_loop_frame3_pre_destretch_separate_all/'
# folder4 = 'small_loop_frame4_pre_destretch_separate_all/'

filename = 'widths_errors.npz'

sample0 = np.load(root+folder0+filename)
sample1 = np.load(root+folder1+filename)
sample2 = np.load(root+folder2+filename)
sample3 = np.load(root+folder3+filename)
sample4 = np.load(root+folder4+filename)

#if gauss1
#arr1 = widthss
#arr2 = widtherrs
#arr3 = startx
#arr4 = starty
#arr5 = endx
#arr6 = endy
#arr7 = r2s
#arr8 = amps
#arr9 = note
#arr10 = time
#arr11 = ylos
#arr12 = yhis
#arr13 = xlos
#arr14 = xhis

def min_max_normalize(data):
  """
  Normalizes a list of numbers to the range [0, 1] using min-max scaling.

  Args:
    data: A list of numbers.

  Returns:
    A new list with normalized values.
  """


  min_val = min(data)
  max_val = max(data)

  if min_val == max_val:
    return [0.0] * len(data)

  normalized_data = [(x - min_val) / (max_val - min_val) for x in data]
  return normalized_data

widths0 = sample0['arr_0']
widths1 = sample1['arr_0']
widths2 = sample2['arr_0']
widths3 = sample3['arr_0']
widths4 = sample4['arr_0']

skewtest0 = skewtest(widths0,nan_policy='omit')
skewtest1 = skewtest(widths1,nan_policy='omit')
skewtest2 = skewtest(widths2,nan_policy='omit')
skewtest3 = skewtest(widths3,nan_policy='omit')
skewtest4 = skewtest(widths4,nan_policy='omit')



amps0 = sample0['arr_7']
amps1 = sample1['arr_7']
amps2 = sample2['arr_7']
amps3 = sample3['arr_7']
amps4 = sample4['arr_7']

ylo0 = sample0['arr_9']
yhi0 = sample0['arr_10']
xlo0 = sample0['arr_11']
xhi0 = sample0['arr_12']

ylo1 = sample1['arr_9']
yhi1 = sample1['arr_10']
xlo1 = sample1['arr_11']
xhi1 = sample1['arr_12']

ylo2 = sample2['arr_9']
yhi2 = sample2['arr_10']
xlo2 = sample2['arr_11']
xhi2 = sample2['arr_12']

ylo3 = sample3['arr_9']
yhi3 = sample3['arr_10']
xlo3 = sample3['arr_11']
xhi3 = sample3['arr_12']

ylo4 = sample4['arr_9']
yhi4 = sample4['arr_10']
xlo4 = sample4['arr_11']
xhi4 = sample4['arr_12']

linelox2 = sample2['arr_2']
lineloy2 = sample2['arr_3']
linehix2 = sample2['arr_4']
linehiy2 = sample2['arr_5']


widtherrs0 = sample0['arr_1']
widtherrs1 = sample1['arr_1']
widtherrs2 = sample2['arr_1']
widtherrs3 = sample3['arr_1']
widtherrs4 = sample4['arr_1']

dkistresolution = 0.016 *727
idx0 = np.where(np.logical_or(widths0<dkistresolution, widtherrs0>15))
idx1 = np.where(np.logical_or(widths1<dkistresolution, widtherrs1>15))
idx2 = np.where(np.logical_or(widths2<dkistresolution, widtherrs2>15))
idx3 = np.where(np.logical_or(widths3<dkistresolution, widtherrs3>15))
idx4 = np.where(np.logical_or(widths4<dkistresolution, widtherrs4>15))

widths0[idx0] = np.nan
amps0[idx0] = np.nan
widtherrs0[idx0] = np.nan

widths1[idx1] = np.nan
amps1[idx1] = np.nan
widtherrs1[idx1] = np.nan


widths2[idx2] = np.nan
amps2[idx2] = np.nan
widtherrs2[idx2] = np.nan


widths3[idx3] = np.nan
amps3[idx3] = np.nan
widtherrs3[idx3] = np.nan


widths4[idx4] = np.nan
amps4[idx4] = np.nan
widtherrs4[idx4] = np.nan


muted = tc.tol_cset('muted')

#os.mkdir('/Users/coletamburri/Desktop/fine_stats/')

fig,ax=plt.subplots(dpi=200)

# ax.errorbar(range(len(widths0)),widths0,widtherrs0,linestyle='',fmt='.',
#              ecolor=muted.indigo,elinewidth=2,capsize=2,zorder=0,label='Frame 0')
# ax.scatter(range(len(widths0)),widths0,c=min_max_normalize(amps0),marker='+',edgecolors='black',\
#            zorder=1,cmap='Greys')

# ax.errorbar(range(len(widths1)),widths1,widtherrs1,linestyle='',fmt='.',
#              ecolor=muted.indigo,elinewidth=2,capsize=2,zorder=0,label='Frame 1')
# ax.scatter(range(len(widths1)),widths1,c=min_max_normalize(amps1),edgecolors='black',\
#            zorder=1,cmap='Blues')    
    
# ax.errorbar(range(len(widths2)),widths2,widtherrs2,linestyle='',fmt='.',
#              ecolor=muted.rose,elinewidth=2,capsize=2,zorder=0,label='Frame 2')
# ax.scatter(range(len(widths2)),widths2,c=min_max_normalize(amps2),marker='s',edgecolors='black',\
#            zorder=1,cmap='Reds')
    
# ax.errorbar(range(len(widths3)),widths3,widtherrs3,linestyle='',fmt='.',
#              ecolor=muted.indigo,elinewidth=2,capsize=2,zorder=0,label='Frame 3')
# ax.scatter(range(len(widths3)),widths3,c=min_max_normalize(amps3),marker='*',edgecolors='black',\
#            zorder=1,cmap='Purples')
    
# ax.errorbar(range(len(widths4)),widths4,widtherrs4,linestyle='',fmt='.',
#              ecolor=muted.indigo,elinewidth=2,capsize=2,zorder=0,label='Frame 4')
# ax.scatter(range(len(widths4)),widths4,c=min_max_normalize(amps4),marker='h',edgecolors='black',\
#            zorder=1,cmap='YlOrRd')

ax.errorbar(range(len(widths0)),widths0,widtherrs0,linestyle='',fmt='.',
              ecolor=muted.indigo,elinewidth=2,capsize=2,zorder=0)
ax.scatter(range(len(widths0)),widths0,marker='H',edgecolors='black',\
            zorder=1,c=muted.indigo,label='Frame 0')

ax.errorbar(range(len(widths1)),widths1,widtherrs1,linestyle='',fmt='.',
              ecolor=muted.rose,elinewidth=2,capsize=2,zorder=0)
ax.scatter(range(len(widths1)),widths1,edgecolors='black',\
            zorder=1,c=muted.rose,label='Frame 1')    
    
ax.errorbar(range(len(widths2)),widths2,widtherrs2,linestyle='',fmt='.',
              ecolor=muted.sand,elinewidth=2,capsize=2,zorder=0)
ax.scatter(range(len(widths2)),widths2,marker='s',edgecolors='black',\
            zorder=1,c=muted.sand,label='Frame 2')
    
ax.errorbar(range(len(widths3)),widths3,widtherrs3,linestyle='',fmt='.',
              ecolor=muted.green,elinewidth=2,capsize=2,zorder=0)
ax.scatter(range(len(widths3)),widths3,marker='*',edgecolors='black',\
            zorder=1,c=muted.green,label='Frame 3')
    
ax.errorbar(range(len(widths4)),widths4,widtherrs4,linestyle='',fmt='.',
              ecolor=muted.teal,elinewidth=2,capsize=2,zorder=0)
ax.scatter(range(len(widths4)),widths4,marker='h',edgecolors='black',\
            zorder=1,c=muted.teal,label='Frame 4')
    
# ax.axhline(np.nanmedian(widths0),linestyle='--',c=muted.indigo,linewidth=4,label='Frame 0 median')
# ax.axhline(np.nanmedian(widths1),linestyle='--',c=muted.rose,linewidth=4,label = 'Frame 1 median') 
# ax.axhline(np.nanmedian(widths2),linestyle='--',c=muted.sand,linewidth=4,label='Frame 2 median')
# ax.axhline(np.nanmedian(widths3),linestyle='--',c=muted.green,linewidth=4,label = 'Frame 3 median') 
# ax.axhline(np.nanmedian(widths4),linestyle='--',c=muted.teal,linewidth=4,label='Frame 4 median')

ax.axhline(np.nanmedian(widths0),linestyle='--',c=muted.indigo,linewidth=3)
ax.axhline(np.nanmedian(widths1),linestyle='--',c=muted.rose,linewidth=3) 
ax.axhline(np.nanmedian(widths2),linestyle='--',c=muted.sand,linewidth=3)
ax.axhline(np.nanmedian(widths3),linestyle='--',c=muted.green,linewidth=3) 
ax.axhline(np.nanmedian(widths4),linestyle='--',c=muted.teal,linewidth=3)



ax.set_ylabel('Width [km]',fontsize=12,font='Tahoma')
ax.set_xticks([])

ax.axvspan(0, 15, alpha=0.5, color='black',zorder=0)
ax.axvspan(30, 45, alpha=0.5, color='black',zorder=0)

ax.set_ylim([0,100])

ax.legend(bbox_to_anchor=(.6,.65))
fig.show()

fig.savefig('/Users/coletamburri/Desktop/scatter_verify'+str(valnum)+'.png')

bins=np.arange(20,95,7.5)

fig,ax=plt.subplots(dpi=200)
ax.hist(widths0,edgecolor='k',alpha=0.3,bins=bins,zorder=1)
ax.hist(widths1,edgecolor='k',alpha=0.3,bins=bins,zorder=1)
ax.hist(widths2,edgecolor='k',alpha=0.3,bins=bins,zorder=0)
ax.hist(widths3,edgecolor='k',alpha=0.3,bins=bins,zorder=1)
ax.hist(widths4,edgecolor='k',alpha=0.3,bins=bins,zorder=0)

ax.axvline(np.nanmedian(widths0),linestyle='--',c=muted.indigo,linewidth=3)
ax.axvline(np.nanmedian(widths1),linestyle='--',c=muted.rose,linewidth=3)
ax.axvline(np.nanmedian(widths2),linestyle='--',c=muted.sand,linewidth=3) 
ax.axvline(np.nanmedian(widths3),linestyle='--',c=muted.green,linewidth=3)
ax.axvline(np.nanmedian(widths4),linestyle='--',c=muted.teal,linewidth=3) 

ax.axvline(np.nanmean(widths0),linestyle='--',c=muted.indigo,linewidth=2)
ax.axvline(np.nanmean(widths1),linestyle='--',c=muted.rose,linewidth=2)
ax.axvline(np.nanmean(widths2),linestyle='--',c=muted.sand,linewidth=2) 
ax.axvline(np.nanmean(widths3),linestyle='--',c=muted.green,linewidth=2)
ax.axvline(np.nanmean(widths4),linestyle='--',c=muted.teal,linewidth=2) 

# jing avg/med - bbso
ax.axvline(124,linestyle='--',c='black',linewidth=1,label='Jing+2016 (Mean)')

# scullion avg/med - crisp
ax.axvline(130,linestyle='-.',c='black',linewidth=1,label='Scullion+2014 (Mode)')

# brooks average/md? - iris? Other structures
ax.axvline(133,linestyle=':',c='black',linewidth=1,label='Brooks+2018 (UFS, Mean)')

# other crisp avg/med post-2014?
# other bbso avg/med post-2016?

# crisp pixel resolution
crispres = 0.059*727
ax.axvline(crispres,linestyle='--',c=muted.olive,linewidth=1,label='SST/CRISP Res.')

# iris pixel resolution
irisres = 0.166*727

ax.axvline(irisres,linestyle='--',c=muted.wine,linewidth=1,label='IRIS Res.')

# # aia pixel resolution
# aiares = 1.5*727
# ax.axvline(aiares,linestyle='--',c=muted.purple,linewidth=1,label='SDO/AIA')
# dkist halpha pixel resolution
dkisthalphares = 0.017*727
ax.axvline(dkisthalphares,linestyle='--',c='black',linewidth=1,label=r'DKIST/VBI H$\alpha$ Res.')

# bbso/nst pixel resolution
nst_visres = 0.03*727 # according to jing+2016
ax.axvline(nst_visres,linestyle='--',c=muted.purple,linewidth=1,label=r'BBSO/VIS H$\alpha$ Res.')
ax.set_ylabel('Num. Occurrences',font='Tahoma',fontsize=12)




ax.set_xlabel('Width [km]',fontsize=12,font='Tahoma')

ax.set_xlabel('Width [km]',fontsize=12,font='Tahoma')
#ax.legend(fontsize=7,bbox_to_anchor=(.85,.45))
fig.show()

fig.savefig('/Users/coletamburri/Desktop/histo_verify'+str(valnum)+'.png')

widthsall = np.concatenate((widths0,widths1,widths2,widths3,widths4))

fig,ax=plt.subplots(dpi=200)
ax.hist(widthsall,edgecolor='k',alpha=0.3,bins=bins,zorder=1,color='black')
ax.axvline(np.nanmedian(widthsall),linestyle='--',c=muted.indigo,linewidth=3,label='All widths median')

# jing avg/med - bbso
ax.axvline(124,linestyle='--',c='black',linewidth=1,label='Jing+2016 (Mean)')

# scullion avg/med - crisp
ax.axvline(130,linestyle='-.',c='black',linewidth=1,label='Scullion+2014 (Mode)')

# brooks average/md? - iris? Other structures
ax.axvline(133,linestyle=':',c='black',linewidth=1,label='Brooks+2018 (UFS, Mean)')

# other crisp avg/med post-2014?
# other bbso avg/med post-2016?

# crisp pixel resolution
crispres = 0.059*727
ax.axvline(crispres,linestyle='--',c=muted.olive,linewidth=1,label='SST/CRISP Res.')

# iris pixel resolution
irisres = 0.166*727

ax.axvline(irisres,linestyle='--',c=muted.wine,linewidth=1,label='IRIS Res.')

# # aia pixel resolution
# aiares = 1.5*727
# ax.axvline(aiares,linestyle='--',c=muted.purple,linewidth=1,label='SDO/AIA')
# dkist halpha pixel resolution
dkisthalphares = 0.017*727
ax.axvline(dkisthalphares,linestyle='--',c='black',linewidth=1,label=r'DKIST/VBI H$\alpha$ Res.')

# bbso/nst pixel resolution
nst_visres = 0.03*727 # according to jing+2016
ax.axvline(nst_visres,linestyle='--',c=muted.purple,linewidth=1,label=r'BBSO/VIS H$\alpha$ Res.')
ax.set_ylabel('Num. Occurrences',font='Tahoma',fontsize=12)




ax.set_xlabel('Width [km]',fontsize=12,font='Tahoma')

ax.set_xlabel('Width [km]',fontsize=12,font='Tahoma')
ax.legend(fontsize=7,bbox_to_anchor=(.85,.45))
fig.show()

fig.savefig('/Users/coletamburri/Desktop/allwidths_verify'+str(valnum)+'alt.png')

print(np.nanmean(widthsall))


# methods for symmetry
#wilcox test
statistic, pvalue = wilcoxon(widthsall - np.nanmedian(widthsall),nan_policy='omit')
print(f"Wilcoxon statistic: {statistic}, p-value: {pvalue}")

alpha = 0.05
if pvalue > alpha:
    print("The distribution is likely symmetric")
else:
    print("The distribution is likely not symmetric")
    
#2sample ks test

#pearson's 1st:
p1 = 3*(np.nanmean(widthsall)-np.nanmedian(widthsall))/np.nanstd(widthsall)
print(p1)

#pearsons's 2nd:
p2 = (np.nanmean(widthsall)-scipy.stats.mode(widthsall,nan_policy='omit')[0])/np.nanstd(widthsall)
print(p2)


#show cuts
# for npz loading
path = '/Users/coletamburri/Desktop/VBI_Destretching/'
folder_vbi = 'AXXJL/' # 8 August X-class flare decay phase
filename = 'AXXJLselection_predestretch.npz'
array = np.load(path+folder_vbi+filename)['first50']

#frame to work with
frame = array[2,:,:]

# Constants
spatial_samp = 0.017 # for vbi red at 656nm
arcsec_to_km = 727 # approximate arcsec to km conversion

# X and Y coordinates of frame
xarr = np.arange(np.shape(frame)[0])
yarr = np.arange(np.shape(frame)[1])

# X and Y coordinates, in KM
xarr_km = xarr*spatial_samp
yarr_km = yarr*spatial_samp

# Meshgrid for plotting
XKM,YKM =np.meshgrid(xarr_km,yarr_km)

normalized = (frame-frame.min()) /(frame.max() -frame.min())
#gamma = 0.7  # Adjust this value to control brightness (gamma < 1 brightens)
#corrected_data = np.power(normalized, gamma)

import matplotlib

# Plot first frame
fig,ax=plt.subplots(dpi=200,figsize=(10,10))
ax.imshow(np.log10(normalized),cmap=matplotlib.colormaps['afmhot'],vmin=np.log10(.2),vmax=np.log10(0.92))
ax.set_aspect('equal')
for i in range(len(ylo1)):
    rect = patches.Rectangle((ylo0[i], xlo0[i]), yhi0[i]-ylo0[i], xhi0[i]-xlo0[i], linewidth=1, edgecolor='r', \
                              facecolor='none')
    ax.add_patch(rect)
ax.set_xticks([])
ax.set_yticks([])
plt.show()

fig.savefig('/Users/coletamburri/Desktop/frame2_context')

#frame to work with
frame = array[1,:,:]

# Constants
spatial_samp = 0.017 # for vbi red at 656nm
arcsec_to_km = 727 # approximate arcsec to km conversion

# X and Y coordinates of frame
xarr = np.arange(np.shape(frame)[0])
yarr = np.arange(np.shape(frame)[1])

# X and Y coordinates, in KM
xarr_km = xarr*spatial_samp
yarr_km = yarr*spatial_samp

# Meshgrid for plotting
XKM,YKM =np.meshgrid(xarr_km,yarr_km)

# Plot first frame
fig,ax=plt.subplots(dpi=200,figsize=(10,10))
ax.imshow(frame,cmap='grey')
ax.set_aspect('equal')
for i in range(len(ylo1)):
    rect = patches.Rectangle((ylo1[i], xlo1[i]), yhi1[i]-ylo1[i], xhi1[i]-xlo1[i], linewidth=1, edgecolor='r', \
                              facecolor='none')
    ax.add_patch(rect)
ax.set_xticks([])
ax.set_yticks([])
plt.show()



# #fig.savefig('/Users/coletamburri/Desktop/fine_stats/frame1_context')

# plotting the frames to select from

fig,ax=plt.subplots(dpi=300);
ax.imshow(np.log10(normalized[xlo0[2]:xhi0[2],ylo0[2]:yhi0[2]]),cmap=matplotlib.colormaps['afmhot'],vmin=np.log10(.2),vmax=np.log10(0.92))

for i in range(15):
    ax.plot([linelox2[i+30], linehix2[i+30]], [lineloy2[i+30], linehiy2[i+30]],c='red')
plt.show()

#plotting the fits









