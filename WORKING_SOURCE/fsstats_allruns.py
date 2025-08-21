#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 09:48:23 2025

@author: coletamburri
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
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
import re

matplotlib.use('Qt5Agg')


largestind = [269,
 261,
 3,
 262,
 263,
 272,
 275,
 274,
 0,
 273,
 2,
 4,
 178,
 1,
 71,
 174,
 271,
 70,
 175,
 268] #defined by largest fried parameters

largestfried = [13.01338884612327,
 12.32584740571319,
 12.186340768244683,
 12.151014480611355,
 11.958786651892,
 11.932966239867623,
 11.831252200753232,
 11.70372880498985,
 11.613607159712098,
 11.55334283180776,
 11.508651467029033,
 11.426981198485024,
 11.41579557969453,
 11.389601873499458,
 11.313371521160983,
 11.159399530029892,
 11.029652443975227,
 10.985195954875286,
 10.922861871485543,
 10.7695324996077]

dkistresolution = 0.016 *727

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

muted = tc.tol_cset('muted')

root = '/Users/coletamburri/Desktop/August_2024_DKIST_Flares/FS_Runs/'

folders = []
samples = []
widths = []
skewtests = []
amps = []
ylos = []
yhis = []
xlos = []
xhis = []
lineloxs = []
lineloys = []
linehixs = []
linehiys = []
widtherrs = []
cents= []
stds = []
slopes = []
intercepts= []
sts= []
ends= []

clrs = ['#E8ECFB', '#D9CCE3', '#D1BBD7', '#CAACCB', '#BA8DB4',
        '#AE76A3', '#AA6F9E', '#994F88', '#882E72', '#1965B0',
        '#437DBF', '#5289C7', '#6195CF', '#7BAFDE', '#4EB265',
        '#90C987', '#CAE0AB', '#F7F056', '#F7CB45', '#F6C141',
        '#F4A736', '#F1932D', '#EE8026', '#E8601C', '#E65518',
        '#DC050C', '#A5170E', '#72190E', '#42150A','#125A56', 
        '#00767B', '#238F9D', '#42A7C6', '#60BCE9','#9DCCEF', 
        '#C6DBED', '#DEE6E7', '#ECEADA', '#F0E6B2','#F9D576']
widths = []
filename2='widths_errors.npz'

for filename in os.listdir(root):
    if re.search('small_loop', filename):
        if os.path.isdir(root+filename+'/'+filename2):
            samp = np.load(root+filename+'/'+filename2)
            widths.append(np.asarray(samp['arr_0'][0]))
        
for i in range(len(largestind)):
#for i in range(len(selind)):
    #if i == 2 or i == 16 or i == 1:
    #    folders.append('/Users/coletamburri/Desktop/single_loop_frame'+str(largestind[i])+'final2/')
    #else:
    folders.append(root+'small_loop_frame'+str(largestind[i])+'final/')
    #folders.append('/Users/coletamburri/Desktop/'+'single_loop_frame'+str(selind[i])+'_testlater_blurred/')
    #folders.append(root+'small_loop_frame'+str(selind[i])+'_validate2/')
    filename='widths_errors.npz'
    samples.append(np.load(folders[i]+filename))
    widths.append(samples[i]['arr_0'])
    #skewtest.append(skewtest(widths[i],nan_policy='omit'))
    amps.append(samples[i]['arr_7'])
    ylos.append(samples[i]['arr_9'])
    yhis.append(samples[i]['arr_10'])
    xlos.append(samples[i]['arr_11'])
    xhis.append(samples[i]['arr_12'])
    lineloxs.append(samples[i]['arr_2'])
    lineloys.append(samples[i]['arr_3'])
    linehixs.append(samples[i]['arr_4'])
    linehiys.append(samples[i]['arr_5'])    
    widtherrs.append(samples[i]['arr_1'])
    if i == 2:
        cents.append(samples[i]['arr_13'])
        stds.append(samples[i]['arr_14'])
        slopes.append(samples[i]['arr_15'])
        intercepts.append(samples[i]['arr_16'])
        sts.append(samples[i]['arr_17'])   
        ends.append(samples[i]['arr_18']) 
for i in range(len(largestind)):
    idx = np.where(np.logical_or(widths[i]<dkistresolution,widtherrs[i]>13))
    widths[i][idx] = np.nan
    amps[i][idx] = np.nan
    widtherrs[i][idx] = np.nan
    
fig,ax=plt.subplots(dpi=200)

for i in range(len(largestind)):
    ax.errorbar(range(len(widths[i])),widths[i],widtherrs[i],c=clrs[2*i],linestyle='',fmt='.',elinewidth=1,capsize=1,zorder=0,markersize=3)
    ax.scatter(range(len(widths[i])),widths[i],3,c=clrs[2*i],\
                zorder=1,label='Frame '+str(largestind[i]))
    #ax.axhline(np.nanmean(widths[i]),c=clrs[i+5],linestyle='--',linewidth=2)
    
ax.set_ylabel('Width [km]',fontsize=12,font='Tahoma')
ax.set_xticks([])

ax.axvspan(0, 15, alpha=0.5, color='black',zorder=0)
ax.axvspan(30, 45, alpha=0.5, color='black',zorder=0)

ax.set_ylim([0,100])

#ax.legend(bbox_to_anchor=(.6,.65))
fig.show()

fig.savefig('/Users/coletamburri/Desktop/scatter_final20.png')

bins=np.arange(0,400,7.333333333333333333333333333333333)

fig,ax=plt.subplots(dpi=200)
ax.set_ylim([0,15*5])

for i in range(len(largestind)):
    if i<5:
        ax.hist(widths[i],edgecolor='k',alpha=1,bins=bins,bottom = 15*i,zorder=1,color=clrs[2*i],
                label='Frame '+str(largestind[i]))
        ax.text(.25-.09,((i+1)/5)-0.03, 'Frame '+str(largestind[i]), transform=ax.transAxes, fontsize=5)
        ax.text(.25-.1,((i+1)/5)-0.06, r'$\mu =$ '+str(round(np.nanmean(widths[i]),2))+ 'km', transform=ax.transAxes, fontsize=5)
        ax.text(.25-.095,((i+1)/5)-0.09, r'$r_0 =$ '+str(round(np.nanmean(largestfried[i]),1))+ 'cm', transform=ax.transAxes, fontsize=5)

        ax.axvline(np.nanmean(widths[i]),ymin=0+(i/5),ymax=((i+1)/5),linestyle='--',c='black',linewidth=1)
    if i<10 and i>4:
        offset = 110
        modified_array = [original_value + offset for original_value in widths[i]]
        ax.hist(modified_array,edgecolor='k',alpha=1,bins=bins,bottom = 15*(i-5),zorder=1,color=clrs[2*i],
                label='Frame '+str(largestind[i]))
        ax.axvline(np.nanmean(modified_array),ymin=0+((i-5)/5),ymax=(((i+1)-5)/5),linestyle='--',c='black',linewidth=1)
        ax.text(.5-.09,((i+1-5)/5)-0.03, 'Frame '+str(largestind[i]), transform=ax.transAxes, fontsize=5)
        ax.text(.5-.1,((i+1-5)/5)-0.06, r'$\mu =$ '+str(round(np.nanmean(widths[i]),2))+ 'km', transform=ax.transAxes, fontsize=5)
        ax.text(.5-.095,((i+1-5)/5)-0.09, r'$r_0 =$ '+str(round(np.nanmean(largestfried[i]),1))+ 'cm', transform=ax.transAxes, fontsize=5)

    if i<15 and i>9:
        offset = 220
        modified_array = [original_value + offset for original_value in widths[i]]
        ax.hist(modified_array,edgecolor='k',alpha=1,bins=bins,bottom = 15*(i-10),zorder=1,color=clrs[2*i],
                label='Frame '+str(largestind[i]))
        ax.axvline(np.nanmean(modified_array),ymin=0+((i-10)/5),ymax=(((i+1)-10)/5),linestyle='--',c='black',linewidth=1)
        ax.text(.75-.09,((i+1-10)/5)-0.03, 'Frame '+str(largestind[i]), transform=ax.transAxes, fontsize=5)
        ax.text(.75-.1,((i+1-10)/5)-0.06, r'$\mu =$ '+str(round(np.nanmean(widths[i]),2))+ 'km', transform=ax.transAxes, fontsize=5)
        ax.text(.75-.095,((i+1-10)/5)-0.09, r'$r_0 =$ '+str(round(np.nanmean(largestfried[i]),1))+ 'cm', transform=ax.transAxes, fontsize=5)

    if i<20 and i>14:
        offset = 330
        modified_array = [original_value + offset for original_value in widths[i]]
        ax.hist(modified_array,edgecolor='k',alpha=1,bins=bins,bottom = 15*(i-15),zorder=1,color=clrs[2*i],
                label='Frame '+str(largestind[i]))
        ax.axvline(np.nanmean(modified_array),ymin=0+((i-15)/5),ymax=(((i+1)-15)/5),linestyle='--',c='black',linewidth=1)
        print(np.nanmean(widths[i]))
        ax.text(1-.09,((i+1-15)/5)-0.03,'Frame '+ str(largestind[i]), transform=ax.transAxes, fontsize=5)
        ax.text(1-.1,((i+1-15)/5)-0.06, r'$\mu =$ '+str(round(np.nanmean(widths[i]),2))+ 'km', transform=ax.transAxes, fontsize=5)
        ax.text(1-.095,((i+1-15)/5)-0.09, r'$r_0 =$ '+str(round(np.nanmean(largestfried[i]),1))+ 'km', transform=ax.transAxes, fontsize=5)

    #ax.axvline(np.nanmean(widths[i]),linestyle='--',c=muted.indigo,linewidth=2)

# jing avg/med - bbso
# ax.axvline(124,linestyle='--',c='black',linewidth=1,label='Jing+2016 (Mean)')
# ax.axvline(124,linestyle='--',c='black',linewidth=1)
# # scullion avg/med - crisp
# #ax.axvline(130,linestyle='-.',c='black',linewidth=1,label='Scullion+2014 (Mode)')
# ax.axvline(130,linestyle='-.',c='black',linewidth=1)
# # brooks average/md? - iris? Other structures
# #ax.axvline(133,linestyle=':',c='black',linewidth=1,label='Brooks+2018 (UFS, Mean)')
# ax.axvline(133,linestyle=':',c='black',linewidth=1)
# # other crisp avg/med post-2014?
# # other bbso avg/med post-2016?
# ax.set_xlim([0,400])
# # crisp pixel resolution
# crispres = 0.059*727
# #ax.axvline(crispres,linestyle='--',c=muted.olive,linewidth=1,label='SST/CRISP Res.')
# ax.axvline(crispres,linestyle='--',c=muted.olive,linewidth=1)
# # iris pixel resolution
# irisres = 0.166*727

# #ax.axvline(irisres,linestyle='--',c=muted.wine,linewidth=1,label='IRIS Res.')
# ax.axvline(irisres,linestyle='--',c=muted.wine,linewidth=1)
# # # aia pixel resolution
# # aiares = 1.5*727
# # ax.axvline(aiares,linestyle='--',c=muted.purple,linewidth=1,label='SDO/AIA')
# # dkist halpha pixel resolution
# dkisthalphares = 0.017*727
# #ax.axvline(dkisthalphares,linestyle='--',c='black',linewidth=1,label=r'DKIST/VBI H$\alpha$ Res.')
# ax.axvline(dkisthalphares,linestyle='--',c='black',linewidth=1)
# # bbso/nst pixel resolution
# nst_visres = 0.03*727 # according to jing+2016
# #ax.axvline(nst_visres,linestyle='--',c=muted.purple,linewidth=1,label=r'BBSO/VIS H$\alpha$ Res.')
# ax.axvline(nst_visres,linestyle='--',c=muted.purple,linewidth=1)
# ax.set_ylabel('Num. Occurrences',font='Tahoma',fontsize=12)


ax.set_xticks([0,40,80,110,150,190,220,260,300,330,370,410],labels=['0','40','80',
                                                                    '0','40','80',
                                                                    '0','40','80',
                                                                    '0','40','80'])

ax.set_yticks([0,10,15,25,30,40,45,55,60,70],labels=['0','10','0','10','0','10','0','10','0','10'])


for i in [110,220,330]:
    ax.axvline(i,c='black')

ax.set_xlabel('Width [km]',fontsize=12,font='Tahoma')

for i in [0,15,30,45,60]:
    ax.axhline(i,c='black')
    
ax.set_xlim([0,440])

#ax.legend(fontsize=5,bbox_to_anchor=(.85,.9))
fig.show()

fig.savefig('/Users/coletamburri/Desktop/histo_final20.png')

widthsall = widths[0].copy()

for i in range(1,len(largestind)):
    widthsall = np.concatenate((widthsall, widths[i]))
    
    
#plotting the fits from one of the frames
# Define the length of the line, in pixels
#for i in range(len()):
    
path = '/Volumes/VBI_External/pid_2_11/'
ind = largestind[2]

folder_vbi = 'AXXJL'
dir_list = os.listdir(path+folder_vbi)
dir_list.sort()
dir_list2 = []

#for i in range(len(dir_list)):
for i in range(300):
    filename = dir_list[i]
    if filename[-5:] == '.fits' and '_I_' in filename:
        dir_list2.append(filename)

dir_list2.sort()

fullhalpha = fits.open(path+folder_vbi+'/'+dir_list2[ind])
#fullhalpha = fits.open(path+folder_vbi+'/'+filename)
frame = fullhalpha[1].data[0]
    
def Gauss_func(x,A,mu,sigma,m,b):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))+ m*x + b

fig,ax=plt.subplots(5,9,dpi=300)

ncol = 45
colormap = tc.tol_cmap(colormap='rainbow_discrete',lut=ncol+4)
cmap_choice = colormap(np.linspace(0,.7,ncol))


#symmetry test
import pandas as pd
from scipy.stats import skew

def check_symmetry(data):
    """
    Checks if a distribution is symmetric using skewness.

    Parameters:
        data (pd.Series or list): The data to check for symmetry.

    Returns:
        str: "Symmetric", "Right-skewed", or "Left-skewed".
    """
    s = skew(data,nan_policy='omit')
    print(s)
    if -0.5 <= s <= 0.5:
        return "Symmetric"
    elif s > 0.5:
        return "Right-skewed"
    else:
        return "Left-skewed"
    



for i in range(3):
    xlo = xlos[2][i]
    xhi = xhis[2][i]
    ylo = ylos[2][i]
    yhi = yhis[2][i]
    
    framezoom = frame[xlo:xhi,ylo:yhi]
    for j in range(15):   
        x0 = lineloxs[2][15*i+j]
        x1 = linehixs[2][15*i+j]
        # Define the length of the line, in pixels
        y0 = lineloys[2][15*i+j]
        y1 = linehiys[2][15*i+j]
        
        
        length = int(np.hypot(x1-x0,y1-y0))
        
        # Define the x and y in # pixels along the cut
        x, y = np.linspace(y0, y1, length), np.linspace(x0, x1, length)
        
        # Find the intensities nearest to the (x,y) coordinates above
        # This essentially finds the intensity profile along the cut.
        zi = framezoom[x.astype(int), y.astype(int)]
        
        # Use skimage to find the intensity profile along the line.
        # skimage.measure.profile_line returns a 1D array of intensity values 
        # in a directly line from (x0,y0) to (x1,y1), of length equal to the 
        # ceil of the computed length of the line (in units of pixels)
        profile = skimage.measure.profile_line(framezoom,[y0,x0],[y1,x1])
        spatial_samp = 0.017
        # Convert the length of the skimage output to arcsec
        xdirection = np.arange(len(profile))*spatial_samp
        popt = [amps[2][15*i+j],cents[0][15*i+j],stds[0][15*i+j],slopes[0][15*i+j],intercepts[0][15*i+j]]
        if np.isnan(sts[0][15*i+j]) == 0:
            xdirection_finer = np.arange(xdirection[int(sts[0][15*i+j])],xdirection[int(ends[0][15*i+j])],.001)
    
            # Plot intensity profile in separate window
           
            #ax.flatten()[(15*i)+j].plot(xdirection[int(sts[0][15*i+j]):int(ends[0][15*i+j])],profile[int(sts[0][15*i+j]):int(ends[0][15*i+j])],'-.',markersize=1,c='black')
            ax.flatten()[(15*i)+j].plot(xdirection,profile,marker='.',markersize=3,c=cmap_choice[15*i+j])
            
            ax.flatten()[(15*i)+j].plot(xdirection_finer,Gauss_func(xdirection_finer,*popt),'-',c='red',alpha=0.9,linewidth=0.9)
            ax.flatten()[(15*i)+j].axes.get_yaxis().set_ticks([])
            ax.flatten()[(15*i)+j].axes.get_xaxis().set_ticks([])
 
        else:
            ax.flatten()[(15*i)+j].plot(xdirection,profile,'--',markersize=1,c=cmap_choice[15*i+j],alpha=0.5)

            ax.flatten()[(15*i)+j].axes.get_yaxis().set_ticks([])
            ax.flatten()[(15*i)+j].axes.get_xaxis().set_ticks([])  
        fig.tight_layout()
     
normalized = (frame-frame.min()) /(frame.max() -frame.min())
fig,ax=plt.subplots(dpi=400);
ax.imshow(np.log10(normalized[xlos[2][0]:xhis[2][0],ylos[2][0]:yhis[2][0]]),cmap=matplotlib.colormaps['afmhot'],vmax=np.log(1.2))#,vmin=np.log10(0.12),vmax=np.log10(0.7),alpha=0.8)

for i in range(15):
    ax.plot([lineloxs[2][i],linehixs[2][i]],[lineloys[2][i], linehiys[2][i]],c=cmap_choice[i])
plt.show()

fig,ax=plt.subplots(dpi=400);
ax.imshow(np.log10(normalized[xlos[2][1]:xhis[2][1],ylos[2][1]:yhis[2][1]]),cmap=matplotlib.colormaps['afmhot'],vmax=np.log(1.2))#,vmin=np.log10(0.12),vmax=np.log10(0.7),alpha=0.8)

for i in range(15):
    ax.plot([lineloxs[2][15+i],linehixs[2][15+i]],[lineloys[2][15+i], linehiys[2][15+i]],c=cmap_choice[15+i])
plt.show()

fig,ax=plt.subplots(dpi=400);
ax.imshow(np.log10(normalized[xlos[2][2]:xhis[2][2],ylos[2][2]:yhis[2][2]]),cmap=matplotlib.colormaps['afmhot'],vmax=np.log(1.2))#,vmin=np.log10(0.12),vmax=np.log10(0.7),alpha=0.8)

for i in range(15):
    ax.plot([lineloxs[2][30+i],linehixs[2][30+i]],[lineloys[2][30+i], linehiys[2][30+i]],c=cmap_choice[30+i])
plt.show()
        
fig,ax=plt.subplots(dpi=400,figsize=(10,10))
ax.imshow(np.log10(normalized),cmap=matplotlib.colormaps['afmhot'],vmin=np.log10(.2),vmax=np.log10(0.92))
ax.set_aspect('equal')
for i in range(len(ylos[1])):
    rect = patches.Rectangle((ylos[2][i], xlos[2][i]), yhis[2][i]-ylos[2][i], xhis[2][i]-xlos[2][i], linewidth=1, edgecolor='r', \
                              facecolor='none')
    ax.add_patch(rect)
ax.set_xticks([])
ax.set_yticks([])
plt.show()

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    