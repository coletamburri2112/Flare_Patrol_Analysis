#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 06:18:46 2024

@author: coletamburri
"""


import numpy as np
import dkistpkg_ct as DKISTanalysis
import matplotlib
import matplotlib.pyplot as plt

import pandas as pd


from nltk.cluster import KMeansClusterer
import nltk

# loads file containing times and 3D spectra (time, dispersion, spatial)
nsteps = 91
line = 0 #  0 for caii/hepsilon, 1 for hbeta


#filename = '/Users/coletamburri/Desktop/August_2024_DKIST_Flares/8AugXclass_Hbeta.npz'
#filename = '/Users/coletamburri/Desktop/August_2024_DKIST_Flares/8AugXclass_caII_hep.npz'
#filename = '/Users/coletamburri/Desktop/ViSPselection11August24Mclass.npz'
#filename = '/Users/coletamburri/Desktop/Misc_DKIST/11August2024_Cclass_imp_CaII.npz'
#filename = '/Users/coletamburri/Desktop/Misc_DKIST/CaII_Hep_Cclass_11Aug2024.npz'
#coord_filename='/Users/coletamburri/Desktop/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/ViSP_coalign_result_11Aug_Cclass'
if line ==1:
    filename = '/Users/coletamburri/Desktop/11_Aug_2024_Cclass_flare/Processed_ViSP_VBI_11Aug2024/ViSP_spectra_processed_11Aug24_Hbeta.npz'
if line ==0:
    filename = '/Users/coletamburri/Desktop/11_Aug_2024_Cclass_flare/Processed_ViSP_VBI_11Aug2024/ViSP_spectra_processed_11Aug24_CaII.npz'
res = np.load(filename)

# coordres = np.load(coord_filename)

# #only for ca II 
# vispx = coordres['arr_0']
# vispy = coordres['arr_1']

flare_arr = res['flare']
wave = res['wl']
#times = res['arr_1']



hbeta_low =443
hbeta_high = 730

caII_low = 570
caII_high = 775

hepsilon_low = 775
hepsilon_high = 900

dkist_coord_file = '/Users/coletamburri/Desktop/ViSPcoords.npz'
dkist_coords = np.load(dkist_coord_file)

xarr_caII = dkist_coords['xarr_caII']
yarr_caII = dkist_coords['yarr_caII']

xarr_hbeta = dkist_coords['xarr_hbeta']
yarr_hbeta = dkist_coords['yarr_hbeta']

#cutoff0 = 1.5 # for more than one frame
if line == 1:
    cutoff0=2.5  # for h-beta
if line == 0:
    cutoff0=2
#cutoff0 = 2.2 # factor of minimum- 1 means all pixels, >1 is search for flare #1.2 works for hbeta #
#cutoff0=2.6 # for hepsilon

if line == 1:
    n_clusters0 = 35 # 10 works for hbeta, 6 for Ca II H seems to be all that's needed, 6 also for h-ep
if line == 0:
    n_clusters0 = 35

nframes = 1

if line == 1:
    startspace = 300 # 500 for ca ii
    endspace = 1800 # 1500 for ca ii
if line == 0:
    startspace = 300 # 500 for ca ii
    endspace = 1700 # 1500 for ca ii
nsteps = 91
#start = 148 #148 for saved Ca II H/Hepsilon files
start = 148 #143 for saved Hbeta spectra

#cent = 396.85
#cent=397.01

if line == 1:
    cent=486.1375
if line == 0:
    cent = 396.847

# change based on line

if line == 1:
    linelow = hbeta_low
    linehigh = hbeta_high
if line == 0:
    linelow = caII_low
    linehigh = caII_high
    
obs_avg_line = np.mean(flare_arr[start:start+(nsteps*nframes),linelow:linehigh,startspace:endspace],1)
flare_arr2 = flare_arr[start:start+(nsteps*nframes),:,startspace:endspace]

def normalize(data):
    normarr=(data-np.nanmin(data))/(np.nanmax(data)-np.nanmin(data)) 
    return normarr

# simple to compare clusters for each image frame for chosen line
def kmeans_nltk(start,masknum,nsteps,startspace,endspace,obs_avg,flarearr,
                normalize,num_clusters,cutoff,line_low=linelow,
                line_high=linehigh,metric='euclidean',normflag=0):

    frame_line = obs_avg
    cut = cutoff*np.nanmedian(frame_line)
    masklim = cut

    mask = np.copy(frame_line)
    mask[mask < masklim] = 0
    mask[mask > masklim] = 1

    maskinds = np.where(mask > .5)
    x_mask = maskinds[0]
    y_mask = maskinds[1]

    line_profiles = []
    
    for i in range(len(x_mask)):
        line_profiles.append(flarearr[x_mask[i],line_low:line_high,y_mask[i]])
        
    normprofiles_line = []

    if normflag == 1:
        for i in range(len(x_mask)):
            line_norm = normalize(line_profiles[i])
        
            normprofiles_line.append(line_norm)
    
        arr_normprofs = np.asarray(normprofiles_line)
        
        km = KMeansClusterer(num_clusters, nltk.cluster.util.euclidean_distance,
                             repeats=10)
        clusters = km.cluster(arr_normprofs,assign_clusters=True)
        
    elif normflag == 0:
        for i in range(len(x_mask)):
            line_norm = line_profiles[i]
        
            normprofiles_line.append(line_norm)
    
        arr_normprofs = np.asarray(normprofiles_line)
        
        km = KMeansClusterer(num_clusters, nltk.cluster.util.euclidean_distance,
                             repeats=10)
        clusters = km.cluster(arr_normprofs,assign_clusters=True)
    
    
    return frame_line, mask, km, normprofiles_line, clusters, x_mask, y_mask

frame_line, mask0, km0, normprofiles_line, groups0, x_mask0, y_mask0 = \
    kmeans_nltk(start,0,nsteps,startspace,
                endspace,obs_avg_line,flare_arr2,normalize,n_clusters0,cutoff0,normflag=1)
    
arr_normprofs0 = normprofiles_line
colors = plt.cm.jet(np.linspace(0,1,n_clusters0))

aa_arr = [[1836.46057348, 2326.73014872],
       [  11.84946237, 1299.16774314],
       [2045.05316607, 2326.73014872],
       [  44.67204301, 1299.16774314],
       [1836.46057348, 1492.19051546],
       [  11.84946237,  850.67365452]]


## for transformation to vbi only
#x_mask0[i] = 91-x_mask0[i]

#transformation - variation on the function in dkistpkg_ct, without plotting
#x_mask_t,y_mask_t = DKISTanalysis.vbi_visp_transformation(aa_arr,x_mask0,y_mask0,matplotlib,d1=1)

# order by width of curves

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def find_30p_height(curve,find_nearest):
    p30int = np.max(curve)*.3

    lowind,lowval = find_nearest(curve[0:round(len(curve)/2)],p30int)
    highind,highval = find_nearest(curve[round(len(curve)/2):],p30int)
    highind=highind+round(len(curve)/2)
    dist = highind-lowind
    return dist

def find_relint(curve,find_nearest):
    lowind,lowval = find_nearest(curve[0:round(len(curve)/2)],np.nanmax(curve[0:round(len(curve)/2)]))
    highind,highval = find_nearest(curve[round(len(curve)/2):],np.nanmax(curve[round(len(curve)/2):]))
    highind=highind+round(len(curve)/2)
    relint = highval/lowval
    return relint

def find_weightmean(curve,find_nearest):
    values = np.linspace(0,len(curve),len(curve))
    
    # Corresponding weights for each data point
    # These weights could represent the importance or frequency of each point
    weights = curve
    
    # Calculate the weighted mean using numpy.average()
    weighted_mean = np.average(values, weights=weights)
    return weighted_mean

def blue_to_core(curve,hbeta_low=hbeta_low,hbeta_high=hbeta_high,blue=510,core=566,red=580):
    values = np.linspace(0,len(curve),len(curve))
    
    # Corresponding weights for each data point
    # These weights could represent the importance or frequency of each point
    blue_intensity = curve[blue-hbeta_low]
    core_intensity = curve[core-hbeta_low]
    red_intensity = curve[red-hbeta_low]
    ratio = (core_intensity/blue_intensity)+(red_intensity/core_intensity)
    return ratio

dists=[]
for i in range(len(km0.means())):
    dists.append(find_30p_height(km0.means()[i],find_nearest))

relint=[]
for i in range(len(km0.means())):
    relint.append(find_relint(km0.means()[i],find_nearest))
    
wm=[]
for i in range(len(km0.means())):
    wm.append(find_weightmean(km0.means()[i],find_nearest))
    
bc_int=[]
for i in range(len(km0.means())):
    bc_int.append(blue_to_core(km0.means()[i]))

inds = np.arange(len(km0.means()))

#df = pd.DataFrame({'x':inds,'y':dists}) # by distance
#df = pd.DataFrame({'x':inds,'y':relint}) # by relint

if line == 0:
    df = pd.DataFrame({'x':inds,'y':wm}) # by relint
if line == 1:
    df = pd.DataFrame({'x':inds,'y':bc_int}) # by blue wing to core - 480 to 600

df.sort_values(by=['y'])

sortedinds = df.sort_values(['y'])['x']
sortedinds=np.asarray(sortedinds)

distlocs = []

for i in range(len(groups0)):
    distlocs.append(np.where(sortedinds==groups0[i])[0][0])
    
colors = plt.cm.turbo(np.linspace(0,1,n_clusters0))

fig,ax=plt.subplots(figsize=(1.5,4),dpi=200)

if line == 1:
    xarr_ch = xarr_hbeta
    yarr_ch = yarr_hbeta[startspace:endspace+1]
if line==0:
    xarr_ch = xarr_caII
    yarr_ch = yarr_caII[startspace:endspace+1]


ax.pcolormesh(xarr_ch,yarr_ch,np.transpose(frame_line[:,:]),cmap='grey',alpha=1)
ax.scatter(xarr_ch[x_mask0],yarr_ch[y_mask0],2,color=colors[distlocs],alpha=.6,marker='s')
ax.invert_xaxis()
ax.invert_yaxis()

##for transformation only
# else:
#     ax.pcolormesh(vispx,vispy,np.transpose(frame_line[:-2,:]),cmap='grey',alpha=1)
#     ax.scatter(x_mask_t,y_mask_t,5,color=colors[distlocs],alpha=.6,marker='s')
#     ax.invert_yaxis()
ax.set_ylim([-251,-218])
ax.set_xlim([759,766])
# # for tranformation
# ax.set_ylim([2800,800])
# ax.set_xlim([1750,2100])
ax.tick_params(axis='y', labelrotation=90)
ax.set_ylabel('DKIST HPC-y [arcsec]',fontsize=6)
ax.set_xlabel('DKIST HPC-x [arcsec]',fontsize=6)
ax.tick_params(axis='x',labelsize=6)
ax.tick_params(axis='y',labelsize=6)

fig.tight_layout()
fig.show()

#fig,ax=plt.subplots(3,4,figsize=(5,4),dpi=200)
#fig,ax=plt.subplots(2,3,figsize=(5,4),dpi=200) #if hep
fig,ax=plt.subplots(5,7,figsize=(10,6),dpi=200) #if hep and caii
arr_normprofs0 = normprofiles_line

    
for i in range(len(arr_normprofs0)):
    curve = arr_normprofs0[i]
    group = groups0[i]
    ind = np.where(sortedinds==group)[0][0]
    ax.flatten()[ind].plot(wave[linelow:linehigh],curve,alpha=0.01,color='black')
    ax.flatten()[ind].axvline(cent,linewidth=0.6,c='black')

    
for i in range(n_clusters0):
    ax.flatten()[i].plot(wave[linelow:linehigh],km0.means()[sortedinds[i]],marker='*',color=colors[i],markersize=.1)
    ax.flatten()[group].axvline(cent,linewidth=0.6,c='black')
    ax.flatten()[i].tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)
    ax.flatten()[i].tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False)# labels along the bottom edge are off
    ax.flatten()[i].text(0.95, 0.95, str(i+1), transform=ax.flatten()[i].transAxes, \
         ha='right', va='top', fontsize=8, fontweight='bold')
    ax.flatten()[i].set_ylim([-0.2,1.2])
    ax.flatten()[i].set_xlim([wave[linelow],wave[linehigh]])



fig.show()

