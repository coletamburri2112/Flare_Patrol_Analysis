#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 06:18:46 2024

@author: coletamburri
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import scipy as scp
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from tslearn.clustering import TimeSeriesKMeans
import sklearn

from nltk.cluster import KMeansClusterer
import nltk

from pyclustering.cluster.kmeans import kmeans
from pyclustering.utils.metric import distance_metric
from pyclustering.cluster.center_initializer import random_center_initializer
from pyclustering.cluster.encoder import type_encoding
from pyclustering.cluster.encoder import cluster_encoder

# loads file containing times and 3D spectra (time, dispersion, spatial)
nsteps = 91

filename = '/Users/coletamburri/Desktop/11AugCclass_caii.npz'
res = np.load(filename)

flare_arr = res['arr_0']
#times = res['arr_1']

hbeta_low =480
hbeta_high = 670

caII_low = 570
caII_high = 730

hepsilon_low = 750
hepsilon_high = 900

cutoff0 = 3 # factor of minimum- 1 means all pixels, >1 is search for flare

n_clusters0 = 10

nframes = 10
startspace = 500
endspace = 1500
nsteps = 91
start = 149

# change based on line

linelow = caII_low
linehigh = caII_high

obs_avg_line = np.mean(flare_arr[:,linelow:linehigh,:],1)

def normalize(data):
    normarr=(data-np.nanmin(data))/(np.nanmax(data)-np.nanmin(data)) 
    return normarr

# simple to compare clusters for each image frame for chosen line
def kmeans_nltk(start,masknum,nsteps,startspace,endspace,obs_avg,flarearr,
                normalize,num_clusters,cutoff,line_low=linelow,
                line_high=linehigh,metric='euclidean'):

    frame_line = obs_avg[(start+(masknum)*nsteps):(start+(masknum+1)*nsteps),:]
  

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
        line_profiles.append(flare_arr[(start*(masknum+1))+x_mask[i],
                                       line_low:line_high,
                                       (start*(masknum+1))+y_mask[i]])
        
    normprofiles_line = []

    for i in range(len(x_mask)):
        line_norm = normalize(line_profiles[i])
    
        normprofiles_line.append(line_norm)

    arr_normprofs = np.asarray(normprofiles_line)
    
    km = KMeansClusterer(num_clusters, nltk.cluster.util.euclidean_distance,
                         repeats=10)
    clusters = km.cluster(arr_normprofs,assign_clusters=True)
    
    
    return frame_line, mask, km, normprofiles_line, clusters, x_mask, y_mask

frame_line, mask0, km0, normprofiles_line, groups0, x_mask0, y_mask0 = \
    kmeans_nltk(start,0,nsteps,startspace,
                endspace,obs_avg_line,flare_arr,normalize,n_clusters0,cutoff0)
    
fig,ax=plt.subplots(n_clusters0,1,figsize=(10,10))
arr_normprofs0 = normprofiles_line
colors = plt.cm.jet(np.linspace(0,1,n_clusters0))

for i in range(len(arr_normprofs0)):
    curve = arr_normprofs0[i]
    group = groups0[i]
    
    ax.flatten()[group].plot(curve,alpha=0.01,color='black')

for i in range(n_clusters0):
    ax.flatten()[i].plot(km0.means()[i],marker='*',color=colors[i])
    
fig,ax=plt.subplots(figsize=(5,10))
ax.pcolormesh(np.transpose(frame_line),cmap = 'hot',alpha=0.5)
ax.scatter(x_mask0,y_mask0,2,color=colors[groups0],alpha=1)

fig,ax=plt.subplots(5,4,figsize=(10,10))
arr_normprofs0 = normprofiles_line
colors = plt.cm.jet(np.linspace(0,1,n_clusters0))

for i in range(len(arr_normprofs0)):
    curve = arr_normprofs0[i]
    group = groups0[i]
    
    ax.flatten()[group].plot(curve,alpha=0.01,color='black')

for i in range(n_clusters0):
    ax.flatten()[i].plot(km0.means()[i],marker='*',color=colors[i])


