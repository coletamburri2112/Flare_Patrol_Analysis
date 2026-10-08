#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 10:48:01 2026

@author: coletamburri
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn import metrics
import pandas as pd
import dkistpkg_ct
from KMeansFunctions import normalize


class SpectralClusterer:
    
    def __init__(self, data):
        self.data = np.asarray(data)
        
        # Restrict the input array to be 4D - [map num, wave, ROI x, ROI y]
        if self.data.ndim != 4: 
            raise ValueError('Datacube must be 4D')
        
        # Get the input lengths in each dimension
        self.n_maps = self.data.shape[0]
        self.n_wave = self.data.shape[1]
        self.n_x = self.data.shape[2]
        self.n_y = self.data.shape[3]
            
        # Change the shape of the dataset to be 3D (squish all maps together)
        self.X_pre = np.concatenate(self.data, 1)
        
        # Initialize the results array
        self.results = {}
        
    def prep_data(self, cutoff, line_low = 20, line_high = -20, 
                  normflag = 1):
        
        self.frame_line = np.nanmean(self.X_pre[line_low:line_high,
                                                :,:],0)
        
        cut = cutoff*np.nanmedian(self.frame_line)
        
        self.mask = np.copy(self.frame_line)
        self.mask[self.mask < cut] = 0
        self.mask[self.mask > cut] = 1
        
        maskinds = np.where(self.mask > .5)
        
        self.x_mask = maskinds[0]
        self.y_mask = maskinds[1]
        
        self.line_profiles = []
        
        self.line_profiles = self.X_pre[line_low:line_high,
                                        self.x_mask,
                                        self.y_mask].T
                
        self.Xr = []
        
        if normflag == 1:
            for i in range(len(self.x_mask)):
                line_norm = normalize(self.line_profiles[i])
                self.Xr.append(line_norm)
        else:
            for i in range(len(self.x_mask)):
                self.Xr.append(self.line_profiles[i])
                
        
    def clustering(self, k_values, metricstest = 0, n_init = 10):
        
        for k in k_values:
            km = KMeans(n_clusters = k, n_init = n_init).fit(self.Xr)
    
            labels = km.labels_
            cc = km.cluster_centers_
            
            if metricstest == 1:
                inertia = km.inertia_
                silhouette = metrics.silhouette_score(self.Xr, labels)
                DB_index = metrics.davies_bouldin_score(self.Xr, labels)
            
                self.results[k] = {
                    'km_fit': km,
                    'labels':labels,
                    'centroids':cc,
                    'inertia':inertia,
                    'silhouette':silhouette,
                    'DB index':DB_index
                    }
            else:
                self.results[k] = {
                    'km_fit': km,
                    'labels':labels,
                    'centroids':cc,
                    }
                
        
    def get_labels(self, k):
        return self.results[k]['labels']
    
    def get_centroids(self, k):
        return self.results[k]['centroids']
    
    metrics = {}
    
    def get_metrics(self):
        
        for k, result in self.results.items():
            metrics[k] = {
                'inertia':result['inertia'],
                'silhouette':result['silhouette']}
        
        return metrics
            
