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
import KMeansFunctions as KMF
from matplotlib.collections import LineCollection


class SpectralClusterer:
    
    def __init__(self, data, pid, wave, times):
        self.data = np.asarray(data)
        
        # Restrict the input array to be 4D - [map num, wave, ROI x, ROI y]
        if self.data.ndim != 4 and pid == '4.17': 
            raise ValueError('Datacube must be 4D')
        elif self.data.ndim != 3 and pid == '2.11':
            raise ValueError('Datacube must be 3D')

        if pid == '4.17':
            # Get the input lengths in each dimension
            self.n_maps = self.data.shape[0]
            self.n_wave = self.data.shape[1]
            self.n_x = self.data.shape[2]
            self.n_y = self.data.shape[3]
        
            # Placeholder for wavelength array
            self.wave = np.arange(self.n_wave)
                
            # Change the shape of the dataset to be 3D (squish all maps together)
            self.X_pre = np.concatenate(self.data, 1)
        
        if pid == '2.11':
            self.n_x = self.data.shape[0]
            self.n_wave = self.data.shape[1]
            self.n_y = self.data.shape[2]
            
            self.wave = wave
            
            self.X_pre = self.data
            self.times = times
        
        # Initialize the results array
        self.results = {}
        
    def prep_data(self, cutoff, cent, coreind, line_low = 20, line_high = -20, 
                  normflag = 1, frame_line_metric = 'core', pid = '2.11',
                  startspace = 0, endspace=None, start=0, nsteps=91, nframes=10):
        

        if frame_line_metric == 'core':
            if pid == '4.17':
                self.frame_line = self.X_pre[coreind,:,:]
            elif pid == '2.11':
                self.frame_line = self.X_pre[:,coreind,:]

        elif frame_line_metric == 'mean':
            if pid == '4.17':
                obs_avg = np.mean(self.X_pre[line_low:line_high,
                                                start:start+(nsteps*nframes),
                                                startspace:endspace],0).copy
                
                flare_arr2 = self.X_pre[start:start+(nsteps*nframes),
                                                        :,
                                                        startspace:endspace].copy()
            
                
                self.frame_line = obs_avg
                
            elif pid == '2.11':
                
                obs_avg = np.nanmean(self.X_pre[start:start+(nsteps*nframes),
                                                        line_low:line_high,
                                                        startspace:endspace],1)
                
                flare_arr2 = self.X_pre[start:start+(nsteps*nframes),
                                                        :,
                                                        startspace:endspace].copy()
                
                self.frame_line = obs_avg
                
        self.line_low = line_low
        self.line_high = line_high
        self.cent = cent
        
        cut = cutoff*np.nanmedian(self.frame_line)
        masklim = cut
        self.mask = np.copy(self.frame_line)
        self.mask[self.mask < masklim] = 0
        self.mask[self.mask > masklim] = 1
        
        maskinds = np.where(self.mask > .5)
        
        self.x_mask = maskinds[0]
        self.y_mask = maskinds[1]
        
        self.line_profiles = []
        
        if pid == '4.17':
            self.line_profiles = flare_arr2[line_low:line_high,
                                            self.x_mask,
                                            self.y_mask].T
        elif pid == '2.11':
            self.line_profiles = flare_arr2[self.x_mask,
                                            line_low:line_high,
                                            self.y_mask]                
        self.Xr = []
        

        if normflag == 1:
            for i in range(len(self.x_mask)):
                line_norm = normalize(self.line_profiles[i,:])
                self.Xr.append(line_norm)
        else:
            for i in range(len(self.x_mask)):
                self.Xr.append(self.line_profiles[i,:])
                
        
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
    
    def plot_metrics(self):
        print('need to do!')
    
    def sort(self,k_choice,method=['30p_height','WM','orig_index']):
        wm = []
        dists = []
        # select the exact k value you want
        self.labelsint = self.results[k_choice]['labels']
        counts = np.bincount(self.labelsint)
        self.centroidsint = self.results[k_choice]['centroids']
        
        for i in range(len(self.centroidsint)):
            dists.append(KMF.find_30p_height(self.centroidsint[i],\
                                                         KMF.find_nearest))
            wm.append(KMF.find_weightmean(self.centroidsint[i],\
                                                      KMF.find_nearest))
        
        inds = np.arange(len(self.centroidsint))

        self.df = pd.DataFrame({'orig_index':inds,'WM':wm,'30p_height':dists,'count':counts})
                
        self.sortedwls = np.asarray(self.df.sort_values([method])['30p_height'])
        self.sortedinds = np.asarray(self.df.sort_values([method])['orig_index'])
        self.sortedcounts = self.df.sort_values([method])['count']
        
        self.distlocs = [
            np.where(self.sortedinds == i)[0][0]
            for i in self.labelsint
        ]
        
        maskind = {'x': self.x_mask, 'y': self.y_mask,'dist': self.distlocs}
        self.df_mask = pd.DataFrame(maskind)

    def plot_maps(self, k_choice, startspace, endspace, xarr, yarr, nframes,
                  nsteps):
        self.colors = plt.cm.turbo(np.linspace(0,1,k_choice))
        
        xarr_ch = xarr
        yarr_ch = yarr[startspace:endspace+1]
        
        
        ##PLOT 1 - COLORMAP
        
        fig,ax = plt.subplots(1,10,figsize=(8,3),dpi=200)
        axes=ax.ravel()
        for i in range(nframes):
            lower_threshold=(nsteps*i)
            upper_threshold=(nsteps*(i+1))
            sel = self.df_mask[(self.df_mask['x'] > lower_threshold) & \
                               (self.df_mask['x'] < upper_threshold)]
            xmask_sel = sel['x'].values-lower_threshold
            ymask_sel = sel['y'].values
            distlocs_sel = sel['dist'].values
            
            ax[i].pcolormesh(xarr_ch,yarr_ch,np.transpose(self.frame_line[(nsteps*i):(nsteps*(i+1)),:]),\
                          cmap='grey',alpha=1)
            ax[i].scatter(xarr_ch[xmask_sel],yarr_ch[ymask_sel],.2,color=self.colors[distlocs_sel],\
                       alpha=0.6,marker='s')
            ax[i].set_ylim([-251,-218])
            ax[i].set_xlim([759,766])
            ax[i].tick_params(axis='x', labelrotation=45)
            ax[i].tick_params(axis='y', labelrotation=45)
            
            if i>0:
                ax[i].set_yticks([])
            
            ax[i].tick_params(axis='x',labelsize=6)
            ax[i].tick_params(axis='y',labelsize=6)
            ax[i].set_title(self.times[57+91*i][11:19],fontsize=8)
            ax[i].set_xticks([760,763,766])
            
        ##PLOT 2 - CLUSTERS
             
        fig,ax=plt.subplots(int(k_choice/5),5,figsize=(10,6),dpi=200) 
        axes = ax.ravel()
        group_to_ind = {g: i for i, g in enumerate(self.sortedinds)}
        
        x = self.wave[self.line_low:self.line_high]
        
        lines_per_ax = [[] for _ in axes]
        
        for curve, group in zip(self.Xr, self.labelsint):
            ind = group_to_ind[group]
            pts = np.column_stack([x,curve])
            lines_per_ax[ind].append(pts)
        
        b=0 
        
        for a, lines in zip(axes, lines_per_ax):
            lc = LineCollection(
                lines,
                colors='black',
                linewidths = 0.5,
                alpha=0.01)    
            a.add_collection(lc)
            a.axvline(self.cent, linewidth=0.5, c='black')
            a.autoscale
            b+=1            
            
        for i in range(k_choice):
            axes[i].plot(self.wave[self.line_low:self.line_high],self.centroidsint[self.sortedinds[i]],\
                         marker='*',color=self.colors[i],markersize=0.1)
            axes[group].axvline(self.cent,linewidth=0.5,c='black')
            if i > k_choice-6:
                axes[i].tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=True,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=True,
                labelsize=6)
                axes[i].set_xlabel('Wavelength [nm]',fontsize=6)
        
            else: 
                axes[i].tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False)                
            
            axes[i].tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left=False,      # ticks along the bottom edge are off
            right=False,         # ticks along the top edge are off
            labelleft=False)# labels along the bottom edge are off
            axes[i].text(0.95, 0.95, str(i+1), transform=axes[i].transAxes, \
                 ha='right', va='top', fontsize=6, fontstyle='italic')
            axes[i].text(0.3, 0.95, str(round(100*float(self.sortedcounts[self.sortedinds[i]])/\
                                              len(self.labelsint),2)), transform=axes[i].transAxes,\
                 ha='right', va='top', fontsize=6, fontstyle='italic')
            axes[i].set_ylim([-0.2,1.2])
            axes[i].set_xlim([self.wave[self.line_low],self.wave[self.line_high]])     
            
            c = 299792458 # speed of light in m/s
            cent=self.cent
            def veltrans(x,mu=1):
                return ((((x)/cent)-1)*c/1000)/mu

            def veltrans2(x,mu2=1):
                return ((((x)/cent)-1)*c/1000)/mu2

            def wltrans(x):
                return ((((x/(c/1000))+1)*cent))
            
            if i < 5:
                secaxx = axes[i].secondary_xaxis('top', functions=(veltrans,wltrans))
                secaxx.set_xlabel(r'Velocity $[km\; s^{-1}]$',fontsize=6)
                secaxx.tick_params(axis='both', which='major', labelsize=6)
                
                
            fig.subplots_adjust(hspace=0,wspace=0)
        
        

            
    
        
            
