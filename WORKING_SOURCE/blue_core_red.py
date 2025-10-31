#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 13:49:03 2025

@author: coletamburri
"""
inds = [600-500,650-500,700-500]

xlow =500
xhigh = 1500

savefolder = '/Users/coletamburri/Desktop/blue_core_red_nonorm/'

def min_max_normalize(arr):
    arr_min = np.min(arr)
    arr_max = np.max(arr)
    if arr_max - arr_min == 0:  # Handle cases where all values are the same
        return np.zeros_like(arr)
    normalized_arr = (arr - arr_min) / (arr_max - arr_min)
    return normalized_arr

spectranorm=min_max_normalize(spectra[:,500:900,:])
cmap = 'seismic'
fig,ax=plt.subplots(1,3,dpi=200,figsize=(5,10));
for i in range(3):
    ax.flatten()[i].pcolormesh(np.transpose(spectranorm[148:239,inds[i],xlow:xhigh]),cmap=cmap,vmin=0,vmax=1);
    ax.flatten()[i].invert_xaxis();
    ax.flatten()[i].invert_yaxis();
    ax.flatten()[i].set_xticks([])
    ax.flatten()[i].set_yticks([])
    plt.subplots_adjust(wspace=0)
    fig.suptitle(time[148])
fig.show()
fig.savefig(savefolder+'1.png')
    

    
fig,ax=plt.subplots(1,3,dpi=200,figsize=(5,10));
for i in range(3):
    ax.flatten()[i].pcolormesh(np.transpose(spectranorm[239:330,inds[i],xlow:xhigh]),cmap=cmap,vmin=0,vmax=1);
    ax.flatten()[i].invert_xaxis();
    ax.flatten()[i].invert_yaxis();
    ax.flatten()[i].set_xticks([])
    ax.flatten()[i].set_yticks([])
    plt.subplots_adjust(wspace=0)
    fig.suptitle(time[239])
fig.show()
fig.savefig(savefolder+'2.png')

    
fig,ax=plt.subplots(1,3,dpi=200,figsize=(5,10));
for i in range(3):
    ax.flatten()[i].pcolormesh(np.transpose(spectranorm[330:421,inds[i],xlow:xhigh]),cmap=cmap,vmin=0,vmax=1);
    ax.flatten()[i].invert_xaxis();
    ax.flatten()[i].invert_yaxis();
    ax.flatten()[i].set_xticks([])
    ax.flatten()[i].set_yticks([])
    plt.subplots_adjust(wspace=0)
    plt.subplots_adjust(wspace=0)
    fig.suptitle(time[330])
fig.show()
fig.savefig(savefolder+'3.png')

    
    
fig,ax=plt.subplots(1,3,dpi=200,figsize=(5,10));
for i in range(3):
    ax.flatten()[i].pcolormesh(np.transpose(spectranorm[421:512,inds[i],xlow:xhigh]),cmap=cmap,vmin=0,vmax=1);
    ax.flatten()[i].invert_xaxis();
    ax.flatten()[i].invert_yaxis();
    ax.flatten()[i].set_xticks([])
    ax.flatten()[i].set_yticks([])
    plt.subplots_adjust(wspace=0)
    fig.suptitle(time[421])

fig.show()
fig.savefig(savefolder+'4.png')

    
fig,ax=plt.subplots(1,3,dpi=200,figsize=(5,10));
for i in range(3):
    ax.flatten()[i].pcolormesh(np.transpose(spectranorm[512:603,inds[i],xlow:xhigh]),cmap=cmap,vmin=0,vmax=1);
    ax.flatten()[i].invert_xaxis();
    ax.flatten()[i].invert_yaxis();
    ax.flatten()[i].set_xticks([])
    ax.flatten()[i].set_yticks([])
    plt.subplots_adjust(wspace=0)
    fig.suptitle(time[512])

fig.show()
fig.savefig(savefolder+'5.png')

    
fig,ax=plt.subplots(1,3,dpi=200,figsize=(5,10));
for i in range(3):
    ax.flatten()[i].pcolormesh(np.transpose(spectranorm[603:694,inds[i],xlow:xhigh]),cmap=cmap,vmin=0,vmax=1);
    ax.flatten()[i].invert_xaxis();
    ax.flatten()[i].invert_yaxis();
    ax.flatten()[i].set_xticks([])
    ax.flatten()[i].set_yticks([])
    plt.subplots_adjust(wspace=0)
    fig.suptitle(time[603])
fig.show()
fig.savefig(savefolder+'6.png')

    
fig,ax=plt.subplots(1,3,dpi=200,figsize=(5,10));
for i in range(3):
    ax.flatten()[i].pcolormesh(np.transpose(spectranorm[694:785,inds[i],xlow:xhigh]),cmap=cmap,vmin=0,vmax=1);
    ax.flatten()[i].invert_xaxis();
    ax.flatten()[i].invert_yaxis();
    ax.flatten()[i].set_xticks([])
    ax.flatten()[i].set_yticks([])
    plt.subplots_adjust(wspace=0)
    fig.suptitle(time[694])
fig.show()
fig.savefig(savefolder+'7.png')

    
#later scan

fig,ax=plt.subplots(1,3,dpi=200,figsize=(5,10));
for i in range(3):
    ax.flatten()[i].pcolormesh(np.transpose(spectranorm[786:876,inds[i],xlow:xhigh]),cmap=cmap,vmin=0,vmax=1);
    ax.flatten()[i].invert_xaxis()
    ax.flatten()[i].invert_yaxis()
    ax.flatten()[i].set_xticks([])
    ax.flatten()[i].set_yticks([])
    plt.subplots_adjust(wspace=0)
    fig.suptitle(time[786])
fig.show()
fig.savefig(savefolder+'8.png')

    
fig,ax=plt.subplots(1,3,dpi=200,figsize=(5,10));
for i in range(3):
    ax.flatten()[i].pcolormesh(np.transpose(spectranorm[876:967,inds[i],xlow:xhigh]),cmap=cmap,vmin=0,vmax=1);
    ax.flatten()[i].invert_xaxis()
    ax.flatten()[i].invert_yaxis()
    ax.flatten()[i].set_xticks([])
    ax.flatten()[i].set_yticks([])
    plt.subplots_adjust(wspace=0)
    fig.suptitle(time[786])
fig.show()
fig.savefig(savefolder+'9.png')
