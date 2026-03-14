#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  9 09:00:29 2025

@author: coletamburri
"""

import numpy as np
import matplotlib.pyplot as plt
import tol_colors as tc
import pandas as pd


hbeta=0

#filehbeta = '/Users/coletamburri/Desktop/8_August_2024_Xclass_Flare/ViSPselection8AugXclass_hbeta.npz'
filehbeta = '/Users/coletamburri/Desktop/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/ViSP_spectra_processed_11Aug24_Hbeta.npz'
hbetascaled = '/Users/coletamburri/Desktop/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/ViSP_spectra_processed_11Aug24_Hbeta_scaled.npz'
filehbeta = '/Users/coletamburri/Desktop/11Aug2024_calibrated_with_NEWQS_Hbeta_cut.npz'
file = '/Volumes/ViSP_External/CaII_11Aug204_Cclass_newcalib.npz'
#file = '/Users/coletamburri/Desktop/8_August_2024_Xclass_Flare/ViSPselection8AugXclass.npz'
#file = '/Users/coletamburri/Desktop/Misc_DKIST/ViSPselection11August24Mclass.npz'


caII_low = 570
caII_high = 730
hep_low = 730
hep_high = 900

hep_inner_low = 810
hep_inner_high = 825

caii_inner_low = 620
caii_inner_high = 680

hbeta_low = 500
hbeta_high = 650

n_points = 10
spec=0


if hbeta==1:
    data_hbeta = np.load(filehbeta)
    data_hbeta_scaled = np.load(hbetascaled)
    wl = data_hbeta['wl']
    # scaled = data['scaled']
    #flare = data_hbeta['flare']
    flare=data_hbeta['bkgd']
    #scaled = data_hbeta_scaled['scaled']
    time_hbeta = data_hbeta['time']
else:
    data = np.load(file)

    #raw = data['raw']
    wl = data['wl'] # not for X-class
    #scaled = data['scaled']

    flare = data['flare']
    times = data['time']
    

spectra = flare

if hbeta==1:
    hbeta_avg = np.mean(spectra[:,hbeta_low:hbeta_high,:],1)
    hbeta_redwing = flare[:,650,:]
    hbeta_bluewing = flare[:,570,:]
    choice = hbeta_avg
else:

    caII_avg = np.mean(spectra[:,caII_low:caII_high,:],1)
    # hep_avg = np.mean(spectra[:,hep_low:hep_high,:],1)
    # hep_avg_inner = np.mean(spectra[:,hep_inner_low:hep_inner_high,:],1)
    
    #caii_avg_inner = np.mean(spectra[:,caii_inner_low:caii_inner_high,:],1)
    
    caii_avg_redwing = spectra[:,600,:] #old was 700
    
    caii_avg_bluewing = spectra[:,540,:] #old was 622
    caii_avg_core = spectra[:,570,:] # old was 650

    
    # both_avg = np.mean(spectra[:,caII_low:hep_high,:],1)
    
    # all_avg = np.mean(spectra,1)
    choice=caii_avg_redwing
    

if spec == 1: # just want to look at spectra
    npoints = 3
    fig,ax=plt.subplots();
    ax.pcolormesh(np.transpose(np.mean(spectra,1)))
    fig.show()
    
    xx = plt.ginput(npoints,timeout=120)
    
 
    
    fig,ax=plt.subplots(1,npoints,dpi=200);
    
    for i in range(npoints):
        ax.flatten()[i].pcolormesh(np.transpose(spectra[int(xx[i][0]),:,:]),cmap='Reds')
        ax.flatten()[i].set_xlim([400,915])
        ax.flatten()[i].set_xticks([450,550,650,750,850])
        ax.flatten()[i].set_yticks([])

        labels = [str(round(wl[450],2)),str(round(wl[550],2)),str(round(wl[650],2)),str(round(wl[750],2)),str(round(wl[850],2))]
        ax.flatten()[i].set_xticklabels(labels, rotation=45)
    
else:


    

    
    
    fig,ax=plt.subplots()
    ax.pcolormesh(np.transpose(choice),cmap='magma')
    ax.invert_yaxis()
    fig.show()
    
    aa = plt.ginput(2,timeout=120)
    
    xlo, xhi, ylo, yhi = int(aa[0][0]), int(aa[1][0]), int(aa[0][1]),\
        int(aa[1][1])
        
    
    fig,ax=plt.subplots()
    ax.pcolormesh(np.transpose(choice[xlo:xhi,ylo:yhi]),cmap='magma')
    ax.invert_xaxis()
    ax.invert_yaxis()
    fig.show()
    
    colors = plt.cm.jet(np.linspace(0,1,n_points))
    
    cc = plt.ginput(n_points,timeout = 120)
    
    #fig,ax=plt.subplots(len(cc),1)
    fig,ax=plt.subplots(figsize=(12,6),dpi=100)
    for i in range(len(cc)):
        xsel,ysel = cc[i][0],cc[i][1]
        #ax.flatten()[i].plot(spectra[int(xsel)+xlo,:,int(ysel)+ylo],color=colors[i])
        ax.plot(wl,flare[int(xsel)+xlo,:,int(ysel)+ylo],color=colors[i],linewidth=4)
        #ax.flatten()[i].plot(flare[int(cc[-1][0])+xlo,:,int(cc[-1][1])+ylo],color='black',alpha=0.5)
        
        #ax.flatten()[i].plot(wl,flare[int(xsel)+xlo,:,int(ysel)+ylo],color=colors[i])
        if hbeta == 0:
            ax.axvline(396.847,linewidth=2,linestyle='dashed',color='grey')
            ax.axvline(397.01,linewidth=2,linestyle='dashed',color='#CC6677')
            ax.set_xlim([396.7,397.1])
            ax.set_ylim([-.1e6,7e6])
        ax.grid('on')
        # if hbeta == 1:
        #     ax.flatten()[i].axvline(486.14)
        #     ax.flatten()[i].set_xlim([486.14-0.6,486.14+.6])
        #     ax.flatten()[i].set_ylim([-.1e6,6e6])
    ax.set_xlabel('Wavelength [nm]',fontsize=12)
    ax.set_ylabel(r'Intensity [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$ sr$^{-1}$]',fontsize=12)

    fig.show()
        
    
    # Extract coordinates
    xlo, xhi, ylo, yhi = int(aa[0][0]), int(aa[1][0]), int(aa[0][1]),\
        int(aa[1][1])
        
    # Extract zoomed-in frame from coordinates
    framezoom = choice[xlo:xhi,ylo:yhi]
    
    
    dkist_coord_file = '/Users/coletamburri/Desktop/11_Aug_2024_Cclass_Flare/Processed_ViSP_VBI_11Aug2024/ViSPcoords.npz'
    dkist_coords = np.load(dkist_coord_file)

    xarr_caII = dkist_coords['xarr_caII']
    yarr_caII = dkist_coords['yarr_caII']

    xarr_hbeta = dkist_coords['xarr_hbeta']
    yarr_hbeta = dkist_coords['yarr_hbeta']
    
    if hbeta == 1:
        startspace = 300 # 500 for ca ii
        endspace = 1800 # 1500 for ca ii
    if hbeta == 0:
        startspace = 300 # 500 for ca ii
        endspace = 1700 # 1500 for ca ii
    
    if hbeta == 1:
        xarr_ch = xarr_hbeta[:-1]
        yarr_ch = yarr_hbeta[ylo:yhi]
    if hbeta==0:
        xarr_ch = xarr_caII[:-1]
        yarr_ch = yarr_caII[ylo:yhi]
    
    maskind = {'x': cc[:][0], 'y': cc[:][1]}
    df_mask = pd.DataFrame(maskind)
    
    fig,ax=plt.subplots(2,5,figsize=(20,7),dpi=200)

    for i in range(10):
        ax.flatten()[i].pcolormesh(xarr_ch,yarr_ch,np.transpose(choice[57+(91*i):56+(91*(i+1)),ylo:yhi-1]),cmap='grey',alpha=1)
        lower_threshold=57+(91*i)
        upper_threshold=56+(91*(i+1))
        xsel = cc[i][0]
        ysel = cc[i][1]
        xsel=xsel-lower_threshold+xlo

        ax.flatten()[i].scatter(xarr_ch[int(xsel)],yarr_ch[int(ysel)],40,color=colors[i],alpha=1,marker='x')
        #ax.flatten()[i].invert_xaxis()
        #ax.flatten()[i].invert_yaxis()
        #ax.flatten()[i].set_ylim([yhi,ylo])
        #ax.flatten()[i].set_ylim([-251,-218])
        #ax.flatten()[i].set_xlim([759,766])
        ax.flatten()[i].tick_params(axis='x', labelrotation=45)
        ax.flatten()[i].tick_params(axis='y', labelrotation=45)
        
        ax.flatten()[i].tick_params(axis='x',labelsize=6)
        ax.flatten()[i].tick_params(axis='y',labelsize=6)
        ax.flatten()[i].set_title(times[57+91*i][11:19],fontsize=8)
        #ax.flatten()[i].set_xticks([760,763,766])
    ax.flatten()[0].set_ylabel('DKIST HPC-y [arcsec]',fontsize=6)
    ax.flatten()[5].set_ylabel('DKIST HPC-y [arcsec]',fontsize=6)

    ax.flatten()[7].set_xlabel('DKIST HPC-x [arcsec]',fontsize=6)
    fig.tight_layout(pad=2.0) 
    fig.subplots_adjust(hspace=0.5,wspace=0.3)


    specmap = []
    
        
    fig,ax=plt.subplots()
    ax.pcolormesh(np.transpose(choice[xlo:xhi,ylo:yhi]),cmap='magma')
    for i in range(len(cc)):
        xsel,ysel = cc[i][0],cc[i][1]
        specmap.append(flare[int(xsel)+xlo,:,int(ysel)+ylo])
        ax.scatter(xsel,ysel,color=colors[i])
    ax.invert_xaxis()
    ax.invert_yaxis()
    fig.show()
    
    fig,ax=plt.subplots();
    X=wl
    Y=np.arange(n_points)

    ax.pcolormesh(X,Y,specmap,cmap='Reds')
    ax.set_xlim(396.8,397.1)
    ax.axvline(396.847)
    ax.axvline(397.01)
        
        
    
