#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 09:32:02 2024

@author: coletamburri
Adapted from rhanalyze, c/o Han Uitenbroek
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib


import rhanalyze
from rhanalyze.satlas import satlas

matplotlib.rcParams.update(matplotlib.rcParamsDefault)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def tocgs(w, s):
    clight=2.99792458e10         #speed of light [cm/s]
    joule_2_erg=1e7
    aa_to_cm=1e-8
    s *=joule_2_erg/aa_to_cm # from Watt /(cm2 ster AA) to erg/(s cm2 ster cm)
    s *=(w*aa_to_cm)**2/clight   # to erg/
    return s

def tosi(wav, s):
    clight=2.99792458e8      #speed of light [m/s]                                  
    aa_to_m=1e-10                                                                        
    cm_to_m=1e-2                       
    s /= cm_to_m**2 * aa_to_m # from from Watt /(s cm2 ster AA) to Watt/(s m2 ster m) 
    s *= (wav*aa_to_m)**2 / clight # to Watt/(s m2 Hz ster)
    return s
times = [500] #

for i in range(len(times)):
    time=times[i]+1
    
    if time==3 or time==30 or time == 52 or time == 53:
        continue
    
    base1 = '/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_Oct_2023/RH/'
    base2 = '/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_output_files_npz/'
    
    #define model to read in

    rhd_choice = rhanalyze.rhout('/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_Oct_2023/RH/rhf1d_11Aug_Cclass_A_final_36s/run')

    #rhd_choice = rhanalyze.rhout('/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_Oct_2023/RH/cat_15_8_5e10_wRC_updated_25s_CRD/run')
    #rhd_choice = rhanalyze.rhout('/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_Oct_2023/RH/cat_15_8_5e10_wRC_updated_'+str(time)+'s_CRD/run')
    
    #rhd_choice = rhanalyze.rhout('/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_Oct_2023/RH/rhf1d_15_8_5e10_20_600/run')
    #rhd_choice = rhanalyze.rhout('/Users/coletamburri/Desktop/RH_Versions_and_Tools/RH_Oct_2023/RH/longduration_'+time+'_50broadc_5vt/run')
    # define output fil13
    #filename = base2+'fchroma30_20b_5vt_'+time+'s_H20'
    #filename = base2+'cat_15_8_5e10_wRC_updated_'+str(time)+'s_H6_CRD'
    filename = base2+'11Aug_Cclass_A_final_36s_mu0.5'
    #filename = base2+'longduration_'+time+'s_H6_50broadc'
                                                      
    # generalize - choice of rhd file                                                
    
    
    # source function code
    ## Preload the FTS atlas
    fts = satlas()
    
    # Get parameters of your favorite line and atom
    
    ATOMNO = 1  ## ['0=Ca', 1='H', '2=Na']
    LINENO = 22 ### ['0=Ca IIH';'22-H-epsilon']
    DWING  = 2
    DLAMB  = 2.0
    
    # if Ca II H
    linelow = 396.95
    linehigh = 396.95
    
    # if Hep
    linelow = 396.9
    linehigh = 397.1
    
    #definition of central lambda
    lambda0 = rhd_choice.atoms[ATOMNO].transition[LINENO].lambda0
    print("lambda_0: ", lambda0)
    lambda_blue = lambda0 - DLAMB
    lambda_red  = lambda0 + DLAMB
    
    obs = fts.nmsiatlas(lambda_blue, lambda_red)
    
    #indices defining bounds for line
    index = np.where((rhd_choice.spectrum.waves >= lambda_blue) & \
                     (rhd_choice.spectrum.waves <= lambda_red))
    lblue = index[0][0]
    lred  = index[0][-1]
    icore = np.where(rhd_choice.spectrum.waves == lambda0)
    lcore = icore[0][0]
    lwing = lcore - DWING
    
    ## read the opacity for wavelength index waveno and ray index rayno
    
    #choice of ray
    rayno  = rhd_choice.geometry.Nrays-1
    
    rhd_choice.opacity.read(lcore, rayno)
    
    # Evaluate source function for that wavelength and viewing angle
    # Stored in rhd_choice.opacity.S
    rhd_choice.opacity.Source()
    print("Shape of source function S:", np.shape(rhd_choice.opacity.S))
    
    # Evaluate Planck function for the same wavelength
    # Stored in # Stored in rhd_choice.opacity.Bp
    rhd_choice.opacity.Planck()
    print("Shape of Planck function: Bp", np.shape(rhd_choice.opacity.Bp))
    
    # Stored in rhd_choice.opacity.tau
    # Evaluate optical depth at disk center
    rhd_choice.opacity.get_tau(center=True)
    print("Shape of optical depth tau:", np.shape(rhd_choice.opacity.tau))
    
    #contribution function
    # Define indices for continuum, wing and core
    ldisp = [lblue, lwing, lcore]
    label = ['cont', 'wing', 'core']
    color = ['r', 'b', 'g']
    
    plt.figure(figsize=[3,3],dpi=200)
    
    #intensity
    I = rhd_choice.rays[0].I
    waves = rhd_choice.spectrum.waves
    
    plt.plot(waves[lblue:lred], I[lblue:lred])
    plt.plot(obs[0], obs[1], "y,", label='atlas')
    #plt.xlim([linelow,linehigh])
    
    yminmax = plt.ylim()
    
    #lines defining line wing, core
    for l in range(len(ldisp)):
        plt.plot([waves[ldisp[l]], waves[ldisp[l]]], yminmax, color[l]+'--')
    
    plt.show()
    
    ## Calculate the contribution function for the selected wavelength indices
    plt.figure(figsize=[6,3],dpi=200)
    for l in range(len(ldisp)):
    
        rhd_choice.opacity.read(ldisp[l], rhd_choice.geometry.Nrays-1)
        rhd_choice.opacity.Source()
        rhd_choice.opacity.get_tau(center=True)
    
        ## Derivative of tau wrt height (note negative because index of height is reversed)
        dtaudh = -np.gradient(rhd_choice.opacity.tau, rhd_choice.geometry.height)
        contrib = rhd_choice.opacity.S * np.exp(-rhd_choice.opacity.tau) * dtaudh
    
        plt.plot(rhd_choice.geometry.height/1.0E3, contrib, color[l], label=label[l])
        if l==2:
            locmax = np.argwhere(contrib==np.nanmax(contrib))[0][0]
        #plt.xlim([100,1000])
    
    #here define the location of maximum intensity - for n_e determination
    #plt.axvline(locmax)
        
        
    plt.xlabel('height [km]')
    plt.ylabel('contribution function')
    plt.legend()
    
    plt.show()
    
    KM_TO_M = 1.0E3
    
    height = rhd_choice.geometry.height/ KM_TO_M
    tau500 = rhd_choice.geometry.tau500
    sourcefxn = rhd_choice.opacity.Source()
    
    plt.figure(figsize=[14,6])
    
    plt.subplot(121)
    plt.plot(height, rhd_choice.atmos.T, 'bo')
    plt.plot(height, rhd_choice.atmos.T, 'k', label=rhd_choice.atmos.ID)
    plt.xlabel('height [km]')
    plt.ylabel('T [K]')
    plt.legend()
    
    plt.subplot(122)
    #plt.plot(height, rhd_choice.atmos.nH[:,0], label='nH[0]')
    #plt.plot(height, rhd_choice.atmos.nH[:,5], label='protons')
    plt.plot(height, rhd_choice.atmos.n_elec, label='electrons')
    #plt.xlim([970,1010])
    #plt.xlabel(r'$\tau_{500}$')
    plt.ylabel('log ne [m$^{-3}$]')
    #plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    #plt.axvline(locmax)
    
    
    plt.show()
    
    
    #nearestidx, nearest = find_nearest(height,locmax)
    
    # define the nearest index 
    n_e_maxC = rhd_choice.atmos.n_elec[locmax]
    print(n_e_maxC)
    
    xmu = rhd_choice.geometry.xmu
    
    plt.figure(figsize=[9, 6])
    plt.xlim([396.6, 397.2])
    plt.ylim([-1e-8,5.0e-8])
    
    plt.plot(rhd_choice.spectrum.waves, rhd_choice.rays[0].I,\
             label="{:.2f}".format(rhd_choice.rays[0].muz), linewidth=4.0)
    
    for mu in range(rhd_choice.geometry.Nrays-1, -1, -1):
      xmulabel = "{:.2f}".format(xmu[mu])
      plt.plot(rhd_choice.spectrum.waves, rhd_choice.spectrum.I[mu, :], label=xmulabel)
    
    plt.xlabel('wavelength [nm]')
    plt.ylabel('intensity [J m$^{-2}$ s$^{-1}$-1 Hz$^{-1}$ sr$^{-1}$]')
    
    plt.legend(title='mu')
    plt.show()
    
    mu = 2#ray of choice
    intensity_new = rhd_choice.spectrum.I[mu, :]*1.9e14 #conversion factor
    
    np.savez(filename,wl_rh=rhd_choice.spectrum.waves,int_rh=intensity_new)
    
    
