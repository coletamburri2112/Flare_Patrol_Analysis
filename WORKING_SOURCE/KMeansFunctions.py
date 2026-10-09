#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 10:50:48 2026

@author: coletamburri
"""
import numpy as np
import matplotlib.pyplot as plt

def normalize(data):
    normarr=(data-np.nanmin(data))/(np.nanmax(data)-np.nanmin(data)) 
    return normarr

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def find_30p_height(curve,find_nearest):
    p30int = np.max(curve)*.3
    
    lowind,lowval = find_nearest(curve[0:round(len(curve)/2)],p30int)
    highind,highval = find_nearest(curve[round(len(curve)/2):],p30int)
    highind=highind+round(len(curve)/2)
    dist = (highind+lowind)/2
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

def find_centmin(wl,curve,find_nearest):
    centrevpos = np.where(curve[65:-100]==np.nanmin(curve[65:-100]))
    wlcentrev = wl[centrevpos]
    
    fig,ax=plt.subplots()
    ax.plot(wl,curve[65:-100])
    return centrevpos, wlcentrev

def find_velocity(rest_wl,obs_wl):
    c= 299792458 # m/s
    
    return c*(obs_wl-rest_wl)/(obs_wl)

def veltrans(x,mu=1):
    return ((((x)/cent)-1)*c/1000)/mu

def veltrans2(x,mu2=1):
    return ((((x)/cent)-1)*c/1000)/mu2

def wltrans(x):
    return ((((x/(c/1000))+1)*cent)-cent)