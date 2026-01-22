#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 11:10:28 2026

@author: coletamburri

Function to co-align two ViSP slits since the L1 headers do not have the correct
info, and save to a file.  The final sxes are only relative to the coordinates
stored in the L1 file of arm1 - i.e. the final coordinates will be in the 
reference frame defined by the L1 header of arm1.  Determination of absolute 
coordinates will require co-alignment between ViSP, VBI, SDO.

Options to save (0 or 1) the relative coordinates to a file, and also to 
process (0 or 1) the ViSP data from the two arms in question.  If not, just skips
to the co-alignment step so that trial and error can be used without re-reading
every file every time.  

Variables addshift and multshift should be changed until the crosshairs in the
ViSP data line up - both will need adjustment.

"""

import numpy as np
import astropy
from astropy.io import fits
import matplotlib.pyplot as plt
import dkistpkg_ct as DKIST_analysis

base = '/Volumes/ViSP_External/pid_2_11/'
arm1_folder = 'ARYEE/'
arm2_folder = 'AQWDJ/'
startind = 2548
nslit = 91
arm1 = 'hbeta'
arm2 = 'caIIH'

addshift = 3.134 
multshift = 0.00033

#save to file?
save = 0

#process all variables? only required first time
#otherwise can just iterate to try different values of addshift and multshift
process = 1

if process == 1:
    arm1_processed =  '/Users/coletamburri/Desktop/11_Aug_2024_Cclass_Flare/'+\
        'Processed_ViSP_VBI_11Aug2024/ViSP_spectra_processed_11Aug24_hbeta.npz'
    arm2_processed = '/Users/coletamburri/Desktop/11_Aug_2024_Cclass_Flare/'+\
        'Processed_ViSP_VBI_11Aug2024/ViSP_spectra_processed_11Aug24_CaII.npz'
    
        
    xarr_arm1, yarr_arm1, xarr_arm2, yarr_arm2, dir_list_arm1, dir_list_arm2=\
        DKIST_analysis.prep_arms(base,arm1_folder,arm2_folder,arm1_processed,\
                            arm2_processed)
            
    Xhbeta,Yhbeta = np.meshgrid(xarr_arm1,yarr_arm1)
    XcaIIH,YcaIIH = np.meshgrid(xarr_arm2,yarr_arm2)
    
    
            
    image_arr_arm1, rasterpos_arm1, times_arm1 = \
        DKIST_analysis.multistepprocess(base,arm1_folder,dir_list_arm1,\
                                        startstep=startind,div=1,\
                                            endstep=startind+nslit+1)
            
    image_arr_arm2, rasterpos_arm2, times_arm2 = \
            DKIST_analysis.multistepprocess(base,arm2_folder,dir_list_arm2,\
                                            startstep=startind,div=1,\
                                                endstep=startind+nslit+1)
                
    mean_arm1 = np.mean(image_arr_arm1[:,500:,:],1)
    
    mean_arm2 = np.mean(image_arr_arm2[:,550:800,:],1)
     

# this is the major worker to change - both addition offset and multiplicative
shiftaxis = addshift-(np.arange(len(yarr_arm2[:-1]))*multshift) 
       
fig, ax = plt.subplots(1,4,dpi=100,figsize=(9,7))

ax.flatten()[0].pcolormesh(xarr_arm1,yarr_arm1,np.transpose(mean_arm1),cmap='sdoaia171')
ax.flatten()[1].pcolormesh(xarr_arm2,yarr_arm2[:-1]-shiftaxis,np.transpose(mean_arm2),cmap='grey')

ax.flatten()[2].pcolormesh(xarr_arm1,yarr_arm1,np.transpose(mean_arm1),cmap='sdoaia171')
ax.flatten()[3].pcolormesh(xarr_arm2,yarr_arm2[:-1]-shiftaxis,np.transpose(mean_arm2),cmap='grey')

ax.flatten()[0].set_ylim([-217,-215])
ax.flatten()[1].set_ylim([-217,-215])

ax.flatten()[2].set_ylim([-262,-259])
ax.flatten()[3].set_ylim([-262,-259])
fig.show()

fig, ax = plt.subplots(1,2, sharey=True,dpi=100,figsize=(3,7))
ax.flatten()[0].pcolormesh(xarr_arm1,yarr_arm1,np.transpose(mean_arm1),cmap='sdoaia171')
ax.flatten()[1].pcolormesh(xarr_arm2,yarr_arm2[:-1]-shiftaxis,np.transpose(mean_arm2),cmap='grey')


fig.subplots_adjust(wspace=0)
fig.align_ylabels()

yarr_arm2_final = yarr_arm2-shiftaxis

if save == 1:

    outfile = '/Users/coletamburri/Desktop/ViSPcoords'+arm1+'_'+arm2+'.npz'
    
    np.savez(outfile,xarr_arm1=xarr_arm1,yarr_arm1=yarr_arm1,\
             xarr_arm2=xarr_arm2,yarr_arm2=yarr_arm2_final,\
                 line_arm1=arm1,line_arm2=arm2)

