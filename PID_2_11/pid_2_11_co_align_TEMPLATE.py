#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2 December 2023
Author: Cole Tamburri, University of Colorado Boulder, National Solar 
Observatory, Laboratory for Atmospheric and Space Physics

Description of script: 
    Co-alignment routines for ViSP and VBI to SDO/HMI and SDO/AIA, applied to
    19 August 2022 20:42UT observations of GOES class C6.7 flare.

"""

#package import
import dkistpkg_ct as DKISTanalysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.visualization import ImageNormalize, SqrtStretch

import sunpy.coordinates
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

from astropy.io import fits

muted = DKISTanalysis.color_muted2()

#define paths
pid = 'pid_2_11'

#ViSP
path = '/Volumes/ViSP_External/pid_2_11/'
folder1 = 'AQWDJ'

#VBI
path_vbi = '/Volumes/VBI_Aug_15_Aug_25_22/pid_1_84/'
folder1_vbi = 'BXWNO'
folder2_vbi = 'BYMOL'

dir_list2 = DKISTanalysis.pathdef(path,folder1)

# Stonyhurst lon/lat from JHv
lon = 57.99 #degrees
lat = -15 #degrees

#central wl
wl = 396.8 # Ca II H

startstep=2000 #where does interesting bit begin?
endstep=4500#where does interesting bit end?

#spatial coordinates
hpc1_arcsec, hpc2_arcsec, x_center, y_center, z, rho, mu, \
    doppshnonrel, doppshrel = \
    DKISTanalysis.spatialinit(path,folder1,dir_list2,lon,lat,wl)

# calculate limb darkening
limbdarkening = DKISTanalysis.limbdarkening(wl,mu=mu,nm=True) # for Ca II H

#processing of raster
image_data_arr_arr, rasterpos, times = \
    DKISTanalysis.multistepprocess(path,folder1,dir_list2,startstep=startstep,div=1,endstep=endstep)

hdul1 = fits.open(path+folder1+'/'+dir_list2[0])

#initial scan from ViSP w/ spatial coordinates
spatial_range, dispersion_range = DKISTanalysis.spatialaxis(path,folder1,dir_list2,line='Ca II H',pid='2_11')

spatial_range2, raster_range, slitlen, rastersize= DKISTanalysis.vispranges(hdul1,
                                                        spatial_range,
                                                        nslitpos=91,
                                                        pid='2_11')

x_cent, y_cent, x_delt, y_delt, x_range, y_range, arcsec_slit, nspace = \
    DKISTanalysis.space_range(hdul1,pid='2_11')

#prep image from ViSP initial scan    
image_data_arrs0 = DKISTanalysis.imgprep(path,folder1,dir_list2,startstep,
                                         endstep,pid='2_11')

#plotting of ViSP scan
caiiavgs = DKISTanalysis.line_avg(image_data_arr_arr,600,750,1999,nspace)
hepavgs = DKISTanalysis.line_avg(image_data_arr_arr,750,915,1999,nspace)

DKISTanalysis.pltraster(caiiavgs,raster_range,spatial_range2)

#processing of VBI data
vbi_X, vbi_Y, hdul1_vbi, dat0_vbi = DKISTanalysis.vbi_process(path_vbi,
                                                              folder1_vbi)

X,Y = np.meshgrid(raster_range,spatial_range2)

#ID corresponding points in ViSP and VBI
aa = DKISTanalysis.plt_precoalign(vbi_X,vbi_Y,hdul1_vbi,X,Y,caiiavgs,
                                  matplotlib,dat0_vbi)

#ViSP to VBI
visp_X_new, visp_Y_new = DKISTanalysis.vbi_visp_transformation(aa,X,Y,nspace,4,
                                                               vbi_X,vbi_Y,
                                                               dat0_vbi,
                                                               caiiavgs,
                                                               matplotlib)

#VBI to SDO calibration
start_time = Time('2022-08-19T20:42:30', scale='utc', format='isot')

#timestamps for SDO observations
lowerx = 675
upperx = 725

lowery = -500
uppery = -425

bottom_left = SkyCoord(lowerx*u.arcsec, lowery*u.arcsec, obstime=start_time, 
                       observer="earth", frame="helioprojective")
top_right = SkyCoord(upperx*u.arcsec, uppery*u.arcsec, obstime=start_time, 
                     observer="earth", frame="helioprojective")
    
cutout = a.jsoc.Cutout(bottom_left, top_right=top_right, tracking=True)

#query SDO, determine features to use in transformation
#remember features and order of features in bb; may need to run many times
#to be sure of best choices with VBI
query, bb = DKISTanalysis.query_sdo(start_time, "cole.tamburri@colorado.edu", 
                                    cutout, matplotlib, lowerx,upperx,lowery,uppery,wavelength = 304,
                                    timesamp=2,passband='cont')

#load continuum VBI?
#load VBI (perhaps different bandpass for this co-alignment, e.g. blue continuum
# w/ HMI continuum)
vbi_X2, vbi_Y2, hdul1_vbi2, dat0_vbi2 = DKISTanalysis.vbi_process(path_vbi,
                                                                  folder2_vbi)

# click to same features as in bb (from SDO)
cc = DKISTanalysis.points_vbi(vbi_X2,vbi_Y2,dat0_vbi2,matplotlib)

# transformation function
vbi_X_new, vbi_Y_new, COB2, A2, A1 = DKISTanalysis.vbi_to_sdo(bb,cc,vbi_X,vbi_Y)

# same transformation, with ViSP (using ViSP ---> VBI coordinates from above)
visp_X_new2, visp_Y_new2 = DKISTanalysis.visp_sdo_trans(visp_X_new,visp_Y_new, 
                                                        COB2, A2, A1, 
                                                        nspace = 2544, nwave=4)

# plotting final co-alignment
DKISTanalysis.plt_final_coalign(vbi_X_new, vbi_Y_new, dat0_vbi2, 
                      visp_X_new2, visp_Y_new2, caiiavgs,
                      dat0_vbi)
    
    