#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 07:09:33 2022

@author: owner
"""
# Most recent version: 23 May 2022 by Cole Tamburri
# University of Colorado Boulder
# Advisors: Maria D. Kazachenko and Adam F. Kowalski


# Major working functions for analysis of SDO/EVE 304 Angstrom light curves,
# SDO/AIA 1600 Angstrom images and ribbon masks, and SDO/HMI magnetograms.
# Includes algorithm for determination of separation and elongation of both
# ribbons relative to the polarity inversion line. Flares able to beanalyzed
# are contained in the RibbonDB database (Kazachenko et al. 2017).  Prior to
# running this script, the user should obtain flare light curves and times
# corresponding to the modeled flares in this database, for which the
# impulsiveness index has been determined previously.

# The polarity inversion line is determined by convolving the HMI masks
# associated for each flare with a Gaussian of predetermined length, then
# finding the major region of intersection between these and using the
# resulting heatmap to fit a fourth-order polynomial.  Details of separation
# and elongation methods relative to the PIL are included below.

# Reconnection rates and ribbon areas for both positive and negative ribbons
# are determined, and the values corresponding to the rise phase of the flare
# (with some flare-specific variation) are fit to an exponential model, in
# order to prepare for modeling efforts of flare heating particularly in the
# rise phase.

# Separation and elongation values (perpendicular and parallel PIL-relative
# motion, respectively) are used to find separation and elongation rates,
# through which significant periods of these two phases of ribbon motion can
# be identified.

# Plotting and data presentation routines are also below, which includes an
# animation showing the timing of separation, elongation, and chromospheric
# line light curves.

# Addition of shear quantification code, 29 April 2022

# Addition of four-paneled figure comparing HXR, AIA, EVE, and shear
# quantification values with coordinated timestamps included.

from os.path import dirname, join as pjoin
import scipy.io as sio
from scipy.io import readsav
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animat
import datetime
from scipy.spatial.distance import cdist
import scipy.signal
import matplotlib.dates as mdates
from astropy.convolution import convolve, Gaussian2DKernel
import time as timepkg
from matplotlib import font_manager
from matplotlib.ticker import MaxNLocator
import pandas as pd
from matplotlib import gridspec
import latex
import netCDF4 as nc
from matplotlib import font_manager


def conv_facts():
    """
    Conversion factors for images.

    Returns
    -------
    X : list
        Meshgrid of x values for image coordinates.
    Y : list
        Meshgrid of y values for image coordinates.
    conv_f : float
        Conversion factor between pixels and megameters.
    xarr_Mm : list
        x-coordinates, in megameters.
    yarr_Mm : list
        y-coordinates, in megameters.

    """
    pix_to_arcsec = 0.6  # asec/pix
    arcsec_to_radians = 1 / 206265  # rad/asec
    radians_to_Mm = 149598  # Mm/rad

    conv_f = pix_to_arcsec * arcsec_to_radians * radians_to_Mm  # Mm/pix

    xarr_Mm = np.zeros(800)
    yarr_Mm = np.zeros(800)

    for i in range(800):
        xarr_Mm[i] = (i-400)*conv_f
        yarr_Mm[i] = (i-400)*conv_f

    X, Y = np.meshgrid(xarr_Mm, yarr_Mm)

    return X, Y, conv_f, xarr_Mm, yarr_Mm


def exponential(x, a, b):
    """
    Defines exponential function.

    Parameters
    ----------
    x : float
        Input x value for function.
    a : float
        Amplitude of exponential function.
    b : float
        Second parameter of exponential function.

    Returns
    -------
    float
        Output of exponential function.

    """
    return a * np.exp(b * x)


def exponential_neg(x, a, b):
    """
    Negative amplitude exponential function.

    Parameters
    ----------
    x : float
        Input x value for function.
    a : float
        Amplitude of exponential function.
    b : float
        Second parameter of exponential function.

    Returns
    -------
    float
        Output of exponential function.

    """
    return -a * np.exp(b * x)


def curve_length(curve):
    """
    Sum of Euclidean distances between points
    """
    return np.sum(np.sqrt(np.sum((curve[:-1] - curve[1:])**2, axis=1)))


def datenum_to_datetime(datenum):
    """
    Convert Matlab datenum into Python datetime.

    Parameters
    ----------
    datenum : float
        Datenum value.

    Returns
    -------
    ret : datetime
        Converted datetime value.

    """
    days = datenum % 1
    ret = datetime.datetime.fromordinal(int(datenum)) + \
        datetime.timedelta(days=days) - datetime.timedelta(days=366)

    return ret


def datenum(d):
    """
    Convert from ordinal to datenum.
    """
    return 366 + d.toordinal() + (d - datetime.datetime.
                                  fromordinal(d.toordinal())).\
        total_seconds()/(24*60*60)


def find_nearest(array, value):
    """
    Find nearest value in array to a value.

    Parameters
    ----------
    array : list
        Array of values to search through.
    value : float
        Value to find the nearest element in array closest to.

    Returns
    -------
    float
        Nearest value in array to "value"

    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def format_time():
    """
    Time formatter.

    Returns
    -------
    string
        Formating for times.

    """
    t = datetime.datetime.now()
    s = t.strftime('%Y-%m-%d %H:%M:%S.%f')
    return s[:-3]


def find_nearest_ind(array, value):
    """
    Find index of element in array closest to value.

    Parameters
    ----------
    array : list
        Array of values to search through.
    value : float
        Value to find the nearest element in array closest to.

    Returns
    -------
    idx: int
        Index of nearest value in array to "value"
    """

    array = np.asarray(array)
    idx = np.nanargmin(np.abs(array - value))
    return idx


def load_variables(bestflarefile, year, mo, day, sthr, stmin, arnum, xclnum,
                   xcl,flag=0):
    """
    Load variables from HMI and AIA files.

    Parameters
    ----------
    bestflarefile : string
        Path to file containing information about the best-performing flares.
    year : int
        Year of event.
    mo : int
        Month of event.
    day : int
        Day of event.
    sthr : int
        Hour of event start
    stmin : int
        Minute of event start.
    arnum : int
        Active region number.
    xclnum : int
        X-ray class number.
    xcl : str
        X-ray class.

    Returns
    -------
    sav_data_aia : AttrDict
        Dictionary containing all of the saved parameters in the AIA file.
    sav_data : AttrDict
        Dictionary containing all of the saved parameters in the HMI file.
    best304 : dict
        Dictionary containing the SDO/EVE 304 Angstrom data of the
        best-performing flares in ribbonDB.
    start304 : list
        Array containing the start times for the flares in best304.
    peak304 : list
        Array containing the peak times for the flares in best304.
    end304 : list
        Array containing the end times for the flares in best304.
    eventindices : list
        Indices of best flares in best304.
    times304 : list
        Time points for all flares in best304.
    curves304 : list
        Light curves for all flares in best304.
    aia_cumul8 : list
        Cumulative ribbon masks from AIA.
    aia_step8 : list
        Instantaneous ribbon masks from AIA
    last_cumul8 : list
        The last image in the cumulative mask array.
    hmi_dat : list
        HMI image prior to the flare, assumed to be the same configuration
        throughout the flare.
    last_mask : list
        The last ribbon mask, multiplied by the HMI image for polarity.
    flag: if diff posfile (for high2 processing, Jan 2024)

    """
    data_dir = pjoin(dirname(sio.__file__), 'tests', 'data')

    # Load matlab file, get 304 light curves and start/peak/end times for flare
    best304 = sio.loadmat(bestflarefile)

    start304 = best304['starttimes_corr'][:, 0]
    peak304 = best304['maxtimes_corr'][:, 0]
    end304 = best304['endtimes_corr'][:, 0]
    eventindices = best304['starttimes_corr'][:, 1]
    times304 = best304['event_times']
    curves304 = best304['event_curves']
    sav_fname_aia = pjoin(data_dir, "/Users/coletamburri/Desktop/Impulsiveness_Paper/Imp_Study_Transfer/"
                          "AIA_Files/aia1600blos" + str(year).zfill(4) +
                          str(mo).zfill(2) + str(day).zfill(2) + "_" +
                          str(sthr).zfill(2) + str(stmin).zfill(2) + "_" +
                          str(arnum).zfill(5) + "_"+xcl + str(xclnum) + ".sav")
    sav_data_aia = readsav(sav_fname_aia)
    sav_fname = ("/Users/coletamburri/Desktop/Impulsiveness_Paper/Imp_Study_Transfer/HMI_files/posfile" +
                 str(year).zfill(4) + str(mo).zfill(2) + str(day).zfill(2) +
                 "_" + str(sthr).zfill(2) + str(stmin).zfill(2) + "_" +
                 str(arnum).zfill(5) + "_"+xcl + str(xclnum) +
                 "_cut08_sat5000.00_brad.sav")
    
    if flag == 1:
        print()
        sav_fname = ("/Users/coletamburri/Desktop/data4cole/posfile" +
                     str(year).zfill(4) + str(mo).zfill(2) + str(day).zfill(2) +
                     "_" + str(sthr).zfill(2) + str(stmin).zfill(2) + "_" +
                     str(arnum).zfill(5) + "_"+xcl + str(xclnum) +
                     "_cut08_sat5000.00_blos.sav")
    
    sav_data = readsav(sav_fname)

    aia_cumul8 = sav_data.pos8
    last_cumul8 = aia_cumul8[-1, :, :]  # Last frame
    hmi_dat = sav_data.hmi
    last_mask = last_cumul8 * hmi_dat

    aia_step8 = sav_data.inst_pos8

    return sav_data_aia, sav_data, best304, start304, peak304, end304, \
        eventindices, times304, curves304, aia_cumul8, aia_step8, \
        last_cumul8, hmi_dat, last_mask


def pos_neg_masking(aia_cumul8, aia_step8, hmi_dat, last_mask):
    """
    Masking of positive and negative ribbons according to HMI polarity.

    Parameters
    ----------
    aia_cumul8 : list
        Cumulative ribbon masks.
    aia_step8 : list
        Instantaneous ribbon masks.
    hmi_dat : list
        HMI image prior to the flare, assumed to be the same configuration
        throughout the flare.
    last_mask : list
        The last ribbon mask, multiplied by the HMI image for polarity.

    Returns
    -------
    hmi_cumul_mask1 : list
        Cumulative magnetic field strength masking estimates for all flare
        images.
    hmi_step_mask1 : list
        Instantaneous magnetic field strength masking estimates for all flare
        images.
    hmi_pos_mask_c : list
        Single-frame mask for negative HMI magnetic field, populated with 1.
    hmi_neg_mask_c : list
        Single-frame mask for negative HMI magnetic field, populated with -1.

    """
    hmi_cumul_mask = np.zeros(np.shape(aia_cumul8))
    hmi_cumul_mask1 = np.zeros(np.shape(aia_cumul8))

    # Find HMI masks for all frames - cumulative
    for i in range(len(aia_cumul8)):
        frame = np.squeeze(aia_cumul8[i, :, :])
        hmi_cumul_mask[i, :, :] = frame * hmi_dat

    # Convert to positive and negative polarities
    for i in range(len(hmi_cumul_mask)):
        for j in range(len(hmi_cumul_mask[0])):
            for k in range(len(hmi_cumul_mask[1])):
                if hmi_cumul_mask[i, j, k] > 0:
                    hmi_cumul_mask1[i, j, k] = 1
                elif hmi_cumul_mask[i, j, k] < 0:
                    hmi_cumul_mask1[i, j, k] = -1
                else:
                    hmi_cumul_mask1[i, j, k] = 0

    hmi_step_mask = np.zeros(np.shape(aia_step8))
    hmi_step_mask1 = np.zeros(np.shape(aia_step8))

    # Find HMI masks for all frames - instantaneous
    for i in range(len(aia_step8)):
        frame = np.squeeze(aia_step8[i, :, :])
        hmi_step_mask[i, :, :] = frame * hmi_dat

    # Convert to positive and negative polarities
    for i in range(len(hmi_step_mask)):
        for j in range(len(hmi_step_mask[0])):
            for k in range(len(hmi_step_mask[1])):
                if hmi_step_mask[i, j, k] > 0:
                    hmi_step_mask1[i, j, k] = 1
                elif hmi_step_mask[i, j, k] < 0:
                    hmi_step_mask1[i, j, k] = -1
                else:
                    hmi_step_mask1[i, j, k] = 0

    # Single-frame masks for positive and negative ribbons
    hmi_pos_mask_c = np.zeros(np.shape(hmi_dat))
    hmi_neg_mask_c = np.zeros(np.shape(hmi_dat))

    for i in range(len(hmi_dat)):
        for j in range(len(hmi_dat[0])):
            if last_mask[i, j] > 0:
                hmi_pos_mask_c[i, j] = 1
                hmi_neg_mask_c[i, j] = 0
            elif last_mask[i, j] < 0:
                hmi_pos_mask_c[i, j] = 0
                hmi_neg_mask_c[i, j] = -1
            else:
                hmi_pos_mask_c[i, j] = 0
                hmi_neg_mask_c[i, j] = 0

    return hmi_cumul_mask1, hmi_step_mask1, hmi_pos_mask_c, hmi_neg_mask_c


def spur_removal_sep(hmi_neg_mask_c, hmi_pos_mask_c, pos_crit=3,
                     neg_crit=3, pt_range=[-2, -1, 1, 2], ihi=800, ilo=0,
                     jhi=800, jlo=0, ihi2=800, ilo2=0, jhi2=800,
                     jlo2=0, ihi3=800, jlo3=0):
    """
    Spur removal in ribbon masks for the perpendicular motion identification.
    Removes regions where both negative and positive pixels exist.

    Parameters
    ----------
    hmi_neg_mask_c : list
        Single-frame mask for negative HMI magnetic field, populated with -1.
    hmi_pos_mask_c : list
        Single-frame mask for negative HMI magnetic field, populated with 1.
    pos_crit : int, optional
        Number of positive points surrounding a negative pixel for which the
        negative pixel should be removed. The default is 3.
    neg_crit : int, optional
        Number of positive points surrounding a negative pixel for which the
        negative pixel should be removed. The default is 3.
    pt_range : list, optional
        Pixels to search around each pixel for opposite polarity. The default
        is [-2,-1,1,2].
    ihi : int, optional
        Upper i-limit for allowance of pixel masks, negative. The default is
        800.
    ilo : int, optional
        Lower i-limit for allowance of pixel masks, negative. The default is
        0.
    jhi : int, optional
        Upper j-limit for allowance of pixel masks, negative. The default is
        800.
    jlo : int, optional
        Lower j-limit for allowance of pixel masks, negative. The default is
        0.
    ihi2 : int, optional
        Upper i-limit for allowance of pixel masks, positive. The default is
        800.
    ilo2 : int, optional
        Lower i-limit for allowance of pixel masks, positive. The default is
        0.
    jhi2 : int, optional
        Upper j-limit for allowance of pixel masks, positive. The default is
        800.
    jlo2 : int, optional
        Lower j-limit for allowance of pixel masks, positive. The default is
        0.
    ihi3 : int, optional
        Special limit for highest impulsiveness flare.  The default is 800
    jlo3 : int, optional
        Special limit for highest impulsiveness flare.  The default is 0.

    Returns
    -------
    neg_rem : list
        The negative polarity HMI image, with spurs removed.
    pos_rem : list
        The positive polarity HMI image, with spurs removed.

    """
    neg_rem = np.zeros(np.shape(hmi_neg_mask_c))
    pos_rem = np.zeros(np.shape(hmi_pos_mask_c))

    # If > neg_crit positive pixels surround a negative pixel, remove negative
    # pixel.
    for i in range(len(neg_rem) - 2):
        for j in range(len(neg_rem[0]) - 2):
            n = 0
            if hmi_neg_mask_c[i, j] == -1:
                for k in pt_range:
                    for h in pt_range:
                        if hmi_pos_mask_c[i + k, j - h] == 1:
                            n = n + 1
                if n > neg_crit or i > ihi or i < ilo or j < jlo or j > jhi\
                        or (i > ihi3 and j < jlo3):
                    neg_rem[i, j] = 0
                else:
                    neg_rem[i, j] = -1
            else:
                neg_rem[i, j] = 0

    # If > pos_crit negative pixels surround a positive pixel, remove positive
    # pixel.
    for i in range(len(pos_rem) - 2):
        for j in range(len(pos_rem[0]) - 2):
            n = 0
            if hmi_pos_mask_c[i, j] == 1:
                for k in pt_range:
                    for h in pt_range:
                        if hmi_neg_mask_c[i + k, j - h] == -1:
                            n = n + 1
                if n > pos_crit or j > jhi2 or j < jlo2 or i < ilo2 or\
                        i > ihi2:
                    pos_rem[i, j] = 0
                else:
                    pos_rem[i, j] = 1
            else:
                pos_rem[i,  j] = 0

    return neg_rem, pos_rem


def gauss_conv(pos_rem, neg_rem, sigma=10):
    """
    Convolve HMI images with a Gaussian of specified width.

    Parameters
    ----------
    neg_rem : list
        The negative polarity HMI image, with spurs removed.
    pos_rem : list
        The positive polarity HMI image, with spurs removed.
    sigma : int, optional
        Width of the Gaussian to convolve with images. The default is 10.

    Returns
    -------
    hmi_con_pos_c : list
        Positive HMI, convolved with Gaussian.
    hmi_con_neg_c : list
        Negative HMI, convolved with Gaussian.
    pil_mask_c : list
        PIL mask found by multiplying positive and negative polarity PIL masks.

    """
    gauss_kernel = Gaussian2DKernel(sigma)
    hmi_con_pos_c = convolve(pos_rem, gauss_kernel)
    hmi_con_neg_c = convolve(neg_rem, gauss_kernel)

    # PIL mask is found by intersection of negative and positive HMI masks
    pil_mask_c = hmi_con_pos_c * hmi_con_neg_c

    return hmi_con_pos_c, hmi_con_neg_c, pil_mask_c


def pil_gen(pil_mask_c, hmi_dat, threshperc=0.05, lx=800, ly=800,
            polyor=4):
    """
    Generate PIL polynomial.

    Parameters
    ----------
    pil_mask_c : list
        PIL mask.
    hmi_dat : list
        Array of HMI values associated with the flare.
    threshperc: float, optional
        Percentage of maximum PIL mask value to allow into the polynomial fit.
        The default is 0.05.
    lx : int, optional
        Length of array in x direction. The default is 800.
    ly : list, optional
        Length of array in y direction. The default is 800.
    polyor : int, optional
        Order of fitting polynomial.  The default is 4.

    Returns
    -------
    pil_mask_c : list
        PIL mask.
    ivs : list
        x-values for PIL polynomial.
    dvs : list
        y-values for PIL polynomial.
    hmik : list
        HMI image, divided by 1000 for unit conversion.

    """

    # Make PIL mask positive
    pil_mask_c = -1.0 * pil_mask_c

    # Threshold for fitting of PIL polynomial
    thresh = threshperc * np.amax(pil_mask_c)

    # Isolate pixels certainly within the mask
    xc, yc = np.where(pil_mask_c > thresh)

    # Fitting of fourth-order polynomial to chosen pixels and generation of
    # PIL polynomial arrays
    x = np.linspace(0, lx, lx)
    y = np.linspace(0, ly, ly)
    coeffs = np.polyfit(y[yc], x[xc], polyor)

    ivs = y[yc]

    dvs = 0

    for i in range(len(coeffs)):
        dvs += coeffs[i] * ivs**(polyor - i)

    hmik = hmi_dat/1000

    return pil_mask_c, ivs, dvs, hmik


def mask_sep(aia_step8, hmi_dat):
    """
    Masking of each image for each time step, for use in separation value
    determination.

    Parameters
    ----------
    aia_step8 : list
        Instantaneous AIA ribbon masks, c=8.
    hmi_dat : list
        SDO/HMI magnetic field data for flare.

    Returns
    -------
    aia8_pos : list
        Contains only the positive ribbon masks for each time step.
    aia8_neg : list
        Contains only the negative ribbon masks for each time step.

    """

    aia8 = aia_step8
    aia8_pos = np.zeros(np.shape(aia8))
    aia8_neg = np.zeros(np.shape(aia8))

    # Separate positive and negative ribbons into different arrays
    for i in range(len(aia8)):
        for j in range(len(aia8[0])):
            for k in range(len(aia8[1])):
                if aia8[i, j, k] == 1 and hmi_dat[j, k] > 0:
                    aia8_pos[i, j, k] = 1
                elif aia8[i, j, k] == 1 and hmi_dat[j, k] < 0:
                    aia8_neg[i, j, k] = 1

    return aia8_pos, aia8_neg


def spur_removal_sep2(aia8_pos, aia8_neg, pos_crit=3, neg_crit=3,
                      pt_range=[-2, -1, 1, 2], jhi=800, jlo=0, khi=800,
                      klo=0, jhi2=800, jlo2=0, khi2=800, klo2=0):
    """
    Second step in removal of spurs from mask images for separation code. Limit
    window where ribbons are considered for PIL-relative perpendicular motion.

    Parameters
    ----------
    aia8_pos : list
        Output of mask_sep, containing positive mask isolated.
    aia8_neg : list
        Output of mask_sep, containing negative mask isolated
    pos_crit : int, optional
        Number of points surrounding another which will be allowed in the
        positive ribbon. The default is 3.
    neg_crit : int, optional
        Number of points surrounding another which will be allowed in the
        negative ribbon. The default is 3.
    pt_range : list, optional
        Range of points around which to search for other pixels of the same
        polarity. The default is [-2,-1,1,2].
    jhi : int, optional
        Upper j-limit for allowance of pixel masks, negative. The default is
        800.
    jlo : int, optional
        Lower j-limit for allowance of pixel masks, negative. The default is
        0.
    khi : int, optional
        Upper k-limit for allowance of pixel masks, negative. The default is
        800.
    klo : int, optional
        Lower k-limit for allowance of pixel masks, negative. The default is
        0.
    jhi2 : int, optional
        Upper j-limit for allowance of pixel masks, positive. The default is
        800.
    jlo2 : int, optional
        Lower j-limit for allowance of pixel masks, positive. The default is
        0.
    khi2 : int, optional
        Upper k-limit for allowance of pixel masks, positive. The default is
        800.
    klo2 : int, optional
        Lower k-limit for allowance of pixel masks, positive. The default is
        0.

    Returns
    -------
    pos_rem0 : list
        Masks with spurious pixels removed, positive ribbon.
    neg_rem0 : list
        Masks with spurious pixels removed, negative ribbon.

    """

    neg_rem0 = np.zeros(np.shape(aia8_pos))
    pos_rem0 = np.zeros(np.shape(aia8_neg))

    for i in range(len(neg_rem0)):
        for j in range(len(neg_rem0[0]) - 2):
            for k in range(len(neg_rem0[1]) - 2):
                n = 0
                if aia8_neg[i, j, k] == 1:
                    for h in pt_range:
                        for m in pt_range:
                            if aia8_neg[i, j + h, k + m] == 1:
                                n = n + 1
                    if n > neg_crit and j < jhi and j > jlo and k > klo \
                            and k < khi:
                        neg_rem0[i, j, k] = 1
                    else:
                        neg_rem0[i, j, k] = 0
                else:
                    neg_rem0[i, j, k] = 0

    for i in range(len(pos_rem0)):
        for j in range(len(pos_rem0[0]) - 2):
            for k in range(len(pos_rem0[1]) - 2):
                n = 0
                if aia8_pos[i, j, k] == 1:
                    for h in pt_range:
                        for m in pt_range:
                            if aia8_pos[i, j + h, k + m] == 1:
                                n = n + 1
                    if (n > pos_crit) and k < khi and k > klo and j > jlo and \
                            j < jhi:
                        pos_rem0[i, j, k] = 1
                    else:
                        pos_rem0[i, j, k] = 0
                else:
                    pos_rem0[i, j, k] = 0

    return pos_rem0, neg_rem0


def spur_removal_sepopt3(aia8_pos, aia8_neg, pos_crit=3, neg_crit=3,
                         pt_range=[-2, -1, 1, 2], jhi=800, jlo=0, khi=800,
                         klo=0, jhi2=800, jlo2=0, khi2=800, klo2=0):
    """
    Second step in removal of spurs from mask images for separation code -
    option for highest impulsiveness flare. Limit window where ribbons are
    considered for PIL-relative perpendicular motion.

    Parameters
    ----------
    aia8_pos : list
        Output of mask_sep, containing positive mask isolated.
    aia8_neg : list
        Output of mask_sep, containing negative mask isolated
    pos_crit : int, optional
        Number of points surrounding another which will be allowed in the
        positive ribbon. The default is 3.
    neg_crit : int, optional
        Number of points surrounding another which will be allowed in the
        negative ribbon. The default is 3.
    pt_range : list, optional
        Range of points around which to search for other pixels of the same
        polarity. The default is [-2,-1,1,2].
    jhi : int, optional
        Upper j-limit for allowance of pixel masks, negative. The default is
        800.
    jlo : int, optional
        Lower j-limit for allowance of pixel masks, negative. The default is
        0.
    khi : int, optional
        Upper k-limit for allowance of pixel masks, negative. The default is
        800.
    klo : int, optional
        Lower k-limit for allowance of pixel masks, negative. The default is
        0.
    jhi2 : int, optional
        Upper j-limit for allowance of pixel masks, positive. The default is
        800.
    jlo2 : int, optional
        Lower j-limit for allowance of pixel masks, positive. The default is
        0.
    khi2 : int, optional
        Upper k-limit for allowance of pixel masks, positive. The default is
        800.
    klo2 : int, optional
        Lower k-limit for allowance of pixel masks, positive. The default is
        0.

    Returns
    -------
    pos_rem0 : list
        Masks with spurious pixels removed, positive ribbon.
    neg_rem0 : list
        Masks with spurious pixels removed, negative ribbon.

    """
    neg_rem0 = np.zeros(np.shape(aia8_pos))
    pos_rem0 = np.zeros(np.shape(aia8_neg))

    for i in range(len(neg_rem0)):
        for j in range(len(neg_rem0[0]) - 2):
            for k in range(len(neg_rem0[1]) - 2):
                n = 0
                if aia8_neg[i, j, k] == 1:
                    for h in pt_range:
                        for m in pt_range:
                            if aia8_neg[i, j + h, k + m] == 1:
                                n = n + 1
                    if (n > neg_crit) and (j < jhi and j > jlo and k > klo and
                                           k < khi):
                        neg_rem0[i, j, k] = 1
                    else:
                        neg_rem0[i, j, k] = 0
                    if (j > 400 and k > 400 and k < 425):
                        neg_rem0[i, j, k] = 0
                else:
                    neg_rem0[i, j, k] = 0

    for i in range(len(pos_rem0)):
        for j in range(len(pos_rem0[0]) - 2):
            for k in range(len(pos_rem0[1]) - 2):
                n = 0
                if aia8_pos[i, j, k] == 1:
                    for h in pt_range:
                        for m in pt_range:
                            if aia8_pos[i, j + h, k + m] == 1:
                                n = n + 1
                    if n > pos_crit and k < khi and k > klo and j > jlo and\
                            j < jhi:
                        pos_rem0[i, j, k] = 1
                    else:
                        pos_rem0[i, j, k] = 0
                else:
                    pos_rem0[i, j, k] = 0

    return pos_rem0, neg_rem0


def separation(aia_step8, ivs, dvs, pos_rem0, neg_rem0):
    """
    Algorithm for determination of parallel motion for positive and negative
    ribbons.

    Parameters
    ----------
    aia_step8 : list
        Instantaneous AIA ribbon masks, c=8.
    ivs : list
        x-values for PIL polynomial.
    dvs : list
        y-values for PIL polynomial.
    aia8_pos : list
        Contains only the positive ribbon masks for each time step.
    aia8_neg : list
        Contains only the negative ribbon masks for each time step.

    Returns
    -------
    distpos_med : list
        Parallel distance of positive ribbon from PIL, median of all pixel
        distances.
    distpos_mean : list
        Parallel distance of positive ribbon from PIL, mean of all pixel
        distances.
    distneg_med : list
        Parallel distance of negative ribbon from PIL, median of all pixel
        distances.
    distpos_mean : list
        Parallel distance of negative ribbon from PIL, mean of all pixel
        distances.

    """

    # Create array of PIL mask values
    pil = list(zip(ivs, dvs))

    distpos_med = np.zeros(len(aia_step8))
    distneg_med = np.zeros(len(aia_step8))
    distpos_mean = np.zeros(len(aia_step8))
    distneg_mean = np.zeros(len(aia_step8))

    # Main working function for separation
    for i in range(len(aia_step8)):
        posframe = pos_rem0[i, :, :]
        negframe = neg_rem0[i, :, :]
        xpos, ypos = np.where(posframe == 1)
        xneg, yneg = np.where(negframe == 1)
        pos_ops = list(zip(ypos, xpos))
        neg_ops = list(zip(yneg, xneg))
        if len(pos_ops) > 0:
            # Distance from positive pixels to each of the PIL pixels
            allpos = cdist(pos_ops, pil)
            # Set the minimum for each pixel first
            allpos_min = np.amin(allpos, axis=1)
            # Median and mean of distances
            distpos_med[i] = np.median(allpos_min)
            distpos_mean[i] = np.mean(allpos_min)
        if len(neg_ops) > 0:
            # Same as in positive pixels
            allneg = cdist(neg_ops, pil)
            allneg_min = np.amin(allneg, axis=1)
            distneg_med[i] = np.median(allneg_min)
            distneg_mean[i] = np.mean(allneg_min)

    return distpos_med, distpos_mean, distneg_med, distpos_mean


def mask_elon(aia_cumul8, hmi_dat):
    """
    Masking for elongation algorithm.

    Parameters
    ----------
    aia_cumul8 : list
        Cumulative ribbon masks, c=8.
    hmi_dat : list
        SDO/HMI image data for flare.

    Returns
    -------
    aia8_pos_2 : list
        Contains only the positive cumulative ribbon masks for each time step.
    aia8_neg_2 : list
        Contains only the negative cumulative ribbon masks for each time step.

    """

    aia8_a = aia_cumul8
    aia8_pos_2 = np.zeros(np.shape(aia8_a))
    aia8_neg_2 = np.zeros(np.shape(aia8_a))

    # Separation of cumulative ribbon masks into separate arrays for opposite
    # polarity
    for i, j, k in np.ndindex(aia8_a.shape):
        if aia8_a[i, j, k] == 1 and hmi_dat[j, k] > 0:
            aia8_pos_2[i, j, k] = 1
        elif aia8_a[i, j, k] == 1 and hmi_dat[j, k] < 0:
            aia8_neg_2[i, j, k] = 1

    return aia8_pos_2, aia8_neg_2


def spur_removal_elon(aia8_pos_2, aia8_neg_2, pos_crit=3, neg_crit=3,
                      pt_range=[-2, -1, 1, 2], jhi=800, jlo=0, khi=800,
                      klo=0, jhi2=800, jlo2=0, khi2=800, klo2=0):
    """
    Removal of isolated regions of very few pixels in all time step images.

    Parameters
    ----------
    aia8_pos_2 : list
        Contains only the positive cumulative ribbon masks for each time step.
    aia8_neg_2 : list
        Contains only the negative cumulative ribbon masks for each time step.
    pos_crit : list, optional
        The number of pixels in positive ribbon within a region above which the
        point is allowed to remain in the image. The default is 3.
    neg_crit : list, optional
        The number of pixels in negative ribbon within a region above which the
        point is allowed to remain in the image. The default is 3.
    pt_range : list, optional
        Pixels to search around each pixel for the same polarity. The default
        is [-2,-1,1,2].
    jhi : int, optional
        Upper j-limit for allowance of pixel masks, negative. The default is
        800.
    jlo : int, optional
        Lower j-limit for allowance of pixel masks, negative. The default is
        0.
    khi : int, optional
        Upper k-limit for allowance of pixel masks, negative. The default is
        800.
    klo : int, optional
        Lower k-limit for allowance of pixel masks, negative. The default is
        0.
    jhi2 : int, optional
        Upper j-limit for allowance of pixel masks, positive. The default is
        800.
    jlo2 : int, optional
        Lower j-limit for allowance of pixel masks, positive. The default is
        0.
    khi2 : int, optional
        Upper k-limit for allowance of pixel masks, positive. The default is
        800.
    klo2 : int, optional
        Lower k-limit for allowance of pixel masks, positive. The default is
        0.

    Returns
    -------
    neg_rem1 : list
        Vetted positive ribbon with the above criteria for each pixel.
    pos_rem1 : list
        Vetted negative ribbon with the above criteria for each pixel.

    """
    neg_rem1 = np.zeros(np.shape(aia8_pos_2))
    pos_rem1 = np.zeros(np.shape(aia8_neg_2))

    # If neg_crit number of pixels not exceeded in a certain region, remove
    # central pixel - for negative mask, then positive mask, at each time step
    for i in range(len(neg_rem1)):
        for j in range(len(neg_rem1[0]) - 2):
            for k in range(len(neg_rem1[1]) - 2):
                n = 0
                if aia8_neg_2[i, j, k] == 1:
                    for h in pt_range:
                        for m in pt_range:
                            if aia8_neg_2[i, j + h, k + m] == 1:
                                n = n + 1
                    if (n > neg_crit) and k > klo and k < khi and j > jlo and\
                            j < jhi:
                        neg_rem1[i, j, k] = 1
                    else:
                        neg_rem1[i, j, k] = 0
                else:
                    neg_rem1[i, j, k] = 0

    for i in range(len(pos_rem1)):
        for j in range(len(pos_rem1[0]) - 2):
            for k in range(len(pos_rem1[1]) - 2):
                n = 0
                if aia8_pos_2[i, j, k] == 1:
                    for h in pt_range:
                        for m in pt_range:
                            if aia8_pos_2[i, j + h, k + m] == 1:
                                n = n + 1
                    if n > pos_crit and k > klo2 and k < khi2 and j > jlo2 and\
                            j < jhi2:
                        pos_rem1[i, j, k] = 1
                    else:
                        pos_rem1[i, j, k] = 0
                else:
                    pos_rem1[i, j, k] = 0

    return neg_rem1, pos_rem1


def lim_pil(ivs, dvs):
    """
    Limt of the inversion line within a certain number of pixels from the
    median image value.

    Parameters
    ----------
    ivs : list
        x-values for PIL polynomial.
    dvs : list
        y-values for PIL polynomial.

    Returns
    -------
    ivs_lim : list
        Vetted x-values for PIL polynomial.
    dvs_lim : list
        Vetted y-values for PIL polynomial.
    med_x : int
        Median x pixel in image.
    med_y : int
        Median y pixel in image.

    """
    med_x = np.median(ivs)
    med_y = np.median(dvs)

    ivs_lim = []
    dvs_lim = []

    for i in range(len(ivs)):
        if not (ivs[i] < (med_x - 200)) and not (ivs[i] > (med_x + 200)):
            ivs_lim.append(ivs[i])
            dvs_lim.append(dvs[i])

    return ivs_lim, dvs_lim, med_x, med_y


def rib_lim_elon(aia8_pos_2, aia8_neg_2, pos_rem1, neg_rem1, med_x, med_y,
                 ylim0_pos, ylim1_pos, ylim0_neg, ylim1_neg, xlim0_pos,
                 xlim1_pos, xlim0_neg, xlim1_neg):
    """
    Limiting of ribbons for processing with elongation algorithm.

    Parameters
    ----------
    aia8_pos_2 : list
        Contains only the positive cumulative ribbon masks for each time step.
    aia8_neg_2 : list
        Contains only the negative cumulative ribbon masks for each time step.
    neg_rem1 : list
        Vetted positive ribbon with the above criteria for each pixel.
    pos_rem1 : list
        Vetted negative ribbon with the above criteria for each pixel.
    med_x : int
        Median x pixel in image.
    med_y : int
        Median y pixel in image.
    ylim0_pos : int
        Lower y-limit for positive ribbon.
    ylim1_pos : int
        Upper y-limit for positive ribbon
    ylim0_neg : int
        Lower y-limit for negative ribbon
    ylim1_neg : int
        Upper y-limit for negative ribbon.
    xlim0_pos : int
        Lower x-limit for positive ribbon
    xlim1_pos : int
        Upper x-limit for positive ribbon
    xlim0_neg : int
        Lower x-limit for negative ribbon
    xlim1_neg : int
        Upper x-limit for negative ribbon

    Returns
    -------
    aia_pos_rem : list
        Isolated positive ribbon masks.
    aia_neg_rem : list
        Isolated negative ribbon masks.

    """
    aia_pos_rem = np.zeros(np.shape(aia8_pos_2))
    aia_neg_rem = np.zeros(np.shape(aia8_neg_2))

    # Limit the negative ribbon image to a certain region in the image
    for i in range(len(aia8_neg_2)):
        for j in range(ylim0_neg, ylim1_neg):
            for k in range(xlim0_neg, xlim1_neg):
                if neg_rem1[i, j, k] > 0:
                    aia_neg_rem[i, j, k] = 1

    # Limit the positive ribbon image to a certain region in the image
    for i in range(len(aia8_pos_2)):
        for j in range(ylim0_pos, ylim1_pos):
            for k in range(xlim0_pos, xlim1_pos):
                if pos_rem1[i, j, k] > 0:
                    aia_pos_rem[i, j, k] = 1

    return aia_pos_rem, aia_neg_rem


def split_rib(aia_pos_rem, aia_neg_rem, split_pos, split_neg):
    rib_pos_1 = np.zeros(np.shape(aia_pos_rem))
    rib_pos_2 = np.zeros(np.shape(aia_pos_rem))

    for i in range(len(aia_pos_rem)):
        for j in range(len(aia_pos_rem[0])):
            for k in range(len(aia_pos_rem[1])):
                if aia_pos_rem[i, j, k] == 1 and k < split_pos:
                    rib_pos_1[i, j, k] = 1
                elif aia_pos_rem[i, j, k] == 1 and k > split_pos:
                    rib_pos_2[i, j, k] = 1

    rib_neg_1 = np.zeros(np.shape(aia_neg_rem))
    rib_neg_2 = np.zeros(np.shape(aia_neg_rem))

    for i in range(len(aia_neg_rem)):
        for j in range(len(aia_neg_rem[0])):
            for k in range(len(aia_neg_rem[1])):
                if aia_neg_rem[i, j, k] == 1 and k < split_neg:
                    rib_neg_1[i, j, k] = 1
                elif aia_neg_rem[i, j, k] == 1 and k > split_neg:
                    rib_neg_2[i, j, k] = 1

    return rib_pos_1, rib_pos_2, rib_neg_1, rib_neg_2


def find_rib_coordinates(aia_pos_rem, aia_neg_rem):
    """
    Find coordinates of extreme limits of positive and negative ribbons.

    Parameters
    ----------
    aia_pos_rem : list
        Isolated positive ribbon masks.
    aia_neg_rem : list
        Isolated negative ribbon masks.

    Returns
    -------
    lr_coord_neg : list
        Extreme limits of negative ribbon for each time step.
    lr_coord_pos : list
        Extreme limits of positive ribbon for each time step.

    """
    lr_coord_pos = np.zeros([len(aia_pos_rem), 4])
    lr_coord_neg = np.zeros([len(aia_neg_rem), 4])

    for i in range(len(aia_pos_rem)):
        left_x = 0
        left_y = 0
        right_x = 0
        right_y = 0
        for k in range(len(aia_pos_rem[1])):
            for j in range(len(aia_pos_rem[0])):
                # Extreme limit to the left of the ribbon, in positive ribbon
                if aia_pos_rem[i, j, k] == 1:
                    left_x = k
                    left_y = j
                    break
            if left_x != 0:
                break
        for k in range(len(aia_pos_rem[1]) - 1, 0, -1):
            for j in range(len(aia_pos_rem[0])):
                # Extreme limit to the right of the ribbon, in positive ribbon
                if aia_pos_rem[i, j, k] == 1:
                    right_x = k
                    right_y = j
                    break
            if right_x != 0:
                break
        lr_coord_pos[i, :] = [left_x, left_y, right_x, right_y]

    for i in range(len(aia_neg_rem)):
        left_x = 0
        left_y = 0
        right_x = 0
        right_y = 0
        for k in range(len(aia_neg_rem[1])):
            for j in range(len(aia_neg_rem[0])):
                # Extreme limit to the left of the ribbon, in negative ribbon
                if aia_neg_rem[i, j, k] == 1:
                    left_x = k
                    left_y = j
                    break
            if left_x != 0:
                break
        for k in range(len(aia_neg_rem[1]) - 1, 0, -1):
            for j in range(len(aia_neg_rem[0])):
                # Extreme limit to the right of the ribbon, in negative ribbon
                if aia_neg_rem[i, j, k] == 1:
                    right_x = k
                    right_y = j
                    break
            if right_x != 0:
                break
        lr_coord_neg[i, :] = [left_x, left_y, right_x, right_y]

    return lr_coord_neg, lr_coord_pos


def sort_pil(ivs_lim, dvs_lim):
    """
    Sort PIL coordinates in ascending order.

    Parameters
    ----------
    ivs_lim : list
        Vetted x-values for PIL polynomial.
    dvs_lim : list
        Vetted y-values for PIL polynomial.

    Returns
    -------
    ivs_sort : list
        Sorted x-values for PIL polynomial.
    dvs_sort : list
        Sorted y-values for PIL polynomial.
    sortedpil : list
        Sorted ordered pairs for PIL polynomial.

    """
    pil_sort = np.vstack((ivs_lim, dvs_lim)).T
    sortedpil = pil_sort[pil_sort[:, 0].argsort()]
    ivs_sort = sortedpil[:, 0]
    dvs_sort = sortedpil[:, 1]

    return ivs_sort, dvs_sort, sortedpil


def elon_dist_arrays(lr_coord_pos, lr_coord_neg, ivs_lim, dvs_lim, ivs_sort,
                     dvs_sort):
    """
    Create array for distances of limits of ribbon masks from PIL for each
    time step.

    Parameters
    ----------
    lr_coord_neg : list
        Extreme limits of negative ribbon for each time step.
    lr_coord_pos : list
        Extreme limits of positive ribbon for each time step.
    ivs_lim : list
        Vetted x-values for PIL polynomial.
    dvs_lim : list
        Vetted y-values for PIL polynomial.
    ivs_sort : list
        Sorted x-values for PIL polynomial.
    dvs_sort : list
        Sorted y-values for PIL polynomial.

    Returns
    -------
    pil_right_near_pos : list
        Closest PIL point to the "right" edge of positive ribbon for each time
        step.
    pil_left_near_pos : list
        Closest PIL point to the "left" edge of positive ribbon for each time
        step.
    pil_right_near_neg : list
        Closest PIL point to the "right" edge of negative ribbon for each time
        step.
    pil_left_near_neg : list
        Closest PIL point to the "left" edge of negative ribbon for each time
        step.

    """
    left_pil_dist_pos = np.zeros([len(lr_coord_pos), len(ivs_sort)])
    right_pil_dist_pos = np.zeros([len(lr_coord_pos), len(ivs_sort)])
    left_pil_dist_neg = np.zeros([len(lr_coord_neg), len(ivs_sort)])
    right_pil_dist_neg = np.zeros([len(lr_coord_neg), len(ivs_sort)])
    pil_left_near_neg = np.zeros([len(left_pil_dist_neg), 3])
    pil_right_near_neg = np.zeros([len(right_pil_dist_neg), 3])
    pil_left_near_pos = np.zeros([len(left_pil_dist_pos), 3])
    pil_right_near_pos = np.zeros([len(right_pil_dist_pos), 3])

    # The following loops determine the distance from all limit pixels to all
    # pixels corresponding to the PIL and stores in arrays.  The first three
    # loops correspond to the positive ribbon, for all time steps.
    for i in range(len(lr_coord_pos)):
        left_x, left_y, right_x, right_y = lr_coord_pos[i]
        for j in range(len(ivs_sort)):
            left_pil_dist_pos[i, j] = \
                np.sqrt((left_x - ivs_sort[j])**2 + (left_y - dvs_sort[j])**2)
            right_pil_dist_pos[i, j] = \
                np.sqrt((right_x - ivs_sort[j])**2 + (right_y -
                                                      dvs_sort[j])**2)

    # The minimum of the distances to each of the extreme limits is found.
    for i in range(len(left_pil_dist_pos)):
        ind = np.where(left_pil_dist_pos[i] == np.min(left_pil_dist_pos[i]))
        pil_left_near_pos[i, :] \
            = [ivs_lim[ind[0][0]], dvs_sort[ind[0][0]], ind[0][0]]

    for j in range(len(right_pil_dist_pos)):
        ind = np.where(right_pil_dist_pos[j] == np.min(right_pil_dist_pos[j]))
        pil_right_near_pos[j, :] \
            = [ivs_lim[ind[0][0]], dvs_sort[ind[0][0]], ind[0][0]]

    # Identical method as above, for the negative ribbon.
    for i in range(len(lr_coord_neg)):
        left_x, left_y, right_x, right_y = lr_coord_neg[i]
        for j in range(len(ivs_sort)):
            left_pil_dist_neg[i, j] \
                = np.sqrt((left_x - ivs_sort[j])**2 + (left_y -
                                                       dvs_sort[j])**2)
            right_pil_dist_neg[i, j] \
                = np.sqrt((right_x - ivs_sort[j])**2 + (right_y -
                                                        dvs_sort[j])**2)

    for i in range(len(left_pil_dist_neg)):
        ind = np.where(left_pil_dist_neg[i] == np.min(left_pil_dist_neg[i]))
        pil_left_near_neg[i, :] \
            = [ivs_lim[ind[0][0]], dvs_sort[ind[0][0]], ind[0][0]]

    for j in range(len(right_pil_dist_neg)):
        ind = np.where(right_pil_dist_neg[j] == np.min(right_pil_dist_neg[j]))
        pil_right_near_neg[j, :] \
            = [ivs_lim[ind[0][0]], dvs_sort[ind[0][0]], ind[0][0]]

    return pil_right_near_pos, pil_left_near_pos, pil_right_near_neg, \
        pil_left_near_neg


def elongation(pil_right_near_pos, pil_left_near_pos, pil_right_near_neg,
               pil_left_near_neg, sortedpil):
    """
    Script determining the perpendicular extent of positive and negative
    ribbons for each time step.

    Parameters
    ----------
    pil_right_near_pos : list
        Closest PIL point to the "right" edge of positive ribbon for each time
        step.
    pil_left_near_pos : list
        Closest PIL point to the "left" edge of positive ribbon for each time
        step.
    pil_right_near_neg : list
        Closest PIL point to the "right" edge of negative ribbon for each time
        step.
    pil_left_near_neg : list
        Closest PIL point to the "left" edge of negative ribbon for each time
        step.
    sortedpil : list
        Sorted ordered pairs for PIL polynomial.

    Returns
    -------
    lens_pos : list
        Parallel extent of positive ribbon for each time step.
    lens_neg : list
        Parallel extent of negative ribbon for each time step.

    """
    lens_pos = []
    lens_neg = []

    # The curve length of the PIL between two points - the closest to each of
    # the limits of each of the ribbon - is used as the elongation value for
    # each time step.
    for i in range(len(pil_right_near_pos)):
        leftin = int(pil_left_near_pos[i, 2])
        rightin = int(pil_right_near_pos[i, 2])
        curvei = sortedpil[leftin:rightin, :]
        lens_pos.append(curve_length(curvei))

    for i in range(len(pil_right_near_neg)):
        leftin = int(pil_left_near_neg[i, 2])
        rightin = int(pil_right_near_neg[i, 2])
        curvei = sortedpil[leftin:rightin, :]
        lens_neg.append(curve_length(curvei))

    return lens_pos, lens_neg


def convert_to_Mm(lens_pos, dist_pos, lens_neg, dist_neg, conv_f):
    """
    Converts elongation and separation values, determined beforehand, to units
    of Mm.  Also determines time derivative (rate of change) of perpendicular
    and parallel motion values.

    Parameters
    ----------
    lens_pos : list
        Parallel extent of positive ribbon for each time step.
    dist_pos : list
        Perpendicular motion values for positive ribbon, either median or mean.
    lens_neg : list
        Parallel extent of negative ribbon for each time step.
    dist_neg : list
        Perpendicular motion values for negative ribbon, either median or mean.
    conv_f : float
        Conversion factor from pixels to Mm.

    Returns
    -------
    lens_pos_Mm : list
        Perpendicular extent of positive ribbon for each time step, in Mm.
    lens_neg_Mm : list
        Parallel extent of positive ribbon for each time step in Mm.
    distpos_Mm : list
        Perpendicular extent of negative ribbon for each time step, in Mm.
    distneg_Mm : list
        Parallel extent of positive ribbon for each time step in Mm.
    dneg_len : list
        Time derivative of negative ribbon elongation.
    dpos_len : list
        Time derivative of positive ribbon elongation.
    dneg_dist : list
        Time derivative of negative ribbon separation.
    dpos_dist : list
        Time derivative of positive ribbon separation.

    """
    lens_pos_Mm = np.zeros(np.shape(lens_pos))
    lens_neg_Mm = np.zeros(np.shape(lens_neg))
    distpos_Mm = np.zeros(np.shape(dist_pos))
    distneg_Mm = np.zeros(np.shape(dist_neg))

    for i in range(len(lens_pos)):
        lens_pos_Mm[i] = lens_pos[i] * conv_f
        lens_neg_Mm[i] = lens_neg[i] * conv_f
        distpos_Mm[i] = dist_pos[i] * conv_f
        distneg_Mm[i] = dist_neg[i] * conv_f

    dneg_len = np.diff(lens_neg_Mm) / 24.
    dpos_len = np.diff(lens_pos_Mm) / 24.
    dneg_dist = np.diff(distneg_Mm) / 24.
    dpos_dist = np.diff(distpos_Mm) / 24.

    return lens_pos_Mm, lens_neg_Mm, distpos_Mm, distneg_Mm, dneg_len, \
        dpos_len, dneg_dist, dpos_dist


def prep_304_1600_parameters(sav_data_aia, sav_data, eventindices, flnum,
                             start304, peak304, end304, times304, curves304,
                             outflag=0):
    """
    Preps parameters for 304 Angstrom images, in addition to some datetime
    processing for 1600 Angstrom SDO/AIA images.

    Parameters
    ----------
    sav_data_aia : list
        AIA 1600 images processed from .sav file.
    sav_data : list
        HMI images processed from .sav file.
    eventindices : list
        RibbonDB event indices for pre-determined best-performing flares
        relative to approximate rise and decay phase models.
    flnum : int
        Event index for flare in question.
    start304 : list
        Array containing the start times for the flares in best304.
    peak304 : list
        Array containing the peak times for the flares in best304.
    end304 : list
        Array containing the end times for the flares in best304.
    times304 : list
        Time points for all flares in best304.
    curves304 : list
        Light curves for all flares in best304.
    outflag: int, optional
        Flag if the flare is not in the original list of "best performing".
        The number corresponds to RibbonDB flare number.  The default is 0, in
        which case the flare exists in the database.

    Returns
    -------
    startin : int
        Array index for the start of the flare.
    peakin : int
        Array index for the peak of the flare.
    endin : int
        Array index for the end of the flare.
    times : arr
        Array of times for the flare, from AIA datafile.
    s304 : int
        Nearest index in AIA data to start time of the flare from EUV 304 light
        curve.
    e304 : int
        Nearest index in AIA data to end time of the flare from EUV 304 light
        curve.
    filter_304 : list
        Smoothed 304 light curve using scipy's medfilt function, with kernel
        size of 5.
    med304 : float
        Median of 304 Angstrom light curve.
    std304 : float
        Standard deviation of 304 Angstrom light curve.
    timelab : list
        Preparation of time labels for future plotting of light curves.
    aiadat : list
        AIA data for each time step.
    nt : int
        Number of time steps (or images).
    dn1600 : list
        Datenum values for 1600 Angstrom data.
    time304 : list
        Times corresponding to selected flare from 304 Angstrom data.
    times1600 : list
        Times corresponding to selected flare from 1600 Angstrom data.

    """
    xlo = sav_data_aia.x1los
    xhi = sav_data_aia.x2los
    ylo = sav_data_aia.y1los
    yhi = sav_data_aia.y2los

    aiadat = sav_data_aia.aia1600
    time = sav_data.tim

    nt = len(time)
    nx = aiadat.shape[1]
    ny = aiadat.shape[2]
    t1 = str(sav_data.tim[0])
    t2 = str(sav_data.tim[-1])

    # Conversion of string times into usable floats
    tst = float(t1[14:15:1]) + (float(t1[17:18:1])/60) +\
        (float(t1[20:24:1])/3600)
    tend = float(t2[14:15:1]) + (float(t2[17:18:1])/60) +\
        (float(t2[20:24:1])/3600)
    times = np.linspace(tst, tend, nt)

    x = np.linspace(xlo, xhi, nx)
    y = np.linspace(ylo, yhi, ny)
    x, y = np.meshgrid(x, y)

    times1600 = np.empty(nt, dtype=datetime.datetime)
    sum1600 = np.empty(nt)
    dn1600 = np.empty(nt)

    for i in range(nt):
        timechoi = str(sav_data.tim[i])
        times1600[i] = datetime.datetime.strptime(timechoi[2:21],
                                                  '20%y-%m-%dT%H:%M:%S')
        dn1600[i] = datenum(times1600[i])
        timestep = aiadat[i, :, :]
        sum1600[i] = timestep.sum()

    # if flare not in list
    if outflag == 1242:
        file1242 = '/Users/coletamburri/Desktop/Impulsiveness_Paper/imp_dev/twelvefortytwo.mat'
        ev304 = sio.loadmat(file1242)

        curve304_0 = ev304['smspl']
        time304_0 = ev304['windowthr']
        st304 = ev304['tst']
        peak304 = ev304['maxt']
        end304 = ev304['tend']

        curve304 = []
        time304 = []
        for i in range(len(curve304_0)):
            curve304.append(curve304_0[i][0])
            time304.append(time304_0[0][i])

        startin = np.where(dn1600 == find_nearest(dn1600, st304))
        peakin = np.where(dn1600 == find_nearest(dn1600, peak304))
        endin = np.where(dn1600 == find_nearest(dn1600, end304))

    if outflag == 1953:
        file1953 = '/Users/coletamburri/Desktop/Impulsiveness_Paper/imp_dev/nineteenfiftythree.mat'
        ev304 = sio.loadmat(file1953)

        curve304_0 = ev304['smspl']
        time304_0 = ev304['windowthr']
        st304 = ev304['tst']
        peak304 = ev304['maxt']
        end304 = ev304['tend']

        curve304 = []
        time304 = []
        for i in range(len(curve304_0)):
            curve304.append(curve304_0[i][0])
            time304.append(time304_0[0][i])

        startin = np.where(dn1600 == find_nearest(dn1600, st304))
        peakin = np.where(dn1600 == find_nearest(dn1600, peak304))
        endin = np.where(dn1600 == find_nearest(dn1600, end304))
        
    elif outflag == 0:
        # Find index of nearest index to flare number in 304 flares array
        ind = (np.isclose(eventindices, flnum))
        index = np.where(ind)[0][0]

        # Light curve for selected flare
        curve304 = curves304[index]
        time304 = times304[index]
        # Time indices for 1600A data - time series not identical
        startin = np.where(dn1600 == find_nearest(dn1600, start304[ind][0]))
        peakin = np.where(dn1600 == find_nearest(dn1600, peak304[ind][0]))
        endin = np.where(dn1600 == find_nearest(dn1600, end304[ind][0]))

    # Integrate over all pixels in 1600A line
    for i in range(nt):
        timestep = aiadat[i, :, :]
        sum1600[i] = timestep.sum()

    for i in range(nt):
        timechoi = str(sav_data.tim[i])
        times1600[i] = datetime.datetime.strptime(timechoi[2:21],
                                                  '20%y-%m-%dT%H:%M:%S')

    # Time indices for 304 - nearest to dn1600 points found
    s304 = find_nearest_ind(time304, min(dn1600))
    e304 = find_nearest_ind(time304, max(dn1600))
    filter_304 = scipy.signal.medfilt(curve304, kernel_size=5)

    med304 = np.median(curve304)
    std304 = np.std(curve304)

    # Remove 304 Angstrom pixels below a threshold - these will be outliers.
    # Only applies to one flare studied as of 14 March 2022
    for i in range(len(curve304)):
        if curve304[i] < 0.54:
            curve304[i] = 'NaN'

    timelab = np.empty(nt)

    timelabs = range(0, 24 * len(times), 24)

    for i in range(len(timelabs)):
        timelab[i] = timelabs[i] / 60

    return startin, peakin, endin, times, s304, e304, filter_304, med304, \
        std304, timelab, aiadat, nt, dn1600, time304, times1600


def img_mask(aia8_pos, aia8_neg, aiadat, nt):
    """
    Mapping of positive and negative masks onto AIA data and sum of pixel
    numbers for each ribbon for each time step.

    Parameters
    ----------
    aia8_pos : list
        Contains only the positive ribbon masks for each time step.
    aia8_neg : list
        Contains only the negative ribbon masks for each time step.
    aiadat : list
        AIA 1600 Angstrom data for each time step.
    nt : int
        Number of time steps for the flare.

    Returns
    -------
    posrib : list
        Time step images for positive ribbon.
    negrib : list
        Time step images for negative ribbon.
    pos1600 : list
        Summed pixel numbers in positive ribbon for each time step.
    neg1600 : list
        Summed pixel numbers in negative ribbon for each time step.

    """

    posrib = np.zeros(np.shape(aia8_pos))
    negrib = np.zeros(np.shape(aia8_neg))

    # Positive ribbon masks - actual values from AIA, not 0/1
    for i in range(len(aia8_pos)):
        posrib[i, :, :] = aia8_pos[i, :, :] * aiadat[i, :, :]

    # Negative
    for j in range(len(aia8_neg)):
        negrib[j, :, :] = aia8_neg[j, :, :] * aiadat[j, :, :]

    pos1600 = np.empty(nt)
    neg1600 = np.empty(nt)

    for i in range(nt):
        timesteppos = posrib[i, :, :]
        pos1600[i] = timesteppos.sum()
        timestepneg = negrib[i, :, :]
        neg1600[i] = timestepneg.sum()

    return posrib, negrib, pos1600, neg1600


def load_from_file(flnum, pick=True):
    """
    Option to load separation and elongation values from saved values in file.

    Parameters
    ----------
    flnum : int
        RibbonDB index for selected flare.
    pick : bool, optional
        allow_pickle input for np.load. The default is True.

    Returns
    -------
    dt1600 : list
        1600 Angstrom datetime values for flare.
    pos1600 : list
        Summed pixels in positive ribbon for 1600 Angstrom data.
    neg1600 : list
        Summed pixels in positive ribbon for 1600 Angstrom data.
    time304 : list
        Time series for 304 Angstrom data.
    filter_304 : list
        Smoothed 304 Angstrom light curve.
    distpos_Mm : list
        Separation values for positive ribbon in Mm.
    distneg_Mm : list
        Separation values for negative ribbon in Mm.
    lens_pos_Mm : list
        Elongation values for positive ribbon in Mm.
    lens_neg_Mm : list
        Elongation values for negative ribbon in Mm.
    ivs : list
        x-coordinates for PIL polynomial.
    dvs : list
        y-coordinates for PIL polynomial.

    """

    # Pickle is not ideal, but all data in these files are only variables saved
    # by Cole Tamburri, Spring 2022
    ev = np.load(flnum, allow_pickle=pick)

    dt1600 = ev['dt1600']
    pos1600 = ev['pos1600']
    neg1600 = ev['neg1600']
    time304 = ev['time304']
    filter_304 = ev['filter_304']
    distpos_Mm = ev['distpos_Mm']
    distneg_Mm = ev['distneg_Mm']
    lens_pos_Mm = ev['lens_pos_Mm']
    lens_neg_Mm = ev['lens_neg_Mm']
    ivs = ev['ivs']
    dvs = ev['dvs']

    return dt1600, pos1600, neg1600, time304, filter_304, distpos_Mm, \
        distneg_Mm, lens_pos_Mm, lens_neg_Mm, ivs, dvs


def elon_periods(dpos_len, dneg_len, pos_crit=1, neg_crit=1, zer_pos_c=2,
                 zer_neg_c=2, n_min=1, m_min=1):
    """
    Determination of periods of elongation for positive and negative ribbons
    from time series.

    Parameters
    ----------
    dpos_len : list
        Time derivative of positive ribbon elongation.
    dneg_len : list
        Time derivative of negative ribbon elongation.
    pos_crit : int, optional
        Number of points beyond which an "extended period" of elongation is
        defined, positive ribbon. The default is 1.
    neg_crit : int, optional
        Number of points beyond which an "extended period" of elongation is
        defined, negative ribbon. The default is 1.
    zer_pos_c : int, optional
        Number of zero-derivative points beyond which an "extended period" of
        elongation is said to end, positive ribbon. The default is 2.
    neg_crit : int, optional
        Number of points beyond which an "extended period" of elongation is
        defined, negative ribbon. The default is 1.
    Returns
    -------
    elonperiod_start_pos : list
        Determined start times for elongation in positive ribbon.
    elonperiod_end_pos : list
        Determined end times for elongation in positive ribbon.
    elonperiod_start_neg : list
        Determined start times for elongation in negative ribbon.
    elonperiod_end_neg : list
        Determined end times for elongation in negative ribbon.

    """
    elonfiltpos = dpos_len
    elonfiltneg = dneg_len
    elonperiod_start_pos = []
    elonperiod_end_pos = []
    elonperiod_start_neg = []
    elonperiod_end_neg = []
    n = 0
    m = 0
    zer_n = 0
    zer_m = 0

    for i in range(len(elonfiltpos)):
        if elonfiltpos[i] > 0:
            n += 1
            if n == 1:
                time = i
            # Tripped if extended period of elongation, not already recorded
            if n > pos_crit and time not in elonperiod_start_pos:
                elonperiod_start_pos.append(time)
        elif elonfiltpos[i] <= 0:
            if n > n_min:
                zer_n += 1
                # If rate of change returns to 0 for several points
                if zer_n > zer_pos_c:
                    elonperiod_end_pos.append(i)
                    n = 0
                    zer_n = 0
            else:
                n = 0
                continue

    # Comments identical to above method
    for j in range(len(elonfiltneg)):
        if elonfiltneg[j] > 0:
            m += 1
            if m == 1:
                time = j
            if m > neg_crit and time not in elonperiod_start_neg:
                elonperiod_start_neg.append(time)
        elif elonfiltneg[j] <= 0:
            if m > m_min:
                zer_m += 1
                if zer_m > zer_neg_c:
                    elonperiod_end_neg.append(j)
                    m = 0
                    zer_m = 0
            elif zer_m > 1:
                m = 0
                continue

    # Remove repeated values
    elonperiod_start_pos = list(set(elonperiod_start_pos))
    elonperiod_end_pos = list(set(elonperiod_end_pos))
    elonperiod_start_neg = list(set(elonperiod_start_neg))
    elonperiod_end_neg = list(set(elonperiod_end_neg))

    elonperiod_start_pos.sort()
    elonperiod_end_pos.sort()
    elonperiod_start_neg.sort()
    elonperiod_end_neg.sort()

    return elonperiod_start_pos, elonperiod_end_pos, elonperiod_start_neg, \
        elonperiod_end_neg


def sep_periods(dpos_dist, dneg_dist, start=20, kernel_size=3, pos_crit=3,
                neg_crit=3, zer_pos_c=3, zer_neg_c=3):
    """
    Determination of periods of separation for each ribbon from time
    derivatives of separation data.

    Parameters
    ----------
    dneg_dist : list
        Time derivative of negative ribbon separation.
    dpos_dist : list
        Time derivative of positive ribbon separation.
    start : int
        Index where to start plotting and processing.
    kernel_size : int, optional
        Kernel size for scipy medfilt of separation curves. The default is 3.

    Returns
    -------
    sepperiod_start_pos : list
        Start times for periods of extended separation, positive ribbon.
    sepperiod_end_pos : list
        End times for periods of extended separation, positive ribbon.
    sepperiod_start_neg : list
        Start times for periods of extended separation, negative ribbon.
    sepperiod_end_neg : list
        End times for periods of extended separation, negative ribbon.

    """

    # Separation values are much more variable, so smoothing is necessary
    sepfiltpos = scipy.signal.medfilt(dpos_dist, kernel_size=kernel_size)
    sepfiltneg = scipy.signal.medfilt(dneg_dist, kernel_size=kernel_size)

    sepperiod_start_pos = []
    sepperiod_end_pos = []
    sepperiod_start_neg = []
    sepperiod_end_neg = []
    n = 0
    m = 0
    for i in range(start, len(sepfiltpos)):
        if sepfiltpos[i] > 0:
            n += 1
            if n == 1:
                time = i
            # Append if extended period of separation - more stringent than
            # elongation, justified by above smoothing
            if n > pos_crit and time not in sepperiod_start_pos:
                sepperiod_start_pos.append(time)
        elif sepfiltpos[i] <= 0:
            # Identify if rate of change has not changed for some time
            if n > zer_pos_c:
                sepperiod_end_pos.append(i)
                n = 0
            else:
                n = 0
                continue

    # Same method as above
    for i in range(start, len(sepfiltneg)):
        if sepfiltneg[i] > 0:
            m += 1
            if m == 1:
                time = i
            if m > neg_crit and time not in sepperiod_start_neg:
                sepperiod_start_neg.append(time)
        elif sepfiltneg[i] <= 0:
            if m > zer_neg_c:
                sepperiod_end_neg.append(i)
                m = 0
            else:
                m = 0
                continue

    # Remove repeated values
    sepperiod_start_pos = list(set(sepperiod_start_pos))
    sepperiod_end_pos = list(set(sepperiod_end_pos))
    sepperiod_start_neg = list(set(sepperiod_start_neg))
    sepperiod_end_neg = list(set(sepperiod_end_neg))

    sepperiod_start_pos.sort()
    sepperiod_end_pos.sort()
    sepperiod_start_neg.sort()
    sepperiod_end_neg.sort()

    return sepperiod_start_pos, sepperiod_end_pos, sepperiod_start_neg, \
        sepperiod_end_neg


def prep_times(dn1600, time304):
    """
    Convert datenum to datetime values for animation and presentation.

    Parameters
    ----------
    dn1600 : list
        Datenum values for 1600 Angstrom data.
    time304 : list
        Datenum values for 304 Angstrom data.

    Returns
    -------
    dt1600 : list
        Datetime values for 1600 Angstrom data.
    dt304 : list
        Datetime values for 304 Angstrom data.

    """
    dt1600 = []
    dt304 = []
    for i in range(len(dn1600)):
        dt1600.append(datenum_to_datetime(dn1600[i]))

    for i in range(len(time304)):
        if np.isnan(time304[i]):
            dt304.append(datenum_to_datetime(time304[0]))
        else:
            dt304.append(datenum_to_datetime(time304[i]))

    return dt1600, dt304

# BEGIN PLOTTING ROUTINES #


def lc_plot(times, nt, time304, filter_304, s304, e304, dn1600, pos1600,
            neg1600, lens_pos_Mm, lens_neg_Mm, distpos_Mm, distneg_Mm, aiadat,
            hmi_cumul_mask1, dt304, timelab, conv_f, ivs, dvs, year, mo, day,
            arnum, xcl, xclnum, X, Y, xarr_Mm, yarr_Mm, dt1600, flag=1,
            stsep=25, stelon=1, lolim=0, hilim=1):
    """
    Animation plotting, with 1600 Angstrom, 304 Angstrom, and HMI data.

    Parameters
    ----------
    times : list
        Times corresponding to each AIA time step.
    nt : list
        Number of image times.
    time304 : list
        Series of times for 304 Angstrom data.
    filter_304 : list
        Smoothed 304 Angstrom data.
    s304 : int
        Nearest index in AIA data to start time of the flare from EUV 304 light
        curve.
    e304 : int
        Nearest index in AIA data to end time of the flare from EUV 304 light
        curve.
    dn1600 : list
        Datenum values for 1600 Angstrom data.
    pos1600 : list
        Summed pixel numbers in positive ribbon for each time step.
    neg1600 : list
        Summed pixel numbers in negative ribbon for each time step.
    lens_pos_Mm : list
        Perpendicular extent of positive ribbon for each time step, in Mm.
    lens_neg_Mm : list
        Parallel extent of positive ribbon for each time step in Mm.
    distpos_Mm : list
        Perpendicular extent of negative ribbon for each time step, in Mm.
    distneg_Mm : list
        Parallel extent of positive ribbon for each time step in Mm.
    aiadat : list
        AIA data for each time step.
    hmi_cumul_mask1 : list
        Cumulative magnetic field strength masking estimates for all flare
        images.
    dt304 : list
        Datetime values for 304 Angstrom data.
    timelab : list
        Points for labels in time axis.
    conv_f : float
        Conversion factor from pixels to Mm.
    ivs : list
        x-coordinates for PIL polynomial
    dvs : list
        y-coordinates for PIL polynomial
    year : int
        Year of event.
    mo : int
        Month of event.
    day : int
        Day of event.
    arnum : int
        Active region number.
    xcl : str
        x-ray class identifier (C, M, X)
    xclnum : float
        x-ray class identifier, number.
    X : list
        Meshgrid of x values for image coordinates.
    Y : list
        Meshgrid of y values for image coordinates.
    xarr_Mm : list
        x-coordinates, in megameters.
    yarr_Mm : list
        y-coordinates, in megameters.
    dt1600 : list
        Datetime values for 1600 Angstrom data.
    flag : int, optional
        0 (plot only first five frames) or 1 (plot all frames). The default is
        1.
    stsep: int, optional
        Start index for separation curve.  Default is 25.
    stelon: int, optional
        Start index for elongation curve.  Default is 1.
    lolim: float, optional
        Start value for second y-axis on sep/elon plot.  Default is 0.
    hilim: float, optional
        End value for second y-axis on sep/elon plot.  Default is 1, which
        triggers 140*conv_f.

    Returns
    -------
    col1 : list
        List comprising AIA data plot.
    col2 : list
        List comprising AIA contourmap.
    lc304 : list
        List comprising 304 Angstrom light curve.
    lc1600 : list
        List comprising 1600 Angstrom light curve.
    sep : list
        List comprising positive ribbon separation plot.
    sep2 : list
        List comprising negative ribbon separation plot.
    elon : list
        List comprising positive ribbon elongation plot.
    elon2 : list
        List comprising negative ribbon elongation plot.

    """

    if hilim == 1:
        hilim = 140*conv_f

    # Extremes of chromospheric line light curves
    min304 = min(filter_304[s304: e304])
    max304 = max(filter_304[s304: e304])
    minpos1600 = min(pos1600)
    maxpos1600 = max(pos1600)
    minneg1600 = min(neg1600)
    maxneg1600 = max(neg1600)

    # Normalize for light curve comparison
    norm304 = (filter_304 - min304) / (max304 - min304)
    normpos1600 = (pos1600 - minpos1600) / (maxpos1600 - minpos1600)
    normneg1600 = (neg1600 - minneg1600) / (maxneg1600 - minneg1600)
    scalefac = max(pos1600) / max(neg1600)

    # Initialize plot
    fig = plt.figure(figsize=(25, 12))
    gs = fig.add_gridspec(9, 9)
    ax1 = fig.add_subplot(gs[:, 5:])
    ax2 = fig.add_subplot(gs[0:4, 0:4])
    ax0 = fig.add_subplot(gs[5:, 0:4])

    # Elongation plots
    lns1 = ax0.plot(dn1600[stelon:], lens_pos_Mm[stelon:], '-+', c='red',
                    markersize=10, label='Pos. Elongation')
    lns2 = ax0.plot(dn1600[stelon:], lens_neg_Mm[stelon:], '-+', c='blue',
                    markersize=10, label='Neg. Elongation')
    ax5 = ax0.twinx()
    ax5.cla()
    lns4 = ax5

    # Separation plots
    lns5 = ax0.plot(dt1600[stsep:], distpos_Mm[stsep:], '-.', c='red',
                    markersize=10, label='Pos. Separation')
    ax0.plot(dt1600[stsep:], distneg_Mm[stsep:], '-.', c='blue',
             markersize=10, label='Neg. Separation')

    # Plot 1600 Angstrom pcolormesh images, as well as HMI images
    col1 = ax1.pcolormesh(X, Y, np.log10(aiadat[0, :, :]), cmap='pink',
                          shading='auto')
    col2 = ax1.contour(X, Y, hmi_cumul_mask1[0, :, :], cmap='seismic')

    # Plot 304 Angstrom light curve
    lc304 = ax2.plot(dt304, norm304, color = 'black', linewidth=1,
                     label=r'Norm. 304\AA Light Curve')
    ax3 = ax2.twinx()

    # Plot 1600 Angstrom light curve
    lc1600 = ax3.plot(dt1600, normpos1600, linewidth=3, color='red',
                      label=r'Norm. 1600\AA Light Curve, +')
    lc1600 = ax3.plot(dt1600, normneg1600, linewidth=3, color='blue',
                      label=r'Norm. 1600\AA Light Curve, +')

    ax1.set_title(str(year) + "-" + str(mo) + "-" + str(day) + ", AR" +
                  str(arnum)+"\n" + xcl + str(xclnum) + " Class Flare\n",
                  font='Times New Roman', fontsize=45)
    ax2.set_title(r'304\AA and 1600\AA Light Curves', fontsize=45)

    ax0.set_title('Ribbon Separation and Elongation', fontsize=45)
    ax0.legend(fontsize=15)
    ax0.grid()
    ax2.set_xlim([dn1600[0], dn1600[-1]])
    ax3.set_xlim([dn1600[0], dn1600[-1]])
    ax0.set_xlim([timelab[0], timelab[-1]])

    # Plot PIL on 1600 Angstrom and HMI panel
    ax1.scatter(ivs, dvs, color='k', s=1)
    lines, labels = ax0.get_legend_handles_labels()
    lines2, labels2 = ax5.get_legend_handles_labels()
    ax0.legend(lines + lines2, labels + labels2)

    ax5.set_ylim([lolim, hilim])

    def animate(t):
        ax1.cla()
        ax2.cla()
        ax0.cla()
        ax5 = ax0.twinx()
        ax5.cla()

        # Plot 1600 Angstrom image
        col1 = ax1.pcolormesh(X, Y, np.log10(aiadat[t, :, :]), cmap='pink',
                              shading='auto')

        # HMI contour over 1600 Angstrom image
        col2 = ax1.contour(X, Y, hmi_cumul_mask1[t, :, :], cmap='seismic')
        ax1.set_xlabel(r'Horizontal Distance from Image Center [$Mm$]',
                       fontsize=15)
        ax1.set_ylabel(r'Vertical Distance from Image Center [$Mm$]',
                       fontsize=15)

        # Separation curves
        sep = ax0.plot(dt1600[stsep:], distpos_Mm[stsep:], '-.', c='red',
                       markersize=10, label='Pos. Separation')
        sep2 = ax0.plot(dt1600[stsep:], distneg_Mm[stsep:], '-.', c='blue',
                        markersize=10, label='Neg. Separation')
        ax1.scatter((ivs-400) * conv_f, (dvs-400) * conv_f, color='k', s=1)

        # Elongation curves
        elon = ax5.plot(dt1600[stelon:], lens_pos_Mm[stelon:], '-+', c='red',
                        markersize=10, label='Pos. Elongation')
        elon2 = ax5.plot(dt1600[stelon:], lens_neg_Mm[stelon:], '-+',
                         c='blue', markersize=10, label='Neg. Elongation')
        ax1.set_xlim([-250 * conv_f, 250 * conv_f])
        ax1.set_ylim([-250 * conv_f, 250 * conv_f])

        # 304 Angstrom light curve
        lc304 = ax2.plot(dt304, norm304, '-x', color='black', linewidth=1,
                         label=r'304\AA')
        ax3 = ax2.twinx()

        # 1600 Angstrom light curve
        lc1600 = ax3.plot(dt1600, scipy.signal.medfilt(normpos1600,kernel_size=5), linewidth=3, color='red',
                          label=r'1600\AA, +')
        lc1600 = ax3.plot(dt1600, scipy.signal.medfilt(normneg1600, kernel_size=5), linewidth=3, color='blue',
                          label=r'1600\AA, -')

        ax2.set_xlim([dt1600[0], dt1600[-1]])
        ax2.set_ylim([-0.05, 1.05])
        ax3.set_ylim([-0.05, 1.05])

        myFmt = mdates.DateFormatter('%H:%M')
        ax2.xaxis.set_major_formatter(myFmt)
        ax3.xaxis.set_major_formatter(myFmt)
        ax0.xaxis.set_major_formatter(myFmt)
        ax5.xaxis.set_major_formatter(myFmt)
        textstr = r'1600\AA +/- Factor: ' + str(round(scalefac, 3))
        ax2.text(2 * (max(dt1600) - min(dt1600)) / 5 + min(dt1600), 0.1,
                 textstr, fontsize=12, bbox=dict(boxstyle="square",
                                                 facecolor="white",
                                                 ec="k", lw=1,
                                                 pad=0.3))
        ax2.set_xlabel(['Time since 00:00 UT [min], ' + year + '-' + mo + '-'
                        + day], fontsize=15)
        ax2.set_xlabel(['Time since 00:00 UT [min], ' + year + '-' + mo + '-'
                        + day], fontsize=15)
        ax2.set_ylabel(r'Norm. Integ. Count, 1600\AA', color='purple',
                       fontsize=15)

        lines, labels = ax2.get_legend_handles_labels()
        lines2, labels2 = ax3.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc='lower right')
        ax2.grid(linestyle='dashed')
        ax3.grid(linestyle='dashdot')
        ax2.axvline(dt1600[t], linewidth=4, color='black')
        ax0.axvline(dt1600[t], linewidth=4, color='black')
        ax0.axvline(dt1600[t], linewidth=4, color='black')
        ax1.set_title(str(year) + "-" + str(mo) + "-" + str(day) + ", AR" +
                      str(arnum) + ", " + xcl + str(xclnum) + " Class Flare",
                      fontsize=45)
        ax2.set_title(r'304\AA and 1600\AA Light Curves', fontsize=45)
        ax0.set_xlim([dt1600[0], dt1600[-1]])
        ax0.set_xlabel(['Time since 00:00 UT [min], ' + year + '-' + mo + '-' +
                        day], fontsize=15)
        ax0.set_ylabel(r'Separation [$Mm$]', fontsize=15)
        ax5.set_ylabel(r'Elongation [$Mm$]', fontsize=15)
        ax0.set_title('Ribbon Separation and Elongation', fontsize=45)
        ax0.legend(fontsize=15)
        ax0.grid()
        ax1.text(57, 95, str(dt1600[t].hour).zfill(2) + ':' +
                 str(dt1600[t].minute).zfill(2) + '.' +
                 str(dt1600[t].second).zfill(2) + ' UT', fontsize=20,
                 bbox=dict(boxstyle="square", facecolor="white", ec="k",
                           lw=1, pad=0.3))

        lines, labels = ax0.get_legend_handles_labels()
        lines2, labels2 = ax5.get_legend_handles_labels()
        ax0.legend(lines + lines2, labels + labels2, loc='lower right')

        ax5.set_ylim([lolim, hilim])
        return col1, col2, lc304, lc1600, sep, sep2, elon, elon2

    # Option to only include first few frames for debugging purposes
    if flag == 1:
        ani = animat.FuncAnimation(fig, animate, frames=np.shape(aiadat)[0],
                                   interval=20, repeat_delay=0)
    elif flag == 0:
        ani = animat.FuncAnimation(fig, animate, frames=5, interval=20,
                                   repeat_delay=0)

    ani.save(['/Users/coletamburri/Desktop/' + mo + '_' + day + '_' + year + '.gif'],
             dpi=200)

    return None


def mask_plotting(X, Y, pos_rem, neg_rem, xarr_Mm, yarr_Mm, flnum):
    """
    Plotting of HMI image masks.

    Parameters
    ----------
    X : list
        Meshgrid of x values for image coordinates.
    Y : list
        Meshgrid of y values for image coordinates.
    pos_rem : list
        The positive polarity HMI image, with spurs removed.
    neg_rem : list
        The negative polarity HMI image, with spurs removed.
    xarr_Mm : list
        x-coordinates, in megameters.
    yarr_Mm : list
        y-coordinates, in megameters.
    flnum : int
        RibbonDB index of flare in question.

    Returns
    -------
    None.

    """
    fig1, ax1 = plt.subplots(figsize=(6, 6))

    # Plot positive mask, with pixel vetting
    ax1.pcolormesh(X, Y, pos_rem, cmap='bwr', vmin=-1, vmax=1)
    ax1.set_title('Positive Mask', font="Times New Roman", fontsize=22,
                  fontweight='bold')

    ax1.set_xlim([xarr_Mm[200], xarr_Mm[600]])
    ax1.set_ylim([yarr_Mm[200], yarr_Mm[600]])
    ax1.set_xlabel('Horizontal Distance from Image Center [Mm]', fontsize=17)
    ax1.set_ylabel('Vertical Distance from Image Center [Mm]', fontsize=17)

    ax1.tick_params(labelsize=15)

    fig2, ax2 = plt.subplots(figsize=(6, 6))

    # Plot negative mask, with pixel vetting
    ax2.set_xlabel('Horizontal Distance from Image Center [Mm]', fontsize=17)
    ax2.set_ylabel('Vertical Distance from Image Center [Mm]', fontsize=17)
    ax2.tick_params(labelsize=15)
    ax2.pcolormesh(X, Y, neg_rem, cmap='bwr', vmin=-1, vmax=1)

    ax2.set_title('Negative Mask', font="Times New Roman", fontsize=22,
                  fontweight='bold')
    ax2.set_xlim([xarr_Mm[200], xarr_Mm[600]])
    ax2.set_ylim([yarr_Mm[200], yarr_Mm[600]])

    fig1.savefig(str(flnum) + '_pos_mask.png')
    fig2.savefig(str(flnum) + '_neg_mask.png')

    return None


def convolution_mask_plotting(X, Y, hmi_con_pos_c, hmi_con_neg_c, pil_mask_c,
                              xarr_Mm, yarr_Mm, flnum, xlim=[200, 600],
                              ylim=[200, 600]):
    """
    Plots masks, convolved with Gassian of width specified above.

    Parameters
    ----------
    X : list
        Meshgrid of x values for image coordinates.
    Y : list
        Meshgrid of y values for image coordinates.
    hmi_con_pos_c : list
        Positive HMI, convolved with Gaussian.
    hmi_con_neg_c : list
        Negative HMI, convolved with Gaussian.
    pil_mask_c : list
        PIL mask.
    xarr_Mm : list
        x-coordinates, in megameters.
    yarr_Mm : list
        y-coordinates, in megameters.
    flnum : int
        RibbonDB index of flare in question.
    xlim : list, optional
        Limts of x-coordinates to plot. The default is [200,600].
    ylim : list, optional
        Limits of y-coordinates to plot. The default is [200,600].

    Returns
    -------
    None.

    """
    fig1, ax1 = plt.subplots(figsize=(6, 6))
    ax1.pcolormesh(X, Y, hmi_con_pos_c, cmap='bwr', vmin=-1, vmax=1)
    ax1.set_title('Positive Mask Convolution', font="Times New Roman",
                  fontsize=22, fontweight='bold')
    ax1.set_xlim([xarr_Mm[200], xarr_Mm[600]])
    ax1.set_ylim([yarr_Mm[200], yarr_Mm[600]])
    ax1.set_xlabel('Horizontal Distance from Image Center [Mm]', fontsize=17)
    ax1.set_ylabel('Vertical Distance from Image Center [Mm]', fontsize=17)
    ax1.tick_params(labelsize=15)

    fig2, ax2 = plt.subplots(figsize=(6, 6))
    ax2.tick_params(labelsize=15)
    ax2.pcolormesh(X, Y, hmi_con_neg_c, cmap='bwr', vmin=-1, vmax=1)
    ax2.set_xlabel('Horizontal Distance from Image Center [Mm]', fontsize=17)
    ax2.set_ylabel('Vertical Distance from Image Center [Mm]', fontsize=17)
    ax2.set_title('Negative Mask Convolution', font="Times New Roman",
                  fontsize=22, fontweight='bold')
    ax2.set_xlim([xarr_Mm[xlim[0]], xarr_Mm[xlim[1]]])
    ax2.set_ylim([yarr_Mm[ylim[0]], yarr_Mm[ylim[1]]])

    fig3, ax3 = plt.subplots()
    ax3.pcolormesh(X, Y, pil_mask_c)
    ax3.set_title('Polarity Inversion Line Mask', font="Times New Roman",
                  fontsize=22, fontweight='bold')
    ax3.tick_params(labelsize=15)
    ax3.set_xlim([xarr_Mm[xlim[0]], xarr_Mm[xlim[1]]])
    ax3.set_ylim([yarr_Mm[ylim[0]], yarr_Mm[ylim[1]]])

    fig1.savefig(str(flnum) + '_pos_conv_mask.png')
    fig2.savefig(str(flnum) + '_neg_conv_mask.png')
    fig3.savefig(str(flnum) + '_PIL_conv_mask.png')

    return None


def pil_poly_plot(X, Y, pil_mask_c, hmi_dat, ivs, dvs, conv_f, xarr_Mm,
                  yarr_Mm, flnum, xlim=[200, 600], ylim=[200, 600]):
    """
    Plotting of PIL polynomial over mask.

    Parameters
    ----------
    X : list
        Meshgrid of x values for image coordinates.
    Y : list
        Meshgrid of y values for image coordinates.
    pil_mask_c : list
        PIL mask.
    hmi_dat : list
        HMI data image for flare in question.
    ivs : list
        x-coordinates for PIL polynomial.
    dvs : list
        y-coordinates for PIL polynomial.
    conv_f : float
        Conversion factor from pixels to Mm.
    xarr_Mm : list
        x-coordinates, in megameters.
    yarr_Mm : list
        y-coordinates, in megameters.
    flnum : int
        RibbonDB index of flare in question.
    xlim : list, optional
        Limts of x-coordinates to plot. The default is [200,600].
    ylim : list, optional
        Limits of y-coordinates to plot. The default is [200,600].

    Returns
    -------
    None.

    """
    # Generate the plot
    fig, ax = plt.subplots(figsize=(7, 10))

    # show color mesh
    ax.pcolormesh(X, Y, pil_mask_c, cmap='hot')

    # plot the line
    ax.scatter((ivs - 400) * conv_f, (dvs - 400) * conv_f, color='c', s=1)
    hmik = hmi_dat / 1000
    plt.contour(X, Y, hmik, levels=[-3, -1.8, -.6, .6, 1.8, 3],
                cmap='seismic')

    ax.set_xlim([xarr_Mm[xlim[0]], xarr_Mm[xlim[1]]])
    ax.set_ylim([yarr_Mm[ylim[0]], yarr_Mm[ylim[1]]])
    ax.set_xticks([-80, -60, -40, -20, 0, 20, 40, 60, 80])
    ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80])
    cbar = plt.colorbar(orientation='horizontal')
    tick_font_size = 15
    ax.tick_params(labelsize=tick_font_size)
    cbar.ax.tick_params(labelsize=tick_font_size)
    ax.set_xlabel('Horizontal Distance from Image Center [Mm]', fontsize=15)
    ax.set_ylabel('Vertical Distance from Image Center [Mm]', fontsize=15)

    cbar.ax.set_xlabel('HMI Contours [kG]', font='Times New Roman',
                       fontsize=17, fontweight='bold')
    ax.set_title('PIL Mask and Polynomial', font='Times New Roman',
                 fontsize=45, fontweight='bold')
    fig.savefig(str(flnum) + '_pilpolyplot.png')

    return None


def ribbon_sep_plot(dist_pos, dist_neg, times, flnum, pltstrt, dt1600):
    """
    Plot ribbon separation values throughout flare.

    Parameters
    ----------
    dist_pos : list
        Separation values, positive ribbon.
    dist_neg : list
        Separation values, negative ribbon.
    times : list
        Times corresponding to each AIA time step.
    flnum : int
        Flare index from RibbonDB database.
    pltstrt : int
        Index for where to start displaying separation values.

    Returns
    -------
    None.

    """
    timelab = range(0, 24 * len(times), 24)
    fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(13, 15))

    # Plot separation, positive ribbon
    ax1.plot(timelab[pltstrt:], dist_pos[pltstrt:], '-+', c='red',
             markersize=10, label='median')
    ax1.legend(fontsize=15)
    ax1.grid()
    s = str(dt1600[0])
    ax1.set_xlabel('Time [s since ' + s[5:-7] + ']', font='Times New Roman',
                   fontsize=15)
    ax1.set_ylabel('Cartesian Pixel Distance', font='Times New Roman',
                   fontsize=15)
    ax1.set_title('Positive Ribbon Separation', font='Times New Roman',
                  fontsize=45)

    # Plot separation, negative ribbon
    ax2.plot(timelab[pltstrt:], dist_neg[pltstrt:], '-+', c='red',
             markersize=10, label='median')
    ax2.legend(fontsize=15)
    ax2.grid()
    ax2.set_xlabel('Time [s since ' + s[5:-7] + ']', font='Times New Roman',
                   fontsize=15)
    ax2.set_ylabel('Cartesian Pixel Distance', font='Times New Roman',
                   fontsize=15)
    ax2.set_title('Negative Ribbon Separation', font='Times New Roman',
                  fontsize=45)

    fig.savefig(str(flnum) + 'sep_raw_plt.png')

    return None


def ribbon_elon_plot(lens_pos, lens_neg, times, pltstrt, flnum, dt1600):
    """
    Plot ribbon elongation values throughout flare.

    Parameters
    ----------
    dist_pos : list
        Elongation values, positive ribbon.
    dist_neg : list
        Elongation values, negative ribbon.
    times : list
        Times corresponding to each AIA time step.
    flnum : int
        Flare index from RibbonDB database.
    pltstrt : int
        Index for where to start displaying elongation values.
    dt1600 : list
        Array of datetimes for SDO/AIA 1600 Angstrom data

    Returns
    -------
    None.

    """
    timelab = range(0, 24 * len(times), 24)

    fig, ax1 = plt.subplots(figsize=(13, 7))

    # Plot elongation, positive ribbon
    ax1.plot(timelab[pltstrt:], lens_pos[pltstrt:], '-+', c='red',
             markersize=10, label='+ Ribbon')

    # Plot elongation, negative ribbon
    ax1.plot(timelab[pltstrt:], lens_neg[pltstrt:], '-+', c='blue',
             markersize=10, label='- Ribbon')
    ax1.legend(fontsize=15)
    ax1.grid()
    s = str(dt1600[0])
    ax1.set_xlabel('Time [s since ' + s[5:-7] + ']', font='Times New Roman',
                   fontsize=17)
    ax1.set_ylabel('Cartesian Pixel Distance', font='Times New Roman',
                   fontsize=17)
    ax1.set_title('Ribbon Elongation', font='Times New Roman', fontsize=45)

    fig.savefig(str(flnum) + 'elon_raw_plt.png')

    return None


def elon_period_plot(dpos_len, dneg_len, times, times1600, lens_pos_Mm,
                     lens_neg_Mm, flnum, elonperiod_start_neg,
                     elonperiod_start_pos, elonperiod_end_neg,
                     elonperiod_end_pos, indstart=1):
    """
    Elongation plotting, with periods of extended elongation included.

    Parameters
    ----------
    dpos_len : list
        Time derivative of positive ribbon elongation.
    dneg_len : list
        Time derivative of negative ribbon elongation.
    times : list
        Times corresponding to each AIA time step.
    times1600 : list
        Times corresponding to selected flare from 1600 Angstrom data.
    lens_pos_Mm : list
        Perpendicular extent of positive ribbon for each time step, in Mm.
    lens_neg_Mm : list
        Parallel extent of positive ribbon for each time step in Mm.
    flnum : int
        Flare index from RibbonDB database.
    elonperiod_start_pos : list
        Determined start times for elongation in positive ribbon.
    elonperiod_end_pos : list
        Determined end times for elongation in positive ribbon.
    elonperiod_start_neg : list
        Determined start times for elongation in negative ribbon.
    elonperiod_end_neg : list
        Determined end times for elongation in negative ribbon.
    indstart : int
        Start index for plotting. The default is 1.

    Returns
    -------
    None.

    """
    timelab = np.linspace(0, 24 * len(times), len(times))
    fig, [ax1, ax2, ax3] = plt.subplots(3, 1, figsize=(13, 20))
    ax3.plot(timelab[indstart:-1], dpos_len[indstart:], '-+', c='red',
             markersize=10, label='+ Ribbon')
    ax3.plot(timelab[indstart:-1], dneg_len[indstart:], '-+', c='blue',
             markersize=10, label='- Ribbon')
    ax3.legend(fontsize=15)
    ax3.grid()

    s = str(times1600[0])

    ax3.set_xlabel('Time [s since ' + s[2:13] + ', ' + s[13:-5] + ']',
                   font='Times New Roman', fontsize=17)
    ax3.set_ylabel('Elongation Rate [Mm/sec]', font='Times New Roman',
                   fontsize=17)
    ax3.set_title('Ribbon Elongation Rate', font='Times New Roman',
                  fontsize=45)

    ax1.plot(timelab[indstart:-1], lens_pos_Mm[indstart:-1], '-o', c='red',
             markersize=6, label='mean')
    ax2.plot(timelab[indstart:-1], lens_neg_Mm[indstart:-1], '-o', c='blue',
             markersize=6, label='mean')

    ax1.grid()
    ax1.set_ylabel('Distance [Mm]', font='Times New Roman', fontsize=17)
    ax1.set_title('Ribbon Elongation, Positive Ribbon',
                  font='Times New Roman', fontsize=45)
    ax2.set_ylabel('Distance [Mm]', font='Times New Roman', fontsize=17)
    ax2.set_title('Ribbon Elongation, Negative Ribbon',
                  font='Times New Roman', fontsize=45)
    ax2.grid()

    for i, j in zip(elonperiod_start_pos, elonperiod_end_pos):
        ax1.axvline(timelab[i], c='green')
        ax1.axvline(timelab[j], c='red')
        ax1.axvspan(timelab[i], timelab[j], alpha=0.5, color='pink')
    for k, l in zip(elonperiod_start_neg, elonperiod_end_neg):
        ax2.axvline(timelab[k], c='green')
        ax2.axvline(timelab[l], c='red')
        ax2.axvspan(timelab[k], timelab[l], alpha=0.5, color='cyan')

    fig.savefig(str(flnum) + 'elon_timing_plt.png')

    return None


def sep_period_plot(dpos_dist, dneg_dist, times, distpos_Mm, distneg_Mm, flnum,
                    sepperiod_start_pos, sepperiod_end_pos,
                    sepperiod_start_neg, sepperiod_end_neg, indstrt):
    """
    Separation plots, including periods of extended perpendicular motion.

    Parameters
    ----------
    dpos_dist : list
        Time derivative of positive ribbon separation.
    dneg_dist : list
        Time derivative of negative ribbon separation.
    times : arr
        Array of times for the flare, from AIA datafile.
    distpos_Mm : list
        Perpendicular extent of negative ribbon for each time step, in Mm.
    distneg_Mm : list
        Parallel extent of positive ribbon for each time step in Mm.
    flnum : int
        Flare index from RibbonDB database.
    sepperiod_start_pos : list
        Start times for periods of extended separation, positive ribbon.
    sepperiod_end_pos : list
        End times for periods of extended separation, positive ribbon.
    sepperiod_start_neg : list
        Start times for periods of extended separation, negative ribbon.
    sepperiod_end_neg : list
        End times for periods of extended separation, negative ribbon.
    indstart : int
        Start index for plotting. The default is 1.

    Returns
    -------
    None.

    """
    timelab = range(0, 24 * len(times), 24)

    fig, [ax1, ax2, ax3] = plt.subplots(3, 1, figsize=(13, 20))
    ax3.plot(timelab[indstrt:-1], scipy.signal.medfilt(dpos_dist[indstrt:],
                                                       kernel_size=3), '-+',
             c='red', markersize=10, label='+ Ribbon')
    ax3.plot(timelab[indstrt:-1], scipy.signal.medfilt(dneg_dist[indstrt:],
                                                       kernel_size=3), '-+',
             c='blue', markersize=10, label='- Ribbon')
    ax3.legend(fontsize=15)
    ax3.grid()

    s = str(times[0])

    ax3.set_xlabel('Time [s since ' + s[2:12] + ', ' + s[13:-5] + ']',
                   font='Times New Roman', fontsize=17)
    ax3.set_ylabel('Separation Rate [Mm/sec]', font='Times New Roman',
                   fontsize=17)
    ax3.set_title('Ribbon Separation Rate', font='Times New Roman',
                  fontsize=45)

    ax1.plot(timelab[indstrt:-1], distpos_Mm[indstrt:-1], '-o', c='red',
             markersize=6, label='mean')
    ax2.plot(timelab[indstrt:-1], distneg_Mm[indstrt:-1], '-o', c='blue',
             markersize=6, label='mean')
    ax1.grid()
    ax1.set_ylabel('Distance [Mm]', font='Times New Roman', fontsize=17)
    ax1.set_title('Ribbon Separation, Positive Ribbon',
                  font='Times New Roman', fontsize=45)
    ax2.set_ylabel('Distance [Mm]', font='Times New Roman', fontsize=17)
    ax2.set_title('Ribbon Separation, Negative Ribbon',
                  font='Times New Roman', fontsize=45)
    ax2.grid()

    for i, j in zip(sepperiod_start_pos, sepperiod_end_pos):
        ax1.axvline(timelab[i], c='green')
        ax1.axvline(timelab[j], c='red')
        ax1.axvspan(timelab[i], timelab[j], alpha=0.5, color='pink')
    for k, l in zip(sepperiod_start_neg, sepperiod_end_neg):
        ax2.axvline(timelab[k], c='green')
        ax2.axvline(timelab[l], c='red')
        ax2.axvspan(timelab[k], timelab[l], alpha=0.5, color='cyan')

    fig.savefig(str(flnum) + 'sep_timing_plt.png')

    return None


def flux_rec_mod_process(sav_data, dt1600, pos1600, neg1600):
    """
    Process data in order to produce reconnection flux and rate arrays in later
    functions.

    Parameters
    ----------
    sav_data : AttrDict
        Dictionary containing all of the saved parameters in the HMI file.
    dt1600 : list
        1600 Angstrom datetime values for flare.
    pos1600 : list
        Summed pixel numbers in positive ribbon for each time step.
    neg1600 : list
        Summed pixel numbers in negative ribbon for each time step.

    Returns
    -------
    hmi : list
        HMI data from file.
    aia8_pos : list
        Cumulative positive ribbon pixels for AIA, c=8.
    aia8_neg : list
        Cumulative negative ribbon pixels for AIA, c=8.
    aia8_inst_pos : list
        Instantaneous positive ribbon pixels for AIA, c=8.
    aia8_inst_neg : list
        Instantaneous negative ribbon pixels for AIA, c=8.
    peak_pos : list
        Time of peak counts from 1600 Angstrom data in positive ribbon.
    peak_neg : list
        Time of peak counts from 1600 Angstrom data in negative ribbon.

    """
    # Process data for reconnection flux, reconnection rate,
    # rise phase exponential modeling
    hmi = sav_data.hmi
    aia8 = sav_data.pos8
    aia8_inst = sav_data.inst_pos8
    aia8_pos = np.zeros(np.shape(aia8))
    aia8_neg = np.zeros(np.shape(aia8))
    aia8_inst_pos = np.zeros(np.shape(aia8_inst))
    aia8_inst_neg = np.zeros(np.shape(aia8_inst))
    xsh, ysh, zsh = aia8.shape
    hmi_dat = sav_data.hmi

    for i, j, k in np.ndindex(aia8.shape):
        if aia8[i, j, k] == 1 and hmi_dat[j, k] > 0:
            aia8_pos[i, j, k] = 1
        elif aia8[i, j, k] == 1 and hmi_dat[j, k] < 0:
            aia8_neg[i, j, k] = 1

    for i, j, k in np.ndindex(aia8.shape):
        if aia8_inst[i, j, k] == 1 and hmi_dat[j, k] > 0:
            aia8_inst_pos[i, j, k] = 1
        elif aia8_inst[i, j, k] == 1 and hmi_dat[j, k] < 0:
            aia8_inst_neg[i, j, k] = 1

    peak_pos = dt1600[np.argmax(pos1600)]
    peak_neg = dt1600[np.argmax(neg1600)]

    return hmi, aia8_pos, aia8_neg, aia8_inst_pos, aia8_inst_neg, peak_pos, \
        peak_neg


def inst_flux_process(aia8_inst_pos, aia8_inst_neg, flnum, conv_f, hmi,
                      dt1600, peak_pos, peak_neg):
    """
    Find and plot instantaneous reconnection flux values.

    Parameters
    ----------
    aia8_inst_pos : list
        Instantaneous positive ribbon pixels for AIA, c=8.
    aia8_inst_neg : list
        Instantaneous negative ribbon pixels for AIA, c=8.
    flnum : int
        Flare index from RibbonDB database.
    conv_f : float
        Conversion factor from pixels to Mm.
    hmi : list
        HMI data for flare in question.
    dt1600 : list
        Datetime values for 1600 Angstrom data.
    peak_pos : list
        Time of peak counts from 1600 Angstrom data in positive ribbon.
    peak_neg : list
        Time of peak counts from 1600 Angstrom data in negative ribbon.

    Returns
    -------
    rec_flux_pos_inst : list
        Instantaneous reconnection flux, positive ribbon.
    rec_flux_neg_inst : list
        Instantaneous reconnection flux, negative ribbon.
    pos_pix_inst : list
        Instantaneous ribbon pixel counts, positive ribbon.
    neg_pix_inst : list
        Instantaneous ribbon pixel counts, negative ribbon.
    ds2 : float
        Conversion factor for 2D size of pixel.

    """
    rec_flux_pos_inst = np.zeros(len(aia8_inst_pos))
    rec_flux_neg_inst = np.zeros(len(aia8_inst_neg))
    pos_area_pix_inst = np.zeros(len(aia8_inst_pos))
    neg_area_pix_inst = np.zeros(len(aia8_inst_neg))
    pos_pix_inst = np.zeros(len(aia8_inst_pos))
    neg_pix_inst = np.zeros(len(aia8_inst_neg))

    conv_f_cm = conv_f * 1e6 * 100  # conversion factor in cm
    ds2 = conv_f_cm**2  # for each pixel grid

    for i in range(len(aia8_inst_pos)):
        pos_mask_inst = aia8_inst_pos[i, :, :]
        neg_mask_inst = aia8_inst_neg[i, :, :]

        pos_area_pix_inst[i] = np.sum(pos_mask_inst)
        neg_area_pix_inst[i] = np.sum(neg_mask_inst)

        hmi_pos_inst = pos_mask_inst*hmi  # in G
        hmi_neg_inst = neg_mask_inst*hmi  # in G

        pos_pix_inst[i] = np.sum(hmi_pos_inst)
        neg_pix_inst[i] = np.sum(hmi_neg_inst)
        rec_flux_pos_inst[i] = np.sum(hmi_pos_inst) * ds2
        rec_flux_neg_inst[i] = np.sum(hmi_neg_inst) * ds2

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(dt1600, rec_flux_pos_inst, c='red', label='+')
    ax.scatter(dt1600, rec_flux_neg_inst, c='blue', label='-')
    ax.grid()
    ax.set_xlabel('Time', font='Times New Roman', fontsize=20)
    ax.axvline(peak_pos, c='red', linestyle=':')
    ax.axvline(peak_neg, c='blue', linestyle='-.')
    ax.set_ylabel('Reconnection Flux [Mx]', font='Times New Roman',
                  fontsize=20)
    ax.set_title('Reconnection Flux', font='Times New Roman', fontsize=45)
    ax.legend()

    fig.savefig(str(flnum) + '_inst_flx.png')

    return rec_flux_pos_inst, rec_flux_neg_inst, pos_pix_inst, neg_pix_inst,\
        ds2


def cumul_flux_process(aia8_pos, aia8_neg, conv_f, flnum, peak_pos, peak_neg,
                       hmi, dt1600):
    """
    Determine reconnection flux from cumulative ribbon masks.

    Parameters
    ----------
    aia8_pos : list
        Cumulative positive ribbon pixels for AIA, c=8.
    aia8_neg : list
        Cumulative negative ribbon pixels for AIA, c=8.
    flnum : int
        Flare index from RibbonDB database.
    conv_f : float
        Conversion factor from pixels to Mm.
    peak_pos : list
        Time of peak counts from 1600 Angstrom data in positive ribbon.
    peak_neg : list
        Time of peak counts from 1600 Angstrom data in negative ribbon.
    hmi : list
        HMI data for flare in question.
    dt1600 : list
        Datetime values for 1600 Angstrom data.

    Returns
    -------
    rec_flux_pos : list
        Cumulative reconnection flux, positive ribbon.
    rec_flux_neg : list
        Cumulative reconnection flux, negative ribbon.
    pos_pix : list
        Ribbon counts of cumulative masks, positive ribbon.
    neg_pix : list
        Ribbon counts of cumulative masks, negative ribbon.
    pos_area_pix : list
        Positive cumulative ribbon area.
    neg_area_pix : list
        Negative cumulative ribbon area.
    ds2 : float
        Conversion factor for 2D size of pixel.
    pos_area : list
        Positive cumulative ribbon area, Mm.
    neg_area : list
        Negative cumulative ribbon area, Mm.

    """
    rec_flux_pos = np.zeros(len(aia8_pos))
    rec_flux_neg = np.zeros(len(aia8_neg))
    pos_area_pix = np.zeros(len(aia8_pos))
    neg_area_pix = np.zeros(len(aia8_neg))
    pos_pix = np.zeros(len(aia8_pos))
    neg_pix = np.zeros(len(aia8_neg))
    pos_area = pos_area_pix
    neg_area = neg_area_pix

    conv_f_cm = conv_f * 1e6 * 100  # conversion factor in cm
    ds2 = conv_f_cm**2
    for i in range(len(aia8_pos)):
        pos_mask = aia8_pos[i, :, :]
        neg_mask = aia8_neg[i, :, :]

        pos_area_pix[i] = np.sum(pos_mask)
        neg_area_pix[i] = np.sum(neg_mask)

        hmi_pos = pos_mask*hmi  # in G
        hmi_neg = neg_mask*hmi  # in G

        pos_pix[i] = np.sum(hmi_pos)
        neg_pix[i] = np.sum(hmi_neg)
        rec_flux_pos[i] = np.sum(hmi_pos) * ds2
        rec_flux_neg[i] = np.sum(hmi_neg) * ds2

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(dt1600, rec_flux_pos, c='red', label='+')
    ax.scatter(dt1600, rec_flux_neg, c='blue', label='-')
    ax.grid()
    ax.set_xlabel('Time', font='Times New Roman', fontsize=20)
    ax.axvline(peak_pos, c='red', linestyle=':')
    ax.axvline(peak_neg, c='blue', linestyle='-.')
    ax.set_ylabel('Reconnection Flux [Mx]', font='Times New Roman',
                  fontsize=20)
    ax.set_title('Reconnection Flux', font='Times New Roman', fontsize=45)
    ax.legend()

    fig.savefig(str(flnum) + '_cumul_flx.png')

    return rec_flux_pos, rec_flux_neg, pos_pix, neg_pix, pos_area_pix, \
        neg_area_pix, ds2, pos_area, neg_area


def exp_curve_fit(exp_ind, exp_ind_area, rec_flux_pos, rec_flux_neg,
                  exponential, exponential_neg, pos_area, neg_area):
    """
    Fit exponential curve to flux and ribbon area curves for each ribbon.

    Parameters
    ----------
    exp_ind : int
        The index where to stop the exponential fitting.
    rec_flux_pos : list
        Cumulative reconnection flux, positive ribbon.
    rec_flux_neg : list
        Cumulative reconnection flux, negative ribbon.
    exponential : function
        Exponential function definition.
    exponential_neg : function
        Negative exponential function definition.
    pos_area : list
        Positive cumulative ribbon area, Mm.
    neg_area : list
        Negative cumulative ribbon area, Mm.

    Returns
    -------
    poptposflx : list
        Fitting parameters, positive reconnection flux.
    pcovposflx : list
        Covariance matrix entries, positive reconnection flux.
    poptnegflx : list
        Fitting parameters, negative reconnection flux.
    pcovnegflx : list
        Covariance matrix entries, negative reconnection flux.
    poptpos : list
        Fitting parameters, positive ribbon area.
    poptneg : list
        Fitting parameters, negative ribbon area.
    pcovpos : list
        Covariance matrix entries, positive ribbon area.
    pcovneg : list
        Covariance matrix entries, negative ribbon area.
    rise_pos_flx : list
        Rise phase reconnection flux, positive ribbon.
    rise_neg_flx : list
        Rise phase reconnection flux, negative ribbon.

    """

    # Fit only from start to specified exp_ind; usually the index corresponding
    # to the peak of the light curve, but sometimes not.
    rise_pos_flx = rec_flux_pos[0:exp_ind]
    rise_neg_flx = rec_flux_neg[0:exp_ind]
    rise_pos_area = pos_area[0:exp_ind_area]
    rise_neg_area = neg_area[0:exp_ind_area]

    # Fitting to exponential and negative exponential models
    poptposflx, pcovposflx = \
        scipy.optimize.curve_fit(exponential, range(0, len(rise_pos_flx)),
                                 rise_pos_flx)
    poptnegflx, pcovnegflx = \
        scipy.optimize.curve_fit(exponential_neg, range(0, len(rise_neg_flx)),
                                 rise_neg_flx)
    poptpos, pcovpos = \
        scipy.optimize.curve_fit(exponential, range(0, len(rise_pos_area)),
                                 rise_pos_area)
    poptneg, pcovneg = \
        scipy.optimize.curve_fit(exponential, range(0, len(rise_neg_area)),
                                 rise_neg_area)

    return poptposflx, pcovposflx, poptnegflx, pcovnegflx, poptpos, poptneg, \
        pcovpos, pcovneg, rise_pos_flx, rise_neg_flx


def exp_curve_plt(dt1600, rec_flux_pos, rec_flux_neg, rise_pos_flx,
                  rise_neg_flx, peak_pos, peak_neg, exp_ind, ds2, exponential,
                  exponential_neg, poptposflx, poptnegflx, flnum):
    """
    Plotting of exponential curve with fit.

    Parameters
    ----------
    dt1600 : list
        Datetime values for 1600 Angstrom data.
    rec_flux_pos : list
        Cumulative reconnection flux, positive ribbon.
    rec_flux_neg : list
        Cumulative reconnection flux, negative ribbon.
    rise_pos_flx : list
        Rise phase reconnection flux, positive ribbon.
    rise_neg_flx : list
        Rise phase reconnection flux, negative ribbon.
    peak_pos : list
        Time of peak counts from 1600 Angstrom data in positive ribbon.
    peak_neg : list
        Time of peak counts from 1600 Angstrom data in negative ribbon.
    exp_ind : int
        The index where to stop the exponential fitting.
    ds2 : float
        Conversion factor for 2D size of pixel.
    exponential : function
        Exponential function definition.
    exponential_neg : function
        Negative exponential function definition.
    poptposflx : list
        Fitting parameters, positive reconnection flux.
    poptnegflx : list
        Fitting parameters, negative reconnection flux.
    flnum : int
        Flare index from RibbonDB database.

    Returns
    -------
    None.

    """

    rise_pos_time = dt1600[0:exp_ind]
    rise_neg_time = dt1600[0:exp_ind]

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(dt1600, rec_flux_pos, c='red', label='+')
    ax.scatter(dt1600, rec_flux_neg, c='blue', label='-')
    ax.grid()
    ax.set_xlabel('Time', font='Times New Roman', fontsize=20)
    ax.axvline(peak_pos, c='red', linestyle=':')
    ax.axvline(peak_neg, c='blue', linestyle='-.')
    ax.set_ylabel('Reconnection Flux [Mx]',
                  font='Times New Roman', fontsize=20)
    ax.set_title('Reconnection Flux', font='Times New Roman', fontsize=45)
    ax.plot(rise_pos_time, ds2*exponential(range(0, len(rise_pos_flx)),
                                           *poptposflx), 'r-',
            label='Exponential Model, +')
    ax.plot(rise_neg_time, ds2 * exponential_neg(range(0, len(rise_neg_flx)),
                                                 *poptnegflx), 'b-',
            label='Exponential Model, -')
    ax.axvline(dt1600[exp_ind])
    ax.legend()

    fig.savefig(str(flnum) + '_recflux_model.png')

    # Now plot log-log of just the impulsive phase

    fig2, [ax1, ax2] = plt.subplots(2, 1, figsize=(10, 20))
    ax1.scatter((dt1600), np.log(rec_flux_pos), c='red')
    ax2.scatter((dt1600), -np.log(- rec_flux_neg), c='blue')
    ax1.grid()
    ax2.grid()
    ax1.set_xlabel('Time', font='Times New Roman', fontsize=20)
    ax2.set_xlabel('Time', font='Times New Roman', fontsize=20)

    ax1.plot((rise_pos_time),
             np.log(ds2 * exponential(range(0, len(rise_pos_flx)),
                                      *poptposflx)), 'r-',
             label='Exponential Model, +')

    ax2.plot((rise_neg_time),
             -np.log(-ds2 * exponential_neg(range(0, len(rise_neg_flx)),
                                            *poptnegflx)), 'b-',
             label='Exponential Model, -')

    ax1.set_ylabel(r'Rec. Flx [Mx]', font='Times New Roman', fontsize=20)
    ax1.set_title('Reconnection Flux, Impulsive Phase',
                  font='Times New Roman', fontsize=45)
    ax1.set_xlim(dt1600[0], dt1600[exp_ind])
    ax1.legend(fontsize=15)
    ax2.set_ylabel(r'Rec. Flx [Mx]', font='Times New Roman', fontsize=20)
    ax2.set_title('Reconnection Flux, Impulsive Phase',
                  font='Times New Roman', fontsize=45)
    ax2.set_xlim(dt1600[0], dt1600[exp_ind])
    ax2.legend(fontsize=15)

    fig2.savefig(str(flnum) + '_rec_impphase_model.png')

    return None


def rib_area_plt(dt1600, poptpos, poptneg, flnum, pos_area_pix, neg_area_pix,
                 peak_pos, peak_neg, exp_ind):
    """
    Plotting ribbon areas with fitted models.

    Parameters
    ----------
    dt1600 : list
        Datetime values for 1600 Angstrom data.
    poptpos : list
        Fitting parameters, positive ribbon area.
    poptneg : list
        Fitting parameters, negative ribbon area.
    flnum : int
        Flare index from RibbonDB database.
    pos_area_pix : list
        Positive cumulative ribbon area.
    neg_area_pix : list
        Negative cumulative ribbon area.
    peak_pos : list
        Time of peak counts from 1600 Angstrom data in positive ribbon.
    peak_neg : list
        Time of peak counts from 1600 Angstrom data in negative ribbon.
    exp_ind : int
        The index where to stop the exponential fitting.

    Returns
    -------
    None.

    """
    # Cumulative
    pos_area = pos_area_pix
    neg_area = neg_area_pix
    rise_pos_area = pos_area[0:exp_ind]
    rise_neg_area = neg_area[0:exp_ind]

    # Plot just the ribbon areas, c=8
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(dt1600, pos_area, c='red', label='+')
    ax.scatter(dt1600, neg_area, c='blue', label='-')
    rise_pos_time = dt1600[0:exp_ind]
    rise_neg_time = dt1600[0:exp_ind]
    ax.grid()
    ax.set_xlabel('Time', font='Times New Roman', fontsize=20)
    ax.axvline(peak_pos, c='red', linestyle=':')
    ax.axvline(peak_neg, c='blue', linestyle='-.')
    ax.plot(rise_pos_time, exponential(range(0, len(rise_pos_area)), *poptpos),
            'r-', label='Exponential Model, +')
    ax.plot(rise_neg_time, exponential(range(0, len(rise_neg_area)), *poptneg),
            'b-', label='Exponential Model, -')
    ax.set_ylabel('Ribbon Area [cm^2]', font='Times New Roman', fontsize=20)
    ax.set_title('Ribbon Area', font='Times New Roman', fontsize=45)

    # If end of modeling region is before end of impulsive phase
    ax.axvline(dt1600[exp_ind])
    ax.legend()

    fig.savefig(str(flnum) + '_ribarea_model.png')

    # Just impulsive region, with log-log
    fig2, [ax1, ax2] = plt.subplots(2, 1, figsize=(10, 20))
    ax1.scatter((dt1600), np.log(pos_area), c='red')
    ax2.scatter((dt1600), np.log(neg_area), c='blue')
    ax1.grid()
    ax2.grid()
    ax1.set_xlabel('Time', font='Times New Roman', fontsize=20)
    ax2.set_xlabel('Time', font='Times New Roman', fontsize=20)
    ax1.plot((rise_pos_time), np.log(exponential(range(0, len(rise_pos_area)),
                                                 *poptpos)), 'r-',
             label='Exponential Model, +')
    ax2.plot((rise_neg_time), np.log(exponential(range(0, len(rise_neg_area)),
                                                 *poptneg)), 'b-',
             label='Exponential Model, -')

    ax1.set_ylabel('Ribbon Area [cm^2]', font='Times New Roman',
                   fontsize=20)
    ax1.set_title('Ribbon Area, Impulsive Phase', font='Times New Roman',
                  fontsize=45)
    ax1.set_xlim(dt1600[0], dt1600[exp_ind])
    ax1.legend(fontsize=15)
    ax2.set_ylabel('Ribbon Area [cm^2]', font='Times New Roman', fontsize=20)
    ax2.set_title('Ribbon Area, Impulsive Phase', font='Times New Roman',
                  fontsize=45)
    ax2.set_xlim(dt1600[0], dt1600[exp_ind])
    ax2.legend(fontsize=15)

    fig2.savefig(str(flnum) + '_impphase_model.png')

    return None


def rec_rate(rec_flux_pos, rec_flux_neg, dn1600, dt1600, peak_pos, peak_neg,
             flnum):
    """
    Reconnection rate determination from reconnection flux values.

    Parameters
    ----------
    rec_flux_pos : list
        Cumulative reconnection flux, positive ribbon.
    rec_flux_neg : list
        Cumulative reconnection flux, negative ribbon.
    dn1600 : list
        Datenum values for 1600 Angstrom data.
    dt1600 : list
        Datetime values for 1600 Angstrom data.
    peak_pos : list
        Time of peak counts from 1600 Angstrom data in positive ribbon.
    peak_neg : list
        Time of peak counts from 1600 Angstrom data in negative ribbon.
    flnum : int
        Flare index from RibbonDB database.

    Returns
    -------
    rec_rate_pos : list
        Reconnection rates for positive ribbon.
    rec_rate_neg : list
        Reconnection rates for negative ribbon.

    """
    rec_rate_pos = (np.diff(rec_flux_pos) / (dn1600[1] - dn1600[0]))\
        / 3600 / 24  # Mx/s
    rec_rate_neg = (np.diff(rec_flux_neg) / (dn1600[1] - dn1600[0]))\
        / 3600 / 24  # Mx/s

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(dt1600[1:], rec_rate_pos, c='red', label='+')
    ax.scatter(dt1600[1:], rec_rate_neg, c='blue', label='-')
    ax.grid()
    ax.set_xlabel('Time', font='Times New Roman', fontsize=20)
    ax.axvline(peak_pos, c='red', linestyle=':')
    ax.axvline(peak_neg, c='blue', linestyle='-.')
    ax.set_ylabel('Reconnection Rate [Mx/s]', font='Times New Roman',
                  fontsize=20)
    ax.set_title('Reconnection Rate', font='Times New Roman', fontsize=45)

    fig.savefig(str(flnum) + '_recrate.png')

    return rec_rate_pos, rec_rate_neg


def shear_ribbon_isolation(aia8_neg, aia8_pos, med_x, med_y,
                           pt_range=[-2, -1, 1, 2], poscrit=6, negcrit=6,
                           negylow=400, negyhi=0, negxlow=300,
                           negxhi=400, posylow=0, posyhi=0,
                           posxlow=450, posxhi=0, flag=0):
    """
    Isolates ribbons with the shear algorithm in mind.

    Parameters
    ----------
    aia8_neg : arr
        Negative polarity array from AIA.
    aia8_pos : arr
        Positive polarity array from AIA
    med_x : int
        Median of x dimension.
    med_y : int
        Median of y dimension.
    pt_range : arr, optional
        Range of points to search around each pixel. The default is
        [-2,-1,1,2].
    poscrit : int, optional
        Minimum number of positive points a pixel must be surrounded by in
        order to be counted. The default is 6.
    negcrit : int, optional
        Minimum number of negative points a pixel must be surrounded by in
        order to be counted. The default is 6.
    negylow : int, optional
        Low y-dimension of image to search through, negative ribbon. The
        default is 400.
    negyhi : int, optional
        High y-dimension of image to search through, negative ribbon.
        The default is 0.
    negxlow : int, optional
        Low x-dimension of image to search through, negative ribbon. The
        default is 300.
    negxhi : int, optional
        High x-dimension of image to search through, negative ribbon.
        The default is 400.
    posylow : int, optional
        Low y-dimension of image to search through, positive ribbon. The
        default is 0.
    posyhi : int, optional
        High y-dimension of image to search through, positive ribbon. The
        default is 0.
    posxlow : int, optional
        Low x-dimension of image to search through, positive ribbon. The
        default is 450.
    posxhi : int, optional
        High x-dimension of image to search through, positive ribbon. The
        default is 0.
    flag : int, optional
        If 0, values for limits in x and y are set based on the median values.
        If 1, values for limits in x and y are set via the function inputs.

    Returns
    -------
    aia_neg_rem_shear : arr
        AIA 1600 Angstrom images of negative ribbon, with spur pixels removed.
    aia_pos_rem_shear : arr
        AIA 1600 Angstrom images of positive ribbon, with spur pixels removed.

    """
    neg_rem_shear = np.zeros(np.shape(aia8_pos))
    pos_rem_shear = np.zeros(np.shape(aia8_neg))
    aia_pos_rem_shear = np.zeros(np.shape(aia8_pos))
    aia_neg_rem_shear = np.zeros(np.shape(aia8_neg))

    if flag == 0:
        negylow = int(round(med_y) - 100)
        negyhi = int(round(med_y) + 100)
        negxlow = int(round(med_x) - 100)
        negxhi = int(round(med_x) + 100)

        posylow = int(round(med_y) - 100)
        posyhi = int(round(med_y) + 100)
        posxlow = int(round(med_x) - 100)
        posxhi = int(round(med_x) + 100)

    # Search through negative ribbon and remove spur points.
    for i in range(len(neg_rem_shear)):
        for j in range(len(neg_rem_shear[0])-2):
            for k in range(len(neg_rem_shear[1])-2):
                n = 0
                if aia8_neg[i, j, k] == 1:
                    for h in pt_range:
                        for m in pt_range:
                            if aia8_neg[i, j+h, k+m] == 1:
                                n = n + 1
                    if (n > negcrit):
                        neg_rem_shear[i, j, k] = 1
                    else:
                        neg_rem_shear[i, j, k] = 0
                else:
                    neg_rem_shear[i, j, k] = 0

    # Search through positive ribbon and remove spur points.
    for i in range(len(pos_rem_shear)):
        for j in range(len(pos_rem_shear[0])-2):
            for k in range(len(pos_rem_shear[1])-2):
                n = 0
                if aia8_pos[i, j, k] == 1:
                    for h in pt_range:
                        for m in pt_range:
                            if aia8_pos[i, j+h, k+m] == 1:
                                n = n + 1
                    if (n > poscrit):
                        pos_rem_shear[i, j, k] = 1
                    else:
                        pos_rem_shear[i, j, k] = 0
                else:
                    pos_rem_shear[i, j, k] = 0

    # Limit ribbons to within desired range for analysis
    for i in range(len(aia8_neg)):
        for j in range(negylow, negyhi):
            for k in range(negxlow, negxhi):
                if neg_rem_shear[i, j, k] > 0:
                    aia_neg_rem_shear[i, j, k] = 1

    for i in range(len(aia8_pos)):
        for j in range(posylow, posyhi):
            for k in range(posxlow, posxhi):
                if pos_rem_shear[i, j, k] > 0:
                    aia_pos_rem_shear[i, j, k] = 1

    return aia_neg_rem_shear, aia_pos_rem_shear

# find left and rightmost pixels


def leftrightshear(aia_pos_rem_shear, aia_neg_rem_shear):
    """
    Finds leftmost and rightmost pixels of each ribbon.

    Parameters
    ----------
    aia_pos_rem_shear : arr
        Images of isolated positive ribbons from above function.
    aia_neg_rem_shear : arr
        Images of isolated negative ribbons from above function.

    Returns
    -------
    lr_coord_neg_shear : arr
        Left and right coordinates of negative ribbon.
    lr_coord_pos_shear : arr
        Left and right coordinates of positive ribbon..

    """

    lr_coord_pos_shear = np.zeros([len(aia_pos_rem_shear), 4])
    lr_coord_neg_shear = np.zeros([len(aia_neg_rem_shear), 4])

    # Find left and rightmost pixels for positive ribbon in each frame
    for i in range(len(aia_pos_rem_shear)):
        left_x = 0
        left_y = 0
        right_x = 0
        right_y = 0

        for k in range(len(aia_pos_rem_shear[1])):
            for j in range(len(aia_pos_rem_shear[0])):
                if aia_pos_rem_shear[i, j, k] == 1:
                    left_x = k
                    left_y = j
                    break

            if left_x != 0:
                break

        for k in range(len(aia_pos_rem_shear[1])-1, 0, -1):
            for j in range(len(aia_pos_rem_shear[0])):
                if aia_pos_rem_shear[i, j, k] == 1:
                    right_x = k
                    right_y = j
                    break
            if right_x != 0:
                break

        lr_coord_pos_shear[i, :] = [left_x, left_y, right_x, right_y]

    # The same, for negative ribbon
    for i in range(len(aia_neg_rem_shear)):
        left_x = 0
        left_y = 0
        right_x = 0
        right_y = 0

        for k in range(len(aia_neg_rem_shear[1])):
            for j in range(len(aia_neg_rem_shear[0])):
                if aia_neg_rem_shear[i, j, k] == 1:
                    left_x = k
                    left_y = j
                    break

            if left_x != 0:
                break

        for k in range(len(aia_neg_rem_shear[1])-1, 0, -1):
            for j in range(len(aia_neg_rem_shear[0])):
                if aia_neg_rem_shear[i, j, k] == 1:
                    right_x = k
                    right_y = j
                    break
            if right_x != 0:
                break

        lr_coord_neg_shear[i, :] = [left_x, left_y, right_x, right_y]

    return lr_coord_neg_shear, lr_coord_pos_shear


def sheardists(lr_coord_pos_shear, lr_coord_neg_shear, ivs_sort, dvs_sort):
    """
    Find the position of the PIL nearest to the extremes of each ribbon.

    Parameters
    ----------
    lr_coord_pos_shear : arr
        Left and rightmost pixels for the positive ribbon.
    lr_coord_neg_shear : arr
        Left and rightmost pixels for negative ribbon.
    ivs_sort : arr
        Independent variable indices for PIL.
    dvs_sort : arr
        Dependent variable indices for PIL.

    Returns
    -------
    pil_right_near_pos_shear : arr
        Indices of PIL point nearest positive ribbon, rightmost extreme.
    pil_left_near_pos_shear : arr
        Indices of PIL points nearest positive ribbon, leftmost extreme.
    pil_right_near_neg_shear : arr
        Indices of PIL points nearest negative ribbon, rightmost extreme.
    pil_left_near_neg_shear : arr
        Indices of PIL points nearest negative ribbon, leftmost extreme.

    """
    left_pil_dist_pos_shear = np.zeros(
        [len(lr_coord_pos_shear), len(ivs_sort)])
    right_pil_dist_pos_shear = np.zeros([len(lr_coord_pos_shear),
                                         len(ivs_sort)])

    pil_left_near_pos_shear = np.zeros([len(left_pil_dist_pos_shear), 3])
    pil_right_near_pos_shear = np.zeros([len(right_pil_dist_pos_shear), 3])

    left_pil_dist_neg_shear = np.zeros(
        [len(lr_coord_neg_shear), len(ivs_sort)])
    right_pil_dist_neg_shear = np.zeros([len(lr_coord_neg_shear),
                                         len(ivs_sort)])
    pil_left_near_neg_shear = np.zeros([len(left_pil_dist_neg_shear), 3])
    pil_right_near_neg_shear = np.zeros([len(right_pil_dist_neg_shear), 3])

    # Arrays of distances from all positive ribbon points to all PIL points, in
    # each dimension
    for i in range(len(lr_coord_pos_shear)):
        left_x, left_y, right_x, right_y = lr_coord_pos_shear[i]
        for j in range(len(ivs_sort)):
            left_pil_dist_pos_shear[i, j] = np.sqrt((left_x - ivs_sort[j])**2 +
                                                    (left_y - dvs_sort[j])**2)
            right_pil_dist_pos_shear[i, j] = np.sqrt((right_x - ivs_sort[j])**2
                                                     + (right_y -
                                                        dvs_sort[j])**2)

    # Arrays of distances from all negative ribbon points to all PIL points, in
    # each dimension
    for i in range(len(lr_coord_neg_shear)):
        left_x, left_y, right_x, right_y = lr_coord_neg_shear[i]
        for j in range(len(ivs_sort)):
            left_pil_dist_neg_shear[i, j] = np.sqrt((left_x - ivs_sort[j])**2 +
                                                    (left_y - dvs_sort[j])**2)
            right_pil_dist_neg_shear[i, j] = np.sqrt((right_x - ivs_sort[j])**2
                                                     + (right_y -
                                                        dvs_sort[j])**2)

    # Find smallest distance, and the corresponding PIL point
    for i in range(len(left_pil_dist_pos_shear)):
        ind = np.where(left_pil_dist_pos_shear[i] ==
                       np.min(left_pil_dist_pos_shear[i]))
        pil_left_near_pos_shear[i, :] = [ivs_sort[ind[0][0]],
                                         dvs_sort[ind[0][0]], ind[0][0]]

    for j in range(len(right_pil_dist_neg_shear)):
        ind = np.where(right_pil_dist_neg_shear[j] ==
                       np.min(right_pil_dist_neg_shear[j]))
        pil_right_near_neg_shear[j, :] = [ivs_sort[ind[0][0]],
                                          dvs_sort[ind[0][0]], ind[0][0]]

    # Find smallest distance, and the corresponding PIL point
    for i in range(len(left_pil_dist_neg_shear)):
        ind = np.where(left_pil_dist_neg_shear[i] ==
                       np.min(left_pil_dist_neg_shear[i]))
        pil_left_near_neg_shear[i, :] = [ivs_sort[ind[0][0]],
                                         dvs_sort[ind[0][0]], ind[0][0]]

    for j in range(len(right_pil_dist_pos_shear)):
        ind = np.where(right_pil_dist_pos_shear[j] ==
                       np.min(right_pil_dist_pos_shear[j]))
        pil_right_near_pos_shear[j, :] = [ivs_sort[ind[0][0]],
                                          dvs_sort[ind[0][0]], ind[0][0]]

    return pil_right_near_pos_shear, pil_left_near_pos_shear,\
        pil_right_near_neg_shear, pil_left_near_neg_shear


def guidefieldlen(pil_right_near_pos_shear, pil_left_near_pos_shear,
                  pil_right_near_neg_shear, pil_left_near_neg_shear,
                  sortedpil, curve_length):
    """
    Find length along the axis of the guide field, the PIL-parallel component
    of magnetic field.

    Parameters
    ----------
    pil_right_near_pos_shear : arr
        Indices of PIL point nearest positive ribbon, rightmost extreme.
    pil_left_near_pos_shear : arr
        Indices of PIL points nearest positive ribbon, leftmost extreme.
    pil_right_near_neg_shear : arr
        Indices of PIL points nearest negative ribbon, rightmost extreme.
    pil_left_near_neg_shear : arr
        Indices of PIL points nearest negative ribbon, leftmost extreme.
    sortedpil : arr
        Independent and dependent values of PIL, sorted along PIL.
    curve_length : function
        Function for determination of a curve length between two points.

    Returns
    -------
    guide_right : arr
        Guide field strength, right-hand side.
    guide_left : arr
        Guide field strength, left-hand side.

    """
    guide_left = []
    guide_right = []

    # Length of segment of PIL between the leftmost ends of positive and
    # negative ribbons
    for i in range(len(pil_left_near_pos_shear)):
        posin = int(pil_left_near_pos_shear[i, 2])
        negin = int(pil_left_near_neg_shear[i, 2])
        if posin > negin:
            curvei = sortedpil[negin:posin, :]
        else:
            curvei = -sortedpil[posin:negin, :]
        guide_left.append(curve_length(curvei))

    # Length of segment of PIL between rightmost ends of positive and negative
    # ribbons
    for i in range(len(pil_right_near_pos_shear)):
        posin = int(pil_right_near_pos_shear[i, 2])
        negin = int(pil_right_near_neg_shear[i, 2])
        if posin > negin:
            curvei = sortedpil[negin:posin, :]
        else:
            curvei = -sortedpil[posin:negin, :]
        guide_right.append(curve_length(curvei))

    return guide_right, guide_left


def guidefieldlen_alt(pil_right_near_pos_shear, pil_left_near_pos_shear,
                      pil_right_near_neg_shear, pil_left_near_neg_shear,
                      sortedpil, curve_length, flag='posright'):
    """
    Find length along the axis of the guide field, the PIL-parallel component
    of magnetic field - alternative, tracking opposite ends of ribbons.

    Parameters
    ----------
    pil_right_near_pos_shear : arr
        Indices of PIL point nearest positive ribbon, rightmost extreme.
    pil_left_near_pos_shear : arr
        Indices of PIL points nearest positive ribbon, leftmost extreme.
    pil_right_near_neg_shear : arr
        Indices of PIL points nearest negative ribbon, rightmost extreme.
    pil_left_near_neg_shear : arr
        Indices of PIL points nearest negative ribbon, leftmost extreme.
    sortedpil : arr
        Independent and dependent values of PIL, sorted along PIL.
    flag : str, optional
        Ends of ribbons to check; either right positive/left negative or
        right negative/left positive.

    Returns
    -------
    guide : arr
        Guide field strength.

    """
    guide = []

    if flag == 'posleft':
        for i in range(len(pil_left_near_pos_shear)):
            posin = int(pil_left_near_pos_shear[i, 2])
            negin = int(pil_right_near_neg_shear[i, 2])
            if posin > negin:
                curvei = sortedpil[negin:posin, :]
            else:
                curvei = -sortedpil[posin:negin, :]
            guide.append(curve_length(curvei))
    elif flag == 'posright':
        for i in range(len(pil_right_near_pos_shear)):
            posin = int(pil_right_near_pos_shear[i, 2])
            negin = int(pil_left_near_neg_shear[i, 2])
            if posin > negin:
                curvei = sortedpil[negin:posin, :]
            else:
                curvei = -sortedpil[posin:negin, :]
            guide.append(curve_length(curvei))

    return guide


def gfrcalc(guide_left, guide_right, distneg_med, distpos_med):
    """
    Determines guide field ratio, the ratio of the PIL-parallel component of
    magnetic field to the PIL-perpendicular component, a proxy for the magnetic
    shear. Alternative, with cross-PIL-perpendicular loops.

    Parameters
    ----------
    guide : arr
        Guide-field strength.
    distneg_med : arr
        Distance from negative ribbon to PIL for each time step; from code for
        separation.
    distpos_med : arr
        Distance from positive ribbon to PIL for each time step; from code for
        separation.

    Returns
    -------
    gfr : arr
        Guide field ratio, left edge of ribbons.

    """
    left_gfr = guide_left/(distneg_med+distpos_med)
    right_gfr = guide_right/(distneg_med+distneg_med)

    return left_gfr, right_gfr


def gfrcalc_alt(guide, distneg_med, distpos_med):
    """
    Determines guide field ratio, the ratio of the PIL-parallel component of
    magnetic field to the PIL-perpendicular component, a proxy for the magnetic
    shear.

    Parameters
Iw     ----------
    guide_left : arr
        Guide field strength, right-hand side.
    guide_right : arr
        Guide field strength, left-hand side.
    distneg_med : arr
        Distance from negative ribbon to PIL for each time step; from code for
        separation.
    distpos_med : arr
        Distance from positive ribbon to PIL for each time step; from code for
        separation.

    Returns
    -------
    left_gfr : arr
        Guide field ratio, left edge of ribbons.
    right_gfr : arr
        Guide field ratio, right edge of ribbons.

    """
    gfr = guide/(distneg_med+distpos_med)

    return gfr


def plt_gfr(times, right_gfr, left_gfr, flnum, dt1600, flag=0):
    """
    Plots guide field ratio for right and left edges of ribbons.

    Parameters
    ----------
    times : arr
        Times corresponding to AIA 1600 Angstrom images.
    right_gfr : arr
        Guide field ratio, right edge of ribbons.
    left_gfr : arr
        Guide field ratio, left edge of ribbons.
    flnum : int
        Flare number from RibbonDB database.
    flag : int, optional
        0 if traditional shear determination, 1 if antiparallel.  Default is 0.
    Returns
    -------
    None.

    """
    timelab = range(0, 24*len(times), 24)
    s = str(dt1600[0])
    fig, ax = plt.subplots(figsize=(13, 7))
    if flag == 0:
        ax.plot(timelab, right_gfr, c='red', marker='o',
                label='GFR proxy, right')

        ax.plot(timelab, left_gfr, c='blue', marker='o',
                label='GFR proxy, left')
    elif flag == 1:
        ax.plot(timelab, right_gfr, c='blue', marker='o', label='GFR proxy')
        ax.set_xlabel('Time [s since '+s[5:-7]+']', font='Times New Roman',
                      fontsize=18)
        ax.set_ylabel('GFR Proxy', font='Times New Roman', fontsize=18)

    ax.set_title('Guide Field Ratio', font='Times New Roman', fontsize=20)
    ax.grid()
    ax.legend(fontsize=15)
    fig.savefig(str(flnum) + '_gfr.png')

    return None


def process_fermi(day, month, year, instrument, dayint, moint, yearint, low=0,
                  high=800, ylo=1e-3, yhi=10):
    """
    Processing of Fermi data in the 25-300 keV band (though applicable to
    others), from a .sav file generated using the Fermi OSPEX database.

    Parameters
    ----------
    day : str
        Day of flare, numerical string format (DD).
    month : str
        Month of flare, numerical string format (MM).
    year : str
        Year of flare, numerical string format (YYYY).
    instrument : str
        Corresponding Fermi instrument (n5, typically).
    dayint : int
        Day of flare, integer format.
    moint : int
        Month of flare, integer format.
    yearint : int
        Year of flare, integer format.
    low : int, optional
        Lower limit of range to search for flare in curve. The default is 0.
    high : int, optional
        Upper limit of range to search for flare in curve. The default is 800.

    Returns
    -------
    raw_hxr_sum : arr
        Raw 25-300 keV energy from Fermi.
    cspec_hxr_sum : arr
        Background-subtracted 25-300 keV energy from Fermi.
    fermitimes : arr
        Timestamps from Fermi data file.

    """

    directory = '/Users/coletamburri/Desktop/Impulsiveness_Paper/Imp_Study_Transfer/'\
        'Fermi_Events_sav/'

    filename_cspec = directory + 'fermi_' + instrument + '_cspec_bkgd_' + day \
        + month + year + '.sav'

    cspec_dat = readsav(filename_cspec, python_dict='True')

    bksub_cspec = cspec_dat['lc_bksub'][0][0]
    raw_cspec = cspec_dat['lc_raw'][0][0]
    times = cspec_dat['time']
    energies = cspec_dat['ct_energy']

    hxrinds = np.where(cspec_dat['ct_energy'] < 300.) and \
        np.where(cspec_dat['ct_energy'] > 25.)

    cspec_hxr = bksub_cspec[:, hxrinds]
    raw_hxr = raw_cspec[:, hxrinds]

    cspec_hxr_sum = np.zeros(len(times))
    raw_hxr_sum = np.zeros(len(times))
    
    cspec_hxr_sum = np.sum(cspec_hxr, axis=2)
    raw_hxr_sum = np.sum(raw_hxr, axis=2)

    a = datetime.datetime(1970, 1, 1, 0, 0, 0)
    b = datetime.datetime(1979, 1, 1, 0, 0, 0)

    err1 = (b-a).total_seconds()

    timesadj1 = times + err1

    curr = datetime.datetime.fromtimestamp(min(timesadj1))
    corr = datetime.datetime(yearint, moint, dayint, 0, 0, 0)

    err2 = (corr-curr).seconds
    totsec = (b-a).total_seconds() + err2

    timesadj = times + totsec

    timepkg.ctime(min(timesadj))
    strtimes = []

    for i in timesadj:
        strtimes.append(datetime.datetime.fromtimestamp(i))

    flag = 0

    for i in range(low, high):
        if cspec_hxr_sum[i] > 0.1:
            flag += 1
        if cspec_hxr_sum[i] < 0.1:
            flag = 0

    fermitimes = strtimes

    return raw_hxr_sum, cspec_hxr_sum, fermitimes


def E_field_det(conv_f, distpos, distneg, timelab, hmi_dat, pos_rem, neg_rem,
                flnum, dt1600, times, startind=25):
    """


    Parameters
    ----------
    conv_f : float
        Conversion factor between pixels and Mm.
    distpos : arr
        Positive ribbon separation values.
    distneg : arr
        Negative ribbon separation values.
    timelab : arr
        Times corresponding 1600 A images.
    hmi_dat : arr
        SDO/HMI magnetogram image.
    pos_rem : arr
        Vetted ribbon images, positive polarity.
    neg_rem : arr
        Vetted ribbon images, negative polarity
    flnum : int
        RibbonDB event number
    dt1600 : arr
        Datetimes corresponding to 1600 A images.
    times : arr
        Times corresponding to 1600 A images.
    startind : int, optional
        End of transient period. The default is 25.

    Returns
    -------
    E_pos : arr
        Positive electric field time series.
    E_neg : arr
        Negative electric field time series.
    E_rat : arr
        Ratio of positive to negative electric field.
    time_E : arr
        Times corresponding to electric field stamps.

    """


    
    hmi_pos = np.zeros(np.shape(hmi_dat))
    hmi_neg = np.zeros(np.shape(hmi_dat))
    sum_pos = 0
    count_pos = 0
    sum_neg = 0
    count_neg = 0

    s = str(dt1600[0])

    pos_Mm_dist = conv_f*distpos[startind:]
    neg_Mm_dist = conv_f*distneg[startind:]
    time_E = timelab[startind:]*60  # conversion to seconds

    time_res = time_E[1] - time_E[0]

    slab = s[startind:]

    # velocities, in Mm/s

    pos_v_Mm = scipy.signal.medfilt(np.diff(pos_Mm_dist)/time_res,
                                    kernel_size=5)
    neg_v_Mm = scipy.signal.medfilt(np.diff(neg_Mm_dist)/time_res,
                                    kernel_size=5)

    # convert to m/s for SI units

    pos_v = pos_v_Mm * 1e6
    neg_v = neg_v_Mm * 1e6

    # setting velocities equal to zero because these are not real
    for i in range(len(pos_v)):
        if pos_v[i] < 0:
            pos_v[i] = 0
        if neg_v[i] < 0:
            neg_v[i] = 0

    for i in range(len(hmi_dat)):
        for j in range(len(hmi_dat[0])):
            if hmi_dat[i, j] > 0.0 and pos_rem[i, j] == 1:
                hmi_pos[i, j] = hmi_dat[i, j]

            if hmi_dat[i, j] < 0.0 and neg_rem[i, j] == -1:
                hmi_neg[i, j] = hmi_dat[i, j]

    for i in range(len(hmi_dat)):
        for j in range(len(hmi_dat[0])):
            if hmi_pos[i, j] > 0:
                sum_pos += hmi_pos[i, j]
                count_pos += 1
            if hmi_neg[i, j] < 0:
                sum_neg += hmi_neg[i, j]
                count_neg += 1

    avg_B_pos = sum_pos/count_pos  # in G
    avg_B_neg = sum_neg/count_neg  # in G

    E_pos = avg_B_pos*pos_v*1e-6  # in V/cm
    E_neg = -avg_B_neg*neg_v*1e-6  # in V/cm

    E_rat = E_pos/E_neg

    fig, ax = plt.subplots(figsize=(20, 7))

    ax.plot(time_E[0:-1]/60, E_pos, '--ro', label='$E_{pos}$')
    ax.plot(time_E[0:-1]/60, E_neg, '--bo', label='$E_{neg}$')
    ax.set_xlabel('Time [min since '+s[5:-7]+']', font='Times New Roman',
                  fontsize=18)
    ax.set_xlim([0, time_E[-1]/60])
    ax.grid()
    ax.set_ylabel('Electric Field [V/cm]', font='Times New Roman',
                  fontsize=18)
    ax.set_title('Reconnection Electric Field Strength',
                 font='Times New Roman', fontsize=45)
    ax.legend()

    fig.savefig(str(flnum) + '_E_field.png')

    return E_pos, E_neg, E_rat, time_E


def shear_to_angle(times, flnum, dt1600, left_gfr, right_gfr, flag=0):
    """
    Determination of shear angle time series.

    Parameters
    ----------
    times : arr
        Times corresponding to 1600 angstrom images.
    flnum : int
        Flare number
    dt1600 : arr
        Datetime values corresponding to 1600 angstrom images.
    left_gfr : arr
        Time series of guide field ratio, left of AR.
    right_gfr : arr
        Time series of guide field ratio, right of AR
    flag : int, optional
        0 if traditional shear method, 1 if antiparallel.

    Returns
    -------
    shear_ang : arr
        Array of angles corresponding to determined shear.

    """
    ang_prep = (left_gfr+right_gfr)/2

    shear_ang = np.arctan(1/ang_prep)*180./np.pi

    timelab = range(0, 24*len(times), 24)
    s = str(dt1600[0])
    fig, ax = plt.subplots(figsize=(13, 7))
    if flag == 0:
        ax.plot(timelab, shear_ang, c='green', marker='o',
                label='Shear Angle, Left')
    elif flag == 1:
        ax.plot(timelab, right_gfr, c='black', marker='o', label='GFR proxy')
        ax.set_xlabel('Time [s since '+s[5:-7]+']', font='Times New Roman',
                      fontsize=18)
        ax.set_ylabel('Shear Angle', font='Times New Roman', fontsize=18)

    ax.set_title('Shear Angle', font='Times New Roman', fontsize=20)
    ax.grid()
    ax.legend(fontsize=15)
    fig.savefig(str(flnum) + '_shear_angle.png')

    return shear_ang


def quartermaxtime(transtime, right_gfr, left_gfr, timelab, find_nearest_idx,
                   flag=0):
    """
    Determination of time to quarter-max of shear relative to minimum of shear.

    Parameters
    ----------
    transtime : int
        Number of points in transient.
    right_gfr : arr
        Guide field ratio, right of AR.
    left_gfr : arr
        Guide field ratio, left of AR.
    timelab : arr
        dt1600 time labels.
    find_nearest_idx : function
        Function to find index of array closest to a certain value.
    flag : int, optional
        0 if traditional shear determination, 1 if antiparallel.  Default is 0.

    Returns
    -------
    quartermax_time : float
        Minutes since flare start until the shear has reached 25% of max.

    """
    if flag == 0:
        gfr = (right_gfr+left_gfr)/2
    if flag == 1:
        gfr = right_gfr

    mingfr = np.nanmin(gfr[transtime:])
    maxgfr = np.nanmax(gfr[transtime:])

    minmaxavg = mingfr+.25*(maxgfr-mingfr)

    idx = find_nearest_idx(gfr[transtime:], minmaxavg)

    quartermax_time = timelab[transtime+idx] - timelab[transtime]

    return quartermax_time

def color_muted():
    muted =['#332288',  '#88CCEE', '#44AA99', '#117733','#999933', '#DDCC77',
            '#CC6677',  '#882245','#AA4499', '#DDDDDD']
    # 0=indigo
    # 1=cyan
    # 2=teal
    # 3= green
    # 4=olive
    # 5= sand
    #6 = rose
    #7= wine
    # 8 = purple
    return muted

def color_vibrant():
    vibrant = ['#0077BB',  '#33BBEE', '#009988','#EE7733', '#CC3311', 
               '#EE3377', '#BBBBBB'] 
    # CAT - changed grey to #5A5A5A to make darker; originally #BBBBBB
    return vibrant

def color_medc(printc='no'):
    contrast = [(245,245,245), (238,204,102), (238,153,170), (102,153,204), 
                (153,119,0), (153,68,85), (0,68,136), (0,0,0)]
    
    # 0 is white, 1 is light yellow, 2 is light red, 3 is light blue,
    # 4 is dark yellow, 5 is dark red, 6 is dark blue, 7 is black
    for i in range(len(contrast)):    
        r, g, b = contrast[i]    
        contrast[i] = (r / 245., g / 245., b / 245.)

    return contrast



def plt_fourpanel(times, right_gfr, left_gfr, flnum, dt1600, time304,
                  filter_304, lens_pos_Mm, lens_neg_Mm, distpos_Mm, distneg_Mm,
                  dt304, timelab, conv_f,
                  elonperiod_start_pos, elonperiod_end_pos,
                  elonperiod_start_neg, elonperiod_end_neg,
                  sepperiod_start_pos, sepperiod_end_pos,
                  sepperiod_start_neg, sepperiod_end_neg, exp_ind,
                  s304, e304, pos1600, neg1600, dn1600, indstrt_elon,
                  indstrt_sep, fermitimes, raw_hxr_sum, cspec_hxr_sum,
                  gfr_trans, E_pos, E_neg, time_E, day, mo, year, xcl, xclnum,
                  imp, muted, vibrant, medc, sxr_fn, level,low_hxr=0, high_hxr=800, period_flag=0, fermioff=0, flag=0,
                  tick_space=0,fllab='Low1',locationleg='upper left',ylimflag='off',pilparlim=[-0.5,35],pilperplim=[-0.5,6.5],
                  gfrlim=[-0.5,12.5],Ereclim=[-0.5,8.5],sxroff=0,ax1legloc='upper left'):
    """
    Four-panel plot to compare HXR/1600 Angstrom/304 Angstrom (panel 1),
    ribbon separation (panel 2), ribbon elongation (panel 3), guide field ratio
    proxy (panel 4, a measurement of shear).  the peak times for the Fermi HXR
    25-300 keV band, 304 Angstrom light curve, and 1600 Angstrom light curves
    are shown in each panel, and timestamps accurately lined up to show periods
    of separation/elongation/shear.  There is an option to plot, also, the
    periods of significant separation and elongation, though this may make
    panels 2 and 3 too busy.

    Parameters
    ----------
    times : arr
        Array of times for the flare, from AIA datafile.
    right_gfr : arr
        Guide field ratio, determined from the right edge of the ribbons.
    left_gfr : arr
        Guide field ratio, determined from the left edge of the ribbons.
    flnum : int
        Flare number
    dt1600 : arr
        Datetime values from 1600 Angstrom data file.
    time304 : arr
        Datenum values from MATLAB 304 Angstrom data file
    filter_304 : arr
        Filtered 304 Angstrom data.
    lens_pos_Mm : list
        Perpendicular extent of positive ribbon for each time step, in Mm.
    lens_neg_Mm : list
        Parallel extent of positive ribbon for each time step in Mm.
    distpos_Mm : list
        Perpendicular extent of negative ribbon for each time step, in Mm.
    distneg_Mm : list
        Parallel extent of positive ribbon for each time step in Mm.
    dt304 : list
        Datetime values for 304 Angstrom data.
    timelab : list
        Preparation of time labels for future plotting of light curves.
    conv_f : float
        Conversion factor between pixels and Mm, approximate assuming no
        curvature
    elonperiod_start_pos : list
        Determined start times for elongation in positive ribbon.
    elonperiod_end_pos : list
        Determined end times for elongation in positive ribbon.
    elonperiod_start_neg : list
        Determined start times for elongation in negative ribbon.
    elonperiod_end_neg : list
        Determined end times for elongation in negative ribbon.
    sepperiod_start_pos : list
        Start times for periods of extended separation, positive ribbon.
    sepperiod_end_pos : list
        End times for periods of extended separation, positive ribbon.
    sepperiod_start_neg : list
        Start times for periods of extended separation, negative ribbon.
    sepperiod_end_neg : list
        End times for periods of extended separation, negative ribbon.
    exp_ind : int
        The index where to stop the exponential fitting.
    s304 : int
        Nearest index in AIA data to start time of the flare from EUV 304 light
        curve.
    e304 : int
        Nearest index in AIA data to end time of the flare from EUV 304 light
        curve.
    pos1600 : list
        Summed pixel numbers in positive ribbon for each time step.
    neg1600 : list
        Summed pixel numbers in negative ribbon for each time step.
    dn1600 : list
        Datenum values for 1600 Angstrom data.
    indstrt_elon : int
        Index at which to begin plotting elongation.
    indstrt_sep : int
        Index at which to begin plotting separation.
    fermitimes : arr
        Timestamps from Fermi data file.
    raw_hxr_sum : arr
        Raw 25-300 keV energy from Fermi.
    cspec_hxr_sum : arr
        Background-subtracted 25-300 keV energy from Fermi.
    gfr_trans : arr
        Index for the end of the Guide Field Ratio transient period.
    low_hxr : int, optional
        Low index for flare search in HXR data. The default is 0.
    high_hxr : int, optional
        High index for flare search in HXR data. The default is 800,
        applied to Oct-13-2013 flare.
    period_flag : TYPE, optional
        DESCRIPTION. The default is 0.
    flag : int, optional
        0 (if traditional method of shear) or 1 (if antiparallel method of
        shear). The default is 0.
    tick_space : int, optional
        Flag, which, when nonzero and set to the number of ticks wanted on axis
        labels, triggers that change; 0 otherwise, defaulting to matplotlib.

    Returns
    -------
    None.

    """

    plt.rcParams['text.usetex']=True
    plt.rcParams['font.family']='sans-serif'
    plt.rcParams['font.sans-serif'] = ['Tahoma']
    plt.rc('text.latex', preamble=r'\usepackage[cm]{sfmath}')

    #rc = {"font.family" : "sans-serif", 
    #  "mathtext.fontset" : "cm"}
    #plt.rcParams.update(rc)
    #plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

    
    dt1600_1 = []
    dt304_1 = []

    for i in range(len(dt1600)):
        dt1600_1.append(dt1600[i].strftime("%H:%M"))
    for i in range(len(dt304)):
        dt304_1.append(dt304[i].strftime("%H:%M"))

    min304 = min(filter_304[s304: e304])
    max304 = max(filter_304[s304: e304])
    minpos1600 = min(pos1600)
    maxpos1600 = max(pos1600)
    minneg1600 = min(neg1600)
    maxneg1600 = max(neg1600)

    # Normalize for light curve comparison

    norm304 = (filter_304 - min304) / (max304 - min304)
    normpos1600 = (pos1600 - minpos1600) / (maxpos1600 - minpos1600)
    normneg1600 = (neg1600 - minneg1600) / (maxneg1600 - minneg1600)
    
    
    euvdf = pd.DataFrame({'euv':norm304},columns=['euv'])
    euvdf_fix = euvdf.fillna(method='ffill').fillna(method='bfill')

    if flag == 0:
        GFR = np.mean([right_gfr, left_gfr], axis=0)
    elif flag == 1:
        # if the alterative version of GFR, take only right_gfr (the input
        # should just be gfr)
        GFR = right_gfr

    

    max304_0 = np.nanargmax(filter_304)
    max304t = dt304[max304_0]
    max304 = find_nearest_ind(dt1600, max304t)

    max1600pos = np.argmax(normpos1600)
    max1600neg = np.argmax(normneg1600)
    fermi = cspec_hxr_sum[low_hxr:high_hxr,0]/275.0
    minfermi = min(fermi)
    maxfermi = max(fermi)
    
    normfermi = (fermi - minfermi)/(maxfermi-minfermi)
    
    #fermi = np.log10(cspec_hxr_sum[low_hxr:high_hxr,0]/275.0)#divide by wavelength band for consistency with OSPEX; file gives cts/s/cm^2, this will divide by wavelength range after summing over all
    #fermi = cspec_hxr_sum[low_hxr:high_hxr,0]/275.0
    fermitimessel = fermitimes[low_hxr:high_hxr]
    fermidf = pd.DataFrame({'hxr':fermi},columns=['hxr'])
    fermidf_fix = fermidf.fillna(method='ffill').fillna(method='bfill')

    hxrmax0 = np.argmax(fermidf_fix.hxr)
    
    hxrmaxt = fermitimessel[hxrmax0]
    
    hxrmax = find_nearest_ind(dt1600, hxrmaxt)
        
    # 0 is white, 1 is light yellow, 2 is light red, 3 is light blue,
    # 4 is dark yellow, 5 is dark red, 6 is dark blue, 7 is black
    # processing sxr
    
    
    
    ds = nc.Dataset(sxr_fn)
    sxr18flux = ds['b_flux']
    
    times_thisday = []
    for i in range(len(ds['time'])):
        times_thisday.append(nc.num2date(int(ds['time'][i].data),units=ds['time'].units,only_use_cftime_datetimes=False,only_use_python_datetimes=True))

    sxrlow = find_nearest_ind(times_thisday,dt1600[0])
    sxrhigh = find_nearest_ind(times_thisday,dt1600[-1])
    minsxr = min(sxr18flux[sxrlow:sxrhigh])
    maxsxr = max(sxr18flux[sxrlow:sxrhigh])
    
    normsxr = (sxr18flux - minsxr) / (maxsxr - minsxr)
    
    #fig, ((ax1, ax3),(ax2, ax4)) = plt.subplots(2,2, figsize=(40, 20))
    fig = plt.figure(figsize=(13, 35))
    # set height ratios for subplots
    gs = gridspec.GridSpec(3,1) 
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0],sharex=ax1)
    ax3 = plt.subplot(gs[2,0],sharex=ax1)
    #ax4 = plt.subplot(gs[3,0],sharex=ax1)
    
    lns1 = ax1.plot(dt1600, normpos1600,markersize=5,linewidth=6, color=vibrant[4],
                    label=r'1600 \AA, +')
    lns2 = ax1.plot(dt1600, normneg1600, color='#81C4E7',linewidth=6,
                    label=r'1600 \AA, -')

    lns3 = ax1.plot(dt304[s304:e304], scipy.signal.medfilt(euvdf_fix.euv[s304:e304],
                                                           kernel_size=5),color=medc[1], 
                    linewidth=6,
                    label=r'304 \AA')

    
    if fermioff  == 0 and sxroff == 0:

        # minfer = min(scipy.signal.medfilt(
        #     fermidf_fix.hxr[0:500], kernel_size=9))
        # maxfer = max(scipy.signal.medfilt(
        #     fermidf_fix.hxr[0:500], kernel_size=9))
        # ax1_0 = ax1.twinx()
        # ax1_0.set_ylim(minfer-0.25,maxfer+0.25)
        # lns4 = ax1_0.plot(fermitimes[low_hxr:high_hxr],
        #                   scipy.signal.medfilt(
        #                       fermidf_fix.hxr, kernel_size=21),
        #                   label='HXR', color='#4EB265',linewidth=6)
        # ax1_0.set_ylabel(
        #     r'HXR Flux [$log(cts/s/cm^{2}/keV)$]',
        #     fontsize=45)
        lns6 = ax1.plot(times_thisday,normsxr,linewidth=6,label=r'GOES SXR')
        lns7 = ax1.plot(fermitimes[low_hxr:high_hxr], scipy.signal.medfilt(normfermi,kernel_size=5),color='#4EB265', 
                        linewidth=6,
                        label='Fermi HXR')
        lns = lns1+lns2+lns3+lns6+lns7
        #ax1.axvline(dt1600[hxrmax],color='#4EB265',linestyle='dashed',linewidth=4)
    elif fermioff == 0 and sxroff == 1:
        lns7 = ax1.plot(fermitimes[low_hxr:high_hxr], scipy.signal.medfilt(normfermi,kernel_size=5),color='#4EB265', 
                        linewidth=6,
                        label='Fermi HXR')
        lns = lns1+lns2+lns3+lns7
    elif fermioff == 1 and sxroff == 0:
        lns6 = ax1.plot(times_thisday,normsxr,linewidth=6,label=r'GOES SXR')
        lns = lns1+lns2+lns3+lns6
   
    ax1.axes.xaxis.set_ticklabels([])

    ax1.grid()

    labs = [k.get_label() for k in lns]
    font = font_manager.FontProperties(style='normal', size=35,family='Tahoma')
    ax1.legend(lns, labs, prop=font,loc=ax1legloc,framealpha=0.95)
    ax1.set_ylabel('Normalized Light Curves',
                   fontsize=45)
    ax1.set_ylim([-0.1,1.1])
   
    #ax1.set_title('Flare Energy Dynamics',
    #              fontsize=50)
    ax1.set_xlim([dt1600[0], dt1600[-1]])
    #ax1.axvline(dt1600[max304], color=medc[1],linestyle='dotted',linewidth=4)
    #ax1.axvline(dt1600[max1600pos], color=vibrant[4],linestyle='dashdot',linewidth=4)
    #ax1.axvline(dt1600[max1600neg], color='#81C4E7', linestyle='dashdot',linewidth=4)
    lns1 = ax2.plot(dt1600[gfr_trans:-1], GFR[gfr_trans:-1], c=muted[2],
                    linewidth=6,
                    label='GFR')

    
    ax2.set_ylabel('Guide Field Ratio (GFR)', fontsize=45)
    #ax2.set_title(r'(b) Magnetic Shear and $|E_{rec}|$',
    #              fontsize=50)
    #ax2.set_ylim([-0.2, np.nanmax(GFR[gfr_trans:])+2])
    ax2.set_ylim(gfrlim)
    ax2.grid()
    

    ax2.set_xlim([dt1600[0], dt1600[-1]])
    #ax2.axvline(dt1600[hxrmax],color='#4EB265',linestyle='dashed',linewidth=4)
    #ax2.axvline(dt1600[max304], color=medc[1],linestyle='dotted',linewidth=4)
    #ax2.axvline(dt1600[max1600pos], color=vibrant[4],linestyle='dashdot',linewidth=4)
    #ax2.axvline(dt1600[max1600neg], color='#81C4E7', linestyle='dashdot',linewidth=4)
    font = font_manager.FontProperties(style='normal', size=35)
    #ax2.legend(prop=font, fontsize=20)
    ax2_0 = ax2.twinx()
    lns2 = ax2_0.plot(dt1600[gfr_trans:], E_pos, color=vibrant[4],linewidth=6, label='$|E_{rec, +}|$')
    lns3 = ax2_0.plot(dt1600[gfr_trans:], E_neg, color='#81C4E7',linewidth=6, label='$|E_{rec, -}|$')

    ax2_0.set_xlim([dt1600[0], dt1600[-1]])

    ax2_0.set_ylabel(r'$|E_{rec}|$ $[V\; cm^{-1}]$',
                     fontsize=45)
    
    ax2_0.set_ylim(Ereclim)


    lns = lns1+lns2+lns3
    labs = [k.get_label() for k in lns]
    font = font_manager.FontProperties(style='normal', size=35,family='Tahoma')
    ax2.legend(lns, labs, prop=font, fontsize=35)
    ax2.axes.xaxis.set_ticklabels([])
    ax2_0.axes.xaxis.set_ticklabels([])
    # regions of separation/elongation make a littlfe busy?

    if period_flag == 1:
        for i, j in zip(sepperiod_start_pos, sepperiod_end_pos):
            #ax2.axvline(dt1600[i], c='green')
            #ax2.axvline(dt1600[j], c='red')
            ax2.axvspan(dt1600[i], dt1600[j], alpha=0.5, color='pink')
        for k, l in zip(elonperiod_start_neg, elonperiod_end_neg):
            #ax2.axvline(dt1600[k], c='green')
            #ax2.axvline(dt1600[l], c='red')
            ax2.axvspan(dt1600[k], dt1600[l], alpha=0.5, color='cyan')

    lns1 = ax3.plot(dt1600[indstrt_elon:-1], lens_pos_Mm[indstrt_elon:-1],
             c=vibrant[4], linewidth=6, label = r'$d_{\parallel, +}$')

    lns2 = ax3.plot(dt1600[indstrt_elon:-1], lens_neg_Mm[indstrt_elon:-1],
             c='#81C4E7',linewidth=6, label = r'$d_{\parallel, -}$')

    ax3.grid()
    ax3.set_ylabel(
        r'PIL-Parallel Distance ($d_{\parallel}$) [$Mm$]', fontsize=45)
    #ax3.set_title('PIL-Relative Ribbon Motion',
    #              fontsize=50)

    

    ax3.set_xlim([dt1600[0], dt1600[-1]])
    #ax3.yaxis.set_label_position("right")
    ax4 = ax3.twinx()
    #ax3.set_ylim([0,33])
    ax3.set_ylim(pilparlim)

    #ax3.axvline(dt1600[max304], color=medc[1],linestyle='dotted',linewidth=4)
    #ax3.axvline(dt1600[max1600pos], color=vibrant[4],linestyle='dashdot',linewidth=4)
    #ax3.axvline(dt1600[max1600neg], color='#81C4E7', linestyle='dashdot',linewidth=4)
    
    ax3.axes.xaxis.set_ticklabels([])
    lns3 = ax4.plot(dt1600[indstrt_sep:-1], distpos_Mm[indstrt_sep:-1], c=vibrant[4],
             markersize=6, linestyle='dashed',linewidth=6, label = r'$d_{\perp, +}$')

    lns4 = ax4.plot(dt1600[indstrt_sep:-1], distneg_Mm[indstrt_sep:-1],
             c='#81C4E7', linestyle='dashed',linewidth=6, label = r'$d_{\perp, -}$')

    ax4.set_ylabel(
        r'PIL-Perpendicular Distance ($d_{\perp}$) [$Mm$]', fontsize=45,fontname="Tahoma")
    #ax4.set_title('(d) Distance Perpendicular to PIL',
    #              fontsize=50)
    #ax4.yaxis.set_label_position("right")

    ax4.set_xlim([dt1600[0], dt1600[-1]])
    


    #ax4.set_ylim([0,2.25])
    ax4.set_ylim(pilperplim)
    #ax4.axvline(dt1600[hxrmax],color='#4EB265',linestyle='dashed',linewidth=4)
    #ax4.axvline(dt1600[max304], color=medc[1],linestyle='dotted',linewidth=4)
    #ax4.axvline(dt1600[max1600pos], color=vibrant[4],linestyle='dashdot',linewidth=4)
    #ax4.axvline(dt1600[max1600neg], color='#81C4E7', linestyle='dashdot',linewidth=4)
    #ax4.grid()
    ax3.set_xlabel('Time, '+year+' '+mo+' '+day+r' [$UT$]',
                   fontsize=45)
    lns = lns1+lns2+lns3+lns4
    labs = [k.get_label() for k in lns]
    font = font_manager.FontProperties(style='normal', size=45,family='Tahoma')
    ax3.legend(lns, labs, prop=font, fontsize=35,loc=locationleg)
    #lns = lns1+lns2+lns3+lns4
    #labs = [k.get_label() for k in lns]
    #font = font_manager.FontProperties(style='normal', size=45,family='Tahoma')
    #ax4.legend(lns, labs, prop=font, fontsize=45)
    ax2.set_xlabel('Time, '+year+' '+mo+' '+day+r' [$UT$]',
                   fontsize=45)

    if tick_space > 0:
        ax1.set_xticks(dt1600[2::tick_space])
        #ax2.set_xticks(dt1600[2::tick_space])
        ax3.set_xticks(dt1600[2::tick_space])
        ax4.set_xticks(dt1600[2::tick_space])
        #if fermioff == 0:
        #    ax1_0.set_xticks(dt1600[2::tick_space])
        #ax2_0.set_xticks(dt1600[2::tick_space])
        ax1.set_xticklabels(dt1600_1[2::tick_space])
        #ax2.set_xticklabels(dt1600_1[2::tick_space])
        ax3.set_xticklabels(dt1600_1[2::tick_space])
        ax4.set_xticklabels(dt1600_1[2::tick_space])
        #if fermioff == 0:
        #    ax1_0.set_xticklabels(dt1600_1[2::tick_space])
        #ax2_0.set_xticklabels(dt1600_1[2::tick_space])
       
    #if fermioff  == 0:
    #    ax1_0.set_yticks([-2,-3,-4])
        
    #ax3.yaxis.tick_right()
    #ax4.yaxis.tick_right()


    labsize=45
    ax1.xaxis.set_tick_params(labelsize=labsize)
    ax1.yaxis.set_tick_params(labelsize=labsize)
    ax2.xaxis.set_tick_params(labelsize=labsize)
    ax2.yaxis.set_tick_params(labelsize=labsize)
    ax3.xaxis.set_tick_params(labelsize=labsize)
    ax3.yaxis.set_tick_params(labelsize=labsize)
    ax4.xaxis.set_tick_params(labelsize=labsize)
    ax4.yaxis.set_tick_params(labelsize=labsize)
    #if fermioff == 0:
    #    ax1_0.xaxis.set_tick_params(labelsize=labsize)
    #    ax1_0.yaxis.set_tick_params(labelsize=labsize)
    ax2_0.xaxis.set_tick_params(labelsize=labsize)
    ax2_0.yaxis.set_tick_params(labelsize=labsize)
    fig.align_ylabels([ax1,ax2])
    

    #if fermioff == 0:
    #    fig.align_ylabels([ax1_0,ax2_0])
    fig.align_ylabels([ax3,ax4])

    
    ax1.set_title('Event '+fllab,fontsize=55)
    #ax1.set_title('Event '+fllab+'\n'+'(GOES '+xcl+str(xclnum)+', '+
    #             r'$i =$ '+str(imp)+r' $ln[min^{-1}]$)',fontsize=55)
    
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2_0.get_xticklabels(), visible=False)
    #plt.setp(ax3.get_xticklabels(), visible=False)
    plt.subplots_adjust(hspace=.0)

    if period_flag == 1:
        for i, j in zip(elonperiod_start_pos, elonperiod_end_pos):
            print(0)
            #ax3.axvline(dt1600[i], c='green')
            #ax3.axvline(dt1600[j], c='red')
            #ax3.axvspan(dt1600[i], dt1600[j], alpha=0.5, color='pink')
        for k, l in zip(elonperiod_start_neg, elonperiod_end_neg):
            print(0)#ax3.axvline(dt1600[k], c='green')
            #ax3.axvline(dt1600[l], c='red')
            #ax3.axvspan(dt1600[k], dt1600[l], alpha=0.5, color='cyan')

   
    fig.savefig(str(flnum) + '_summary.png')

    return None
