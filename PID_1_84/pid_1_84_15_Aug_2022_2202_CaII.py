#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 15:34:06 2026

@author: coletamburri
"""

# shift of wavelength range by inspection
end=0

# package initialize
import dkistpkg_ct as DKISTanalysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update(matplotlib.rcParamsDefault)

# color scheme for plotting
muted = DKISTanalysis.color_muted2()

# path and file ID for ViSP data
path = '/Volumes/ViSP_External/pid_1_84/'
path2 = '/Volumes/ViSP_External/pid_1_84/'
folder1 = 'XVOUZY'
folder2 = 'CZGJML'  #need to add data for QS for 11 august...

# list of files in directory for DKIST/ViSP
dir_list2 = DKISTanalysis.pathdef(path,folder1) #flaretime
dir_list3 = DKISTanalysis.pathdef(path2,folder2) #qs

caII_low = 570
caII_high = 775
hep_low = 775
hep_high = 918