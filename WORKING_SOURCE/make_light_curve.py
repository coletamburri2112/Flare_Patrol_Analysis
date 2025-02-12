#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 11:29:43 2025

@author: coletamburri
"""
import numpy as np
import matplotlib.pyplot as plt
import tol_colors as tc

import dkist
import matplotlib.pyplot as plt

from sunpy.net import Fido, attrs as a
import dkist.net
from datetime import date, datetime, timedelta, time
import matplotlib.dates as mdates

import urllib
import numpy as np
import os
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import dates
import fl_funcs as ff
def perdelta(start, end, delta):
    curr = start
    while curr <= end:
        yield curr
        curr += delta

dtfmt = '%Y-%m-%d %H:%M:%S'

sttime = datetime(2024,8,8,20,12,32,333333)
endtime = datetime(2024,8,8,21,5,7,0)


stack=[]
for result in perdelta(sttime , endtime, timedelta(seconds=2.666666)):
    stack.append(result)
    
timeshhmmss = []

# timesdt = []

# for i in stack:
#     timesdt.append(datetime.strptime(i, '%Y-%m-%dT%H:%M:%S.%f'))
    
# timesdthr = []

# for i in timesdt:
#     timesdthr.append(datetime.strftime(i, 
#                                  "%H:%M:%S"))

# times in correct format for plotting
for i in range(len(stack)):
    timeshhmmss.append(stack[i][-15:-7])
    
# new_x = "/Users/coletamburri/Desktop/8_August_2024.nc"
# x = 'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/xrsf-l2-flx1s_science/2022/12/sci_xrsf-l2-flx1s_g16_d20221228_v2-2-0.nc'
# urllib.request.urlretrieve(x, new_x)
fn = "/Users/coletamburri/Desktop/DKIST_Data_Tools_misc/GOES_NETCDF_Flare_Files/8_August_2024.nc"
ds = nc.Dataset(fn)

def datenum_to_datetime(datenum):
    """
    Convert Matlab datenum into Python datetime.
    :param datenum: Date in datenum format
    :return:        Datetime object corresponding to datenum.
    """
    days = datenum % 1
    hours = days % 1 * 24
    minutes = hours % 1 * 60
    seconds = minutes % 1 * 60
    return datetime.fromordinal(int(datenum)) \
           + timedelta(days=int(days)) \
           + timedelta(hours=int(hours)) \
           + timedelta(minutes=int(minutes)) \
           + timedelta(seconds=round(seconds)) \
           - timedelta(days=366)

#x = [datetime.strptime(d,'%Y-%m-%d %H:%M:%S.%f').time() for d in stack]

muted = tc.tol_cset('muted')

data=np.load('/Users/coletamburri/Desktop/August_2024_DKIST_Flares/VBI_X_class/X_class_lc_8Aug2024_save.npz')

fig,ax=plt.subplots(dpi=200)
ax.axvline(datetime(2024,8,8,20,12,32,333333),color='grey')
ax.axvline(datetime(2024,8,8,21,5,7,0),color='grey')
ax.axvspan(datetime(2024,8,8,20,12,32,333333),
           datetime(2024,8,8,21,5,7,0), alpha=0.2,label = 'DKIST/VBI Obs.',color='grey')

#ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S.%f'))
#ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=20))
#ax.xaxis.set_major_locator(x.DayLocator())

#now plot the GOES curve...
t = np.arange(datetime(2024,8,8,0,0,0),
              datetime(2024,8,9,0,0,0), 
              timedelta(seconds=1)).astype(datetime)
ax1 = ax.twinx()
lns3=ax1.plot(t,np.array(ds['xrsb1_flux']),linewidth=3,color='k',label='GOES Flux')
ax.set_xlim(datetime(2024,8,8,18,30),datetime(2024,8,8,22))
ax1.set_xlim(datetime(2024,8,8,18,30),datetime(2024,8,8,22))
ax1.set_ylabel('GOES Flux [$W\;m^{-2}$]')
plt.yscale("log")

#DKIST
lns1=ax.plot(stack,data['arr_0']/1e4,c=muted[1],label='Full VBI FoV');
#ax1=ax.twinx();
#lns2 = ax1.plot(data['arr_1'],c=muted[0],label='Loop-top source');
lns2=ax.plot(stack,data['arr_1']/1e4,c=muted[0],label='Loop-top source');
ax.set_ylabel('Avg. Intensity over FoV [$10^4\;DN\;pix^{-2}$]')
#ax.autofmt_xdate()
#ax.set_xticks(timeshhmmss[0:-1:190],timeshhmmss[0:-1:190],rotation=45)
ax.set_xlabel('Time [UT]')
ax.set_ylim([1.3,2.5])
ax1.set_ylim([3e-6,1e-3])


fmtr = dates.DateFormatter("%H:%M")
# need a handle to the current axes to manipulate it
ax = plt.gca()
# set this formatter to the axis
ax.xaxis.set_major_formatter(fmtr)

lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0,fontsize=10)

   

#now plot SolO/PSP



fig.subplots_adjust(bottom=0.2)





