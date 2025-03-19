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
# for i in range(len(stack)):
#     timeshhmmss.append(stack[i][-15:-7])
    
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
lns3=ax1.plot(t,np.array(ds['xrsb1_flux']),linewidth=3,color='k',label='GOES Flux',zorder=0)
ax.set_xlim(datetime(2024,8,8,18,30),datetime(2024,8,8,22))
ax1.set_xlim(datetime(2024,8,8,18,30),datetime(2024,8,8,22))
ax1.set_ylabel('GOES Flux [$W\;m^{-2}$]')
plt.yscale("log")

#DKIST
lns1=ax.plot(stack,data['arr_0']/1e4,c=muted[1],label='Full VBI FoV');
#ax1=ax.twinx();
#lns2 = ax1.plot(data['arr_1'],c=muted[0],label='Loop-top source');
lns2=ax.plot(stack,data['arr_1']/1e4,c=muted[0],label='Loop-top source',zorder=5);
ax.set_ylabel('Avg. Intensity over FoV [$10^4\;DN\;pix^{-2}$]')
#ax.autofmt_xdate()
#ax.set_xticks(timeshhmmss[0:-1:190],timeshhmmss[0:-1:190],rotation=45)
ax.set_xlabel('Time [UT]')
ax.set_ylim([1.7,2.5])
ax1.set_ylim([3e-6,1e-3])

st284 = datetime(2024, 8, 8, 20, 10, 56)
st171 = datetime(2024, 8, 8, 20, 18, 16)

end284 = datetime(2024, 8, 8, 20, 38, 46)
end171 = datetime(2024, 8, 8, 20, 42, 26)

y284 = 2.2
y171 = 2.1

#lns4 = ax.hlines(y=y284, xmin=st284, xmax=end284, color=muted[5], linewidth=2,label='284 span')
#lns5 = ax.hlines(y=y171, xmin=st171, xmax=end171, color=muted[3], linewidth=2,label='171 span')

data171 = np.load('/Users/coletamburri/Desktop/August_2024_DKIST_Flares/171curve_SUVI_8Aug2024_Flare.npz');

t171 = np.arange(datetime(2024,8,8,20,6,26),
              datetime(2024,8,8,20,46,16), 
              timedelta(seconds=4*60)).astype(datetime)

data284=np.load('/Users/coletamburri/Desktop/August_2024_DKIST_Flares/284curve_SUVI_8Aug2024_Flare.npz')


t284 = np.arange(datetime(2024,8,8,20,6,56),
              datetime(2024,8,8,20,38,46), 
              timedelta(seconds=4*60)).astype(datetime)

datagong=np.load('/Users/coletamburri/Desktop/gong_curve.npz')


tgong = np.arange(datetime(2024,8,8,20,10,42),
              datetime(2024,8,8,20,34,42), 
              timedelta(seconds=60)).astype(datetime)

ax3 = ax.twinx()
lns4 = ax3.plot(t171,data171['vals']*12,'-x', color=muted[5],linewidth=2,label=r'SUVI 171$\AA$')
lns5 = ax3.plot(t284,data284['vals']*1.1,'--',color=muted[7], linewidth=2,label=r'SUVI 284$\AA$')
ax4 = ax.twinx()
lns6 = ax4.plot(tgong,datagong['vals'],'-o',markersize=2,color=muted[8], linewidth=2,label=r'GONG H$\alpha$')

ax3.set_axis_off()
ax4.set_axis_off()
ax4.set_ylim([2000000,3200000])
#ax.set_xlim([datetime(2024,8,8,20,5,32,333333),datetime(2024,8,8,21,15,7,0)])

#ax3.set_ylim([-10000,7000])
fmtr = dates.DateFormatter("%H:%M")
# need a handle to the current axes to manipulate it
ax = plt.gca()
# set this formatter to the axis
ax.xaxis.set_major_formatter(fmtr)

lns = lns1+lns2+lns3+lns4+lns5+lns6
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0,fontsize=10)

# ax.axvline(t171[np.where(data171['vals']==np.max(data171['vals'][:]))],linestyle='--', color=muted[5])
# ax.axvline(t284[np.where(data284['vals']==np.max(data284['vals'][1:]))],linestyle='--', color=muted[7])
# ax.axvline(stack[160],linestyle='--', color=muted[0])


#now plot SolO/PSP



fig.subplots_adjust(bottom=0.2)





