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
    stack.append(str(result))
    
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

#x = [datetime.strptime(d,'%Y-%m-%d %H:%M:%S.%f').time() for d in stack]

muted = tc.tol_cset('muted')

data=np.load('/Users/coletamburri/Desktop/lc_save.npz')

fig,ax=plt.subplots(dpi=200)
#ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S.%f'))
#ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=20))
#ax.xaxis.set_major_locator(x.DayLocator())
ax.plot(timeshhmmss,data['arr_0']/1e4,c=muted[1],label='Full VBI FoV');
#ax1=ax.twinx();
#lns2 = ax1.plot(data['arr_1'],c=muted[0],label='Loop-top source');
ax.plot(timeshhmmss,data['arr_1']/1e4,c=muted[0],label='Loop-top source');
ax.set_ylabel('Avg. Intensity over FoV [$10^4\;DN\;pix^{-2}$]')
#ax.autofmt_xdate()
ax.set_xticks(timeshhmmss[0:-1:190],timeshhmmss[0:-1:190],rotation=45)
ax.set_xlabel('Time [UT]')
fig.subplots_adjust(bottom=0.2)





