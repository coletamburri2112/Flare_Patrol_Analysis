#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 10:41:00 2025

@author: coletamburri
"""

array = mask

psds = []
for i in range(4095):
    ps_x = np.abs(np.fft.fft(array[i,:]))
    psds.append(ps_x)
    
psdsavg = np.mean(psds,0)

freqs = np.fft.fftfreq(array[i,:].size, step)

idx = np.argsort(freqs)

fig,ax=plt.subplots();ax.loglog((1/freqs[idx])*727,psdsavg[idx]);fig.show()


array  = maskneg

psdsneg = []
for i in range(4095):
    ps_x = np.abs(np.fft.fft(array[i,:]))
    psdsneg.append(ps_x)
    
psdsavgneg = np.mean(psdsneg,0)

freqs = np.fft.fftfreq(array[i,:].size, step)

idx = np.argsort(freqs)

fig,ax=plt.subplots();ax.loglog((1/freqs[idx])*727,psdsavgneg[idx]);fig.show()


array  = normalized[0:250,:]

psdsqs = []
for i in range(249):
    ps_x = np.abs(np.fft.fft(array[i,:]))
    psdsqs.append(ps_x)
    
psdsavgqs = np.mean(psdsqs,0)

freqs = np.fft.fftfreq(array[i,:].size, step)

idx = np.argsort(freqs)

fig,ax=plt.subplots();ax.loglog((1/freqs[idx])*727,psdsavgqs[idx]);fig.show()

array  = qsmask

psdsqsmask = []
for i in range(4095):
    ps_x = np.abs(np.fft.fft(array[i,:]))
    psdsqsmask.append(ps_x)
    
psdsavgqsmask = np.mean(psdsqsmask,0)

freqs = np.fft.fftfreq(array[i,:].size, step)

idx = np.argsort(freqs)

fig,ax=plt.subplots();ax.loglog((1/freqs[idx])*727,psdsavgqsmask[idx]);fig.show()

fig,ax=plt.subplots()
ax.loglog((1/freqs[idx])*727,psdsavg[idx],label='ribbon')
ax.loglog((1/freqs[idx])*727,psdsavgneg[idx],label='arcade')
#ax2 = ax.twinx()
#ax2.loglog((1/freqs[idx])*727,psdsavgqs[idx],label='qs',c='red')
ax.loglog((1/freqs[idx])*727,psdsavgqsmask[idx],label='qs',c='red')
ax.legend()
ax.axvline(24)
fig.show()