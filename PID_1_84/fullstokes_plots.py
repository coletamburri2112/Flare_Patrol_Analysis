#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 15:32:43 2025

@author: coletamburri
"""

import numpy as np
import matplotlib.pyplot as plt

fig,ax=plt.subplots(2,2)
for i in range(4):
    ax.flatten()[i].pcolormesh(image_data_arr_arr[i,:,:])
    
fig,ax=plt.subplots(2,2)
for i in range(4):
    
    ax.flatten()[i].plot(image_data_arr_arr[i,:,1377])
    
ax.flatten()[0].set_title('I')
ax.flatten()[1].set_title('Q')
ax.flatten()[2].set_title('U')
ax.flatten()[3].set_title('V')
    
