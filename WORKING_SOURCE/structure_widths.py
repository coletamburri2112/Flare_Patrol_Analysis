#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 12:25:21 2025

@author: coletamburri
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import skimage
import scipy

# Function definitions for Gaussian fitting
def Gauss_func(x,A,mu,sigma,m,b):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))+ m*x + b

def double_gaussian( x, c1, mu1, sigma1, c2, mu2, sigma2 ,m,b):
    res =   (c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) )) \
          + (c2 * np.exp( - (x - mu2)**2.0 / (2.0 * sigma2**2.0) )) \
          + (m * x + b)
    return res
        
# Switches
gauss2 = 1 # double-gaussian models?
save = 0 # save output arrays?
filename = '/Users/coletamburri/Desktop/large_loops.npz' # filename for output
numareas = 1 # number of areas to look at
numcuts = 1 # number of strands of interest per area

# Constants
spatial_samp = 0.017 # for vbi red at 656nm
arcsec_to_km = 725 # approximate arcsec to km conversion

# Arrays for coordinates of start and end
startx = []
starty = []

endx = []
endy = []

# Arrays for widths, errors associated with Gaussian models
if gauss2 == 1:
    width1s = []
    width2s = []
    widtherr1s = []
    widtherr2s = []
else:
    widths = []
    widtherrs = []

# Clear plot memory
plt.close('all')

# Define path for input flare file, directory of files
path = '/Users/coletamburri/Desktop/VBI_Destretching/'
folder_vbi = 'AXXJL' # 8 August X-class flare decay phase
#folder_vbi = 'BDJKM' # 11 August M-class flare decay phase
filename='postdestretch_histomatch_dataCube.fits'
dir_list = os.listdir(path+folder_vbi)

# Define H-alpha file (all data) and frame to work with
fullhalpha = fits.open(path+folder_vbi+'/'+dir_list[1])
#fullhalpha = fits.open(path+folder_vbi+'/'+filename)
first_frame = fullhalpha[0].data[0,:,:]

# X and Y coordinates of frame
xarr = np.arange(np.shape(first_frame)[0])
yarr = np.arange(np.shape(first_frame)[1])

# X and Y coordinates, in KM
xarr_km = xarr*spatial_samp
yarr_km = yarr*spatial_samp

# Meshgrid for plotting
XKM,YKM =np.meshgrid(xarr_km,yarr_km)

# Plot first frame
fig,ax=plt.subplots(dpi=300,figsize=(10,10))
ax.pcolormesh(first_frame,cmap='grey')
ax.set_aspect('equal')
ax.invert_yaxis()

# User input, click two points at the upper left and lower right corners
# of a rectangle to zoom in on, respectively.  Do this for the number of 
# features defined by numareas.
cc = plt.ginput(numareas*2,timeout = 40)
    
# Begin loops - first, define the number of areas to search through
for i in range(0,2*numareas,2):
    
    # Extract coordinates
    ylo, yhi, xlo, xhi = int(cc[i][0]), int(cc[i+1][0]), int(cc[i][1]),\
        int(cc[i+1][1])

    # Extract zoomed-in frame from coordinates
    framezoom = first_frame[xlo:xhi,ylo:yhi]
    
    # Plot zoomed-in
    fig,ax=plt.subplots()
    ax.pcolormesh(framezoom,cmap='grey')
    ax.invert_yaxis()
    ax.set_aspect('equal')
    plt.show()
    
    # Point-and-click to define a line perpendicular to the desired feature;
    # Do this for the number of features defined by numcuts
    aa = plt.ginput(numcuts*2,timeout =40)
    
    for j in range(0,2*numcuts,2):
                
        # Extract line coordinates
        x0, y0, x1, y1 = aa[j][0], aa[j][1], aa[j+1][0], aa[j+1][1]
        
        # Define the length of the line, in pixels
        length = int(np.hypot(x1-x0,y1-y0))
        
        # Define the x and y in # pixels along the cut
        x, y = np.linspace(y0, y1, length), np.linspace(x0, x1, length)
        
        # Find the intensities nearest to the (x,y) coordinates above
        # This essentially finds the intensity profile along the cut.
        zi = framezoom[x.astype(int), y.astype(int)]
        
        # Plot the frame in one panel with the selected line, and the intensity
        # profile in the second panel along the selected line/
        fig, axes = plt.subplots(nrows=2)
        axes[0].imshow(framezoom)
        axes[0].plot([x0, x1], [y0, y1], 'ro-')
        axes[0].axis('image')
        axes[1].plot(zi)
        plt.show()
        
        # Use skimage to find the intensity profile along the line.
        # skimage.measure.profile_line returns a 1D array of intensity values 
        # in a directly line from (x0,y0) to (x1,y1), of length equal to the 
        # ceil of the computed length of the line (in units of pixels)
        profile = skimage.measure.profile_line(framezoom,[y0,x0],[y1,x1])
        
        # Convert the length of the skimage output to arcsec
        xdirection = np.arange(len(profile))*spatial_samp
        
        # Plot intensity profile in separate window
        fig,ax=plt.subplots()
        ax.plot(profile,'-x')
        plt.show()
        
        # Define the limits of the Gaussian (or 2-Gaussian) fitting.
        bb = plt.ginput(2,timeout = 20)
        
        # Extract the start and end of the fitting window
        st=int(bb[0][0])
        end=int(bb[1][0])+1
        
        # Append the start and end coordinates for future use.
        startx.append(x0)
        starty.append(y0)
        endx.append(x1)
        endy.append(y1)
        
        if gauss2==0:
            p0=[-6000, 0.3, 0.1, 0, 35000]
            
            popt,pcov = scipy.optimize.curve_fit(Gauss_func,xdirection[st:end+1],profile[st:end+1],p0=p0)
            
            perr = np.sqrt(np.diag(pcov))
            
            amp, cent, std, slope, intercept = popt
            amp_err,cent_err,std_err,slope_err,intercept_err = np.sqrt(np.diag(pcov))
            
            fwhm=2*np.sqrt(2*np.log(2))*std
            
            #propagation of error
            fwhm_err = np.abs(fwhm*np.sqrt((std_err/std)**2))
            
            width = np.abs(fwhm*arcsec_to_km)
            
            width_err = width*np.sqrt((fwhm_err/fwhm)**2)
            
            xdirection_finer = np.arange(xdirection[st],xdirection[end],.001)
            
            fig,ax=plt.subplots(figsize=(8,5),dpi=200)
            ax.plot(xdirection_finer, Gauss_func(xdirection_finer,*popt), '-',\
                     label=str(round(width,2))+ '$\;\pm\;$'+\
                         str(round(width_err,2))+'$\;km$',c='#882255')
            ax.scatter(xdirection[st:end], profile[st:end], label='Flux across cut',\
                       c='#009988')
            ax.set_xlabel('Position along cut')
            ax.set_ylabel('Flux')
            ax.legend(loc=2)
            
            widths.append(width) 
            widtherrs.append(width_err)
            
            
            
        elif gauss2 == 1:
            p0=[-6000, 0.25, 0.1,-6000,.35,0.1, 0, 35000]
            
            popt,pcov = scipy.optimize.curve_fit(double_gaussian,xdirection[st:end+1],profile[st:end+1],p0=p0)
            
            perr = np.sqrt(np.diag(pcov))
            
            amp1, cent1, std1, amp2, cent2, std2, slope, intercept = popt
            amp_err1,cent_err1,std_err1,amp_err2,cent_err2,\
                std_err2,slope_err,intercept_err = np.sqrt(np.diag(pcov))
            
            fwhm1=2*np.sqrt(2*np.log(2))*std1
            fwhm2=2*np.sqrt(2*np.log(2))*std2
            
            #propagation of error
            fwhm_err1 = np.abs(fwhm1*np.sqrt((std_err1/std1)**2))
            width1 = np.abs(fwhm1*arcsec_to_km)
            width_err1 = width1*np.sqrt((fwhm_err1/fwhm1)**2)
            
            fwhm_err2 = np.abs(fwhm2*np.sqrt((std_err2/std2)**2))
            width2 = np.abs(fwhm2*arcsec_to_km)
            width_err2 = width2*np.sqrt((fwhm_err2/fwhm2)**2)
            
            xdirection_finer = np.arange(xdirection[st],xdirection[end],.001)
            fig,ax=plt.subplots(figsize=(8,5),dpi=200)
            
            ax.scatter(xdirection[st:end], profile[st:end], label='Flux across cut',\
                       c='#009988')
            ax.plot(xdirection_finer, double_gaussian(xdirection_finer,*popt), '-',\
                    c='#882255')
        
            ax.plot(xdirection_finer, Gauss_func(xdirection_finer,\
                                                  *[amp1,cent1,std1,slope,intercept]),\
                     '--',c='#222255',markersize=5,label=str(round(width1,2))+\
                              '$\;\pm\;$'+str(round(width_err1,2))+'$\;km$')
            ax.plot(xdirection_finer, Gauss_func(xdirection_finer,\
                                                  *[amp2,cent2,std2,slope,intercept]),\
                     '--',c='#663333',markersize=5,label=str(round(width2,2))+ '$\;\pm\;$'+str(round(width_err2,2))+\
                         '$\;km$')
            ax.set_xlabel('Position along cut [arcsec]')
            ax.set_ylabel('Flux')
            ax.legend()
            
            if width1 < width2:
                smallerw = width1
                biggerw = width2
                smallerwerr = width_err1
                biggerwerr = width_err2
            elif width2 < width1:
                smallerw = width2
                biggerw = width1
                smallerwerr = width_err2
                biggerwerr = width_err1
                
            width1s.append(smallerw) 
            width2s.append(biggerw) 
            
            widtherr1s.append(smallerwerr)
            widtherr2s.append(biggerwerr)
            
if save == 1:
    if gauss2 == 1:
        np.savez(filename,'widths1','widths2','widtherr1s','widtherr2s',\
                 'startx','starty','endx','endy')
    elif gauss2 == 0:
        np.savez(filename,'widths','widtherrs','startx','starty','endx','endy')        
    
