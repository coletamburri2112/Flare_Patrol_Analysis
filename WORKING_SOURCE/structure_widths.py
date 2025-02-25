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
import tol_colors as tc

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
save = 1 # save output arrays?
directory = '/Users/coletamburri/Desktop/double_loop_frame0_pre_destretch/'
time = '2024-08-08T20:12:32.333333'
if os.path.isdir(directory) == 0:
    os.mkdir(directory)
filenamesave = directory+'widths_errors.npz' # filename for output
numareas = 1 # number of areas to look at
numcuts = 20 # number of strands of interest per area
ampdir = 'neg'
note = []

#Determine mu
d = 151.68e9 # distance to the sun on 8 August source: https://theskylive.com/planetarium?objects=sun-moon-mercury-venus-mars-jupiter-saturn-uranus-neptune-pluto&localdata=40.01499%7C-105.27055%7CBoulder%20CO%20(US)%7CAmerica%2FDenver%7C0&obj=sun&h=14&m=01&date=2024-08-08#ra|9.242130505796545|dec|15.985314118209057|fov|50
solrad = 695700000

# Coordinates from DKIST are not correct, but define them anyways as a starting
# point.  Will co-align later in routine.

hpc1_arcsec = -175
hpc2_arcsec= 375

# image center
x_center = d*np.cos(hpc1_arcsec/206265)*np.sin(hpc2_arcsec/206265) # m
y_center = d*np.sin(hpc1_arcsec/206265) # m
z = solrad - d*np.cos(hpc1_arcsec/206265)*np.cos(hpc2_arcsec/206265) # m

# to mu value
rho = np.sqrt(x_center**2+y_center**2)
mu = np.sqrt(1-(rho/solrad)**2)

# Constants
spatial_samp = 0.017 # for vbi red at 656nm
arcsec_to_km = 727 # approximate arcsec to km conversion

# Arrays for coordinates of start and end
startx = []
starty = []

endx = []
endy = []

l = 0 # initialization of cut number

# Arrays for widths, errors associated with Gaussian models
if gauss2 == 1:
    width1s = []
    width2s = []
    widtherr1s = []
    widtherr2s = []
    r2s = []
    amp1s = []
    amp2s = []
else:
    widths = []
    widtherrs = []
    r2s = []
    amps = []

    

# # Clear plot memory
# plt.close('all')
# # for fits loading
# # Define path for input flare file, directory of files
# path = '/Users/coletamburri/Desktop/VBI_Destretching/'
# folder_vbi = 'AXXJL' # 8 August X-class flare decay phase
# #folder_vbi = 'BDJKM' # 11 August M-class flare decay phase
# filename='postdestretch_histomatch_dataCube.fits'
# dir_list = os.listdir(path+folder_vbi)

# # Define H-alpha file (all data) and frame to work with
# fullhalpha = fits.open(path+folder_vbi+'/'+dir_list[1])
# #fullhalpha = fits.open(path+folder_vbi+'/'+filename)
# first_frame = fullhalpha[0].data[0,:,:]

# for npz loading
path = '/Users/coletamburri/Desktop/VBI_Destretching/'
folder_vbi = 'AXXJL/' # 8 August X-class flare decay phase
filename = 'AXXJLselection_predestretch.npz'
array = np.load(path+folder_vbi+filename)['first50'] #first50 or brightening

#frame to work with
frame = array[0,:,:]

# X and Y coordinates of frame
xarr = np.arange(np.shape(frame)[0])
yarr = np.arange(np.shape(frame)[1])

# X and Y coordinates, in KM
xarr_km = xarr*spatial_samp
yarr_km = yarr*spatial_samp

# Meshgrid for plotting
XKM,YKM =np.meshgrid(xarr_km,yarr_km)

# Plot first frame
fig,ax=plt.subplots(dpi=400,figsize=(10,10))
ax.pcolormesh(frame,cmap='grey')
ax.set_aspect('equal')
ax.invert_yaxis()

# User input, click two points at the upper left and lower right corners
# of a rectangle to zoom in on, respectively.  Do this for the number of 
# features defined by numareas.
cc = plt.ginput(numareas*2,timeout = 120)

#arrays for storgage of zoom-in boxes
ylos=[]
yhis=[]
xlos=[]
xhis=[]
    
# Begin loops - first, define the number of areas to search through
for i in range(0,2*numareas,2):
    plt.close('all')
    
    # Extract coordinates
    ylo, yhi, xlo, xhi = int(cc[i][0]), int(cc[i+1][0]), int(cc[i][1]),\
        int(cc[i+1][1])
        
    ylos.append(ylo)
    yhis.append(yhi)
    xlos.append(xlo)
    xhis.append(xhi)

    # Extract zoomed-in frame from coordinates
    framezoom = frame[xlo:xhi,ylo:yhi]
    
    # Plot zoomed-in
    fig,ax=plt.subplots(dpi=300)
    ax.pcolormesh(framezoom,cmap='grey')
    ax.invert_yaxis()
    ax.set_aspect('equal')
    plt.show()
    
    # Point-and-click to define a line perpendicular to the desired feature;
    # Do this for the number of features defined by numcuts
    aa = plt.ginput(numcuts*2,timeout =120)
    
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
        
        # Use skimage to find the intensity profile along the line.
        # skimage.measure.profile_line returns a 1D array of intensity values 
        # in a directly line from (x0,y0) to (x1,y1), of length equal to the 
        # ceil of the computed length of the line (in units of pixels)
        profile = skimage.measure.profile_line(framezoom,[y0,x0],[y1,x1])
        
        # Convert the length of the skimage output to arcsec
        xdirection = np.arange(len(profile))*spatial_samp
        
        # Plot intensity profile in separate window
        fig,ax=plt.subplots(dpi=300)
        ax.plot(profile,'-x')
        plt.show()
        
        # Define the limits of the Gaussian (or 2-Gaussian) fitting.
        bb = plt.ginput(2,timeout = 120)
        
        # Extract the start and end of the fitting window
        st=int(bb[0][0])
        end=int(bb[1][0])+1
        
        # Append the start and end coordinates for future use.
        startx.append(x0)
        starty.append(y0)
        endx.append(x1)
        endy.append(y1)
        l+=1
        
        # Perform the fitting
        if gauss2==0: # If 2-Gauss model
            
            inf=np.inf
            # Initial guesses - can be pretty bad

            if ampdir == 'neg':
                p0=[-6000, np.median(xdirection[st:end]), 0.1, 0, 35000]
            elif ampdir == 'pos':
                p0=[6000, np.median(xdirection[st:end]), 0.1, 0, 35000]
            
            # Perform the fit, extract parameters (popt) and cov. matrix (pcov)
            try:
                popt,pcov = scipy.optimize.curve_fit(Gauss_func,\
                                                     xdirection[st:end+1],\
                                                         profile[st:end+1],p0=p0,
                                                         maxfev=200000)
            except RuntimeError:
                amps.append(np.nan)
                widths.append(np.nan) 
                widtherrs.append(np.nan)
                r2s.append(np.nan)
                print('RunTime Error!')
                
            residuals = profile[st:end+1] - Gauss_func(xdirection[st:end+1],\
                                                           *popt)
            ss_res = np.sum(residuals**2) # sum of squares
            ss_tot = np.sum((profile[st:end+1]-\
                                np.mean(profile[st:end+1]))**2)#total sum of sqs
                
            r_squared = 1 - (ss_res / ss_tot)
            
            # Extract fit values and errors
            amp, cent, std, slope, intercept = popt
            amp_err,cent_err,std_err,slope_err,\
                intercept_err = np.sqrt(np.diag(pcov))
            
            # Define the FWHM using the standard deviation
            fwhm=2*np.sqrt(2*np.log(2))*std
            
            # Propagation of error --> to FWHM
            fwhm_err = np.abs(fwhm*np.sqrt((std_err/std)**2))
            
            # Define the line width in KM
            width = np.abs(fwhm*arcsec_to_km)
            
            # Error in KM
            width_err = width*np.sqrt((fwhm_err/fwhm)**2)
            
            # Finer grid of x values
            xdirection_finer = np.arange(xdirection[st],xdirection[end],.001)
            
            # Plot the result, with both the original data and fitted model
            fig,ax=plt.subplots(figsize=(8,5),dpi=300)
            ax.plot(xdirection_finer, Gauss_func(xdirection_finer,*popt), '-',\
                     label=str(round(width,2))+ '$\;\pm\;$'+\
                         str(round(width_err,2))+'$\;km$',c='#882255')
            ax.scatter(xdirection[st:end], profile[st:end], \
                       label='Flux across cut',c='#009988')
            ax.set_xlabel('Position along cut')
            ax.set_ylabel('Flux')
            ax.legend(loc=2)
            
            
            # Plot the frame in one panel with the selected line, and the intensity
            # profile in the second panel along the selected line/
            fig, axes = plt.subplots(nrows=2,dpi=200)
            axes[0].imshow(framezoom,cmap='magma')
            axes[0].plot([x0, x1], [y0, y1], 'ro-')
            axes[0].axis('image')
            axes[1].scatter(xdirection[st:end], profile[st:end],10,c='red')
            axes[1].plot(xdirection_finer, Gauss_func(xdirection_finer,*popt),c='#882255')
            axes[0].set_title(str(l)+', w = '+str(int(round(width,2)))+'km')
            if save == 1:
                fig.savefig(directory+'cutdescrip'+str(l)+'_'+str(int(round(width,2)))+'km.png')
            
            # Append width and error values to arrays for output
            amps.append(amp)
            widths.append(width) 
            widtherrs.append(width_err)
            r2s.append(r_squared)
            
        # In the case of single gaussian fitting
        elif gauss2 == 1:
            inf=np.inf
            
            if ampdir == 'neg':
                # Initial parameter guesses
                p0=[-20000, np.median(xdirection[st:end])*2/3, 0.1,-20000,
                    np.median(xdirection[st:end])*4/3,0.1, 0, 35000]
            elif ampdir == 'pos':
                p0=[6000, .25, 0.1,6000,.35,0.1, 0, 35000]
            
            # Fitting - popt is output params, pcov is covariance matrix
            popt,pcov = scipy.optimize.curve_fit(double_gaussian,\
                                                 xdirection[st:end+1],\
                                                     profile[st:end+1],p0=p0,
                                                     maxfev=200000)
                
            residuals = profile[st:end+1] - double_gaussian(xdirection[st:end+1],\
                                                           *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((profile[st:end+1]-\
                                np.mean(profile[st:end+1]))**2)#total sum of sqs
                
            r_squared = 1 - (ss_res / ss_tot)
            
            # Extract fit and error values
            amp1, cent1, std1, amp2, cent2, std2, slope, intercept = popt
            amp_err1,cent_err1,std_err1,amp_err2,cent_err2,\
                std_err2,slope_err,intercept_err = np.sqrt(np.diag(pcov))
            
            # FWHM of each Gaussian in the model
            fwhm1=2*np.sqrt(2*np.log(2))*std1
            fwhm2=2*np.sqrt(2*np.log(2))*std2
            
            # Propagation of error to FWHM and width (in km) for both 2-Gs
            fwhm_err1 = np.abs(fwhm1*np.sqrt((std_err1/std1)**2))
            width1 = np.abs(fwhm1*arcsec_to_km)
            width_err1 = width1*np.sqrt((fwhm_err1/fwhm1)**2)
            
            fwhm_err2 = np.abs(fwhm2*np.sqrt((std_err2/std2)**2))
            width2 = np.abs(fwhm2*arcsec_to_km)
            width_err2 = width2*np.sqrt((fwhm_err2/fwhm2)**2)
            
            # Finer grid resolution
            xdirection_finer = np.arange(xdirection[st],xdirection[end],.001)
            
            # Plot result with data, fitted model, component Gaussians
            fig,ax=plt.subplots(figsize=(8,5),dpi=200)
            ax.scatter(xdirection[st:end], profile[st:end], \
                       label='Flux across cut',c='#009988')
            ax.plot(xdirection_finer, double_gaussian(xdirection_finer,*popt),\
                    '-',c='#882255')
            ax.plot(xdirection_finer, Gauss_func(xdirection_finer,\
                                                  *[amp1,cent1,std1,slope,\
                                                    intercept]),'--',\
                    c='#222255',markersize=5,label=str(round(width1,2))+\
                              '$\;\pm\;$'+str(round(width_err1,2))+'$\;km$')
            ax.plot(xdirection_finer, Gauss_func(xdirection_finer,\
                                                  *[amp2,cent2,std2,slope,\
                                                    intercept]),'--',\
                    c='#663333',markersize=5,label=str(round(width2,2))+\
                        '$\;\pm\;$'+str(round(width_err2,2))+\
                         '$\;km$')
            ax.set_xlabel('Position along cut [arcsec]')
            ax.set_ylabel('Flux')
            #ax[0].set_title(str(l)+', w2 = '+str(int(round(width1,2)))+\
            #                  ', w2 = '+str(int(round(width2,2)))+'km')
            ax.legend()
            
            if save == 1:
                fig.savefig(directory+'cutdescrip'+str(l)+'_'+str(int(round(width1,2)))+'_km,_'+str(int(round(width1,2)))+'_km.png')
            
            # Properly order component Gaussians and append to arrays
            # Same for errors
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
            r2s.append(r_squared)
            amp1s.append(amp1)
            amp2s.append(amp2)
            
# remove those that have failed - either for error or incorrect amplitude
if gauss2 == 0:
    for i in range(len(amps)):
        if (ampdir == 'neg' and amps[i] > 0) or (ampdir == 'pos' and amps[i] < 0):
            amps[i] = np.nan
            widths[i] = np.nan
            widtherrs[i] = np.nan
            r2s[i] = np.nan
            note.append(str(i)+' -- wrong amp')
        elif widtherrs[i] > 100:
            amps[i] = np.nan
            widths[i] = np.nan
            widtherrs[i] = np.nan
            r2s[i] = np.nan
            note.append(str(i)+' -- large error')
elif gauss2 == 1:
    for i in range(len(amp1s)):
        if (ampdir == 'neg' and (amp1s[i] > 0 or amp2s[i] > 0)) or \
            (ampdir == 'pos' and (amp1s[i] < 0 or amp2s[i] < 0)):
            amp1s[i] = np.nan
            width1s[i] = np.nan
            widtherr1s[i] = np.nan
            amp2s[i] = np.nan
            width2s[i] = np.nan
            widtherr2s[i] = np.nan
            r2s[i] = np.nan
            note.append(str(i)+' -- wrong amp')
        if widtherr1s[i] > 100 or widtherr2s[i] > 100:
            amp1s[i] = np.nan
            width1s[i] = np.nan
            widtherr1s[i] = np.nan
            amp2s[i] = np.nan
            width2s[i] = np.nan
            widtherr2s[i] = np.nan
            r2s[i] = np.nan
            note.append(str(i)+' -- large error')

            
# Save results depending on "save" switch
if save == 1:
    if gauss2 == 1:
        np.savez(filenamesave,width1s,width2s,widtherr1s,widtherr2s,\
                 startx,starty,endx,endy,r2s,amp1s,amp2s,note,time,ylos,yhis,xlos,xhis)
    elif gauss2 == 0:
        np.savez(filenamesave,widths,widtherrs,startx,starty,endx,endy,r2s,amps,
                 note,time,ylos,yhis,xlos,xhis)    
        
muted = tc.tol_cset('muted')

fig,ax=plt.subplots()
ax.errorbar(range(len(widths)),widths,yerr=widtherrs,linestyle='',fmt='o',\
            markersize=4,color=muted.indigo,ecolor=muted.rose,elinewidth=2,\
                capsize=3)
fig.show()
    
