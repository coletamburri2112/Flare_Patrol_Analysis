#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 17:31:01 2025

@author: coletamburri
"""

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
import matplotlib

from scipy.signal import convolve2d
from scipy.signal import convolve
from scipy.ndimage import gaussian_filter

# allow for pop-up... need pyqt5 installed
matplotlib.use('Qt5Agg')

# Function definitions for Gaussian fitting
def Gauss_func(x,A,mu,sigma,m,b):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))+ m*x + b

def double_gaussian( x, c1, mu1, sigma1, c2, mu2, sigma2 ,m,b):
    res =   (c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) )) \
          + (c2 * np.exp( - (x - mu2)**2.0 / (2.0 * sigma2**2.0) )) \
          + (m * x + b)
    return res
        
#### USER-DEFINED FUNCTIONS
directory = '/Users/coletamburri/Desktop/small_loop_cutline3_0/' #folder to save output images, array to
#time = '2024-08-08T20:12:32.333333' #time for observations... not really necessary
if os.path.isdir(directory) == 0:
    os.mkdir(directory)
filenamesave = directory+'widths_errors.npz' # filename for output
numareas = 3 # number of areas to look at
numcuts = 15 # number of strands of interest per area
ampdir = 'neg'
note = []

##### OPTION FOR NPZ LOADING#### 
path = '/Users/coletamburri/Desktop/VBI_Destretching/' # path for VBI file
folder_vbi = 'AXXJL/' # 8 August X-class flare decay phase is AXXJL
filename = 'AXXJLselection_predestretch.npz'
array = np.load(path+folder_vbi+filename)['first50'] #first50 or brightening - choose frame selection

##### OPTION FOR FITS LOADING#### 
# # Define path for input flare file, directory of files
# path = '/Users/coletamburri/Desktop/VBI_Destretching/'
# folder_vbi = 'AXXJL' # 8 August X-class flare decay phase
# #folder_vbi = 'BDJKM' # 11 August M-class flare decay phase
# filename='postdestretch_histomatch_dataCube.fits'
# dir_list = os.listdir(path+folder_vbi)

# # Define H-alpha file (all data) and frame to work with
# fullhalpha = fits.open(path+folder_vbi+'/'+dir_list[1])
# #fullhalpha = fits.open(path+folder_vbi+'/'+filename)
# array = fullhalpha[0].data[3,:,:]


#frame to work with
frame = array[3,:,:]

# Switches
gauss2 = 0 # double-gaussian models?
save = 1 # save output arrays?
blur = 0 # test at "worse" resolution? e.g. for comparison to BBSO

# Scaling definitions

# Determine mu
d = 151.68e9 # distance to the sun on 8 August source: https://theskylive.com/planetarium?objects=sun-moon-mercury-venus-mars-jupiter-saturn-uranus-neptune-pluto&localdata=40.01499%7C-105.27055%7CBoulder%20CO%20(US)%7CAmerica%2FDenver%7C0&obj=sun&h=14&m=01&date=2024-08-08#ra|9.242130505796545|dec|15.985314118209057|fov|50
solrad = 695700000

# Coordinates from DKIST are not exactly correct, but use them anyways as a starting
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


# Clear plot memory
plt.close('all')


halpha_samp = 0.017 #arcsec/pixel
resolution_aia_var = .6**2 #arcsec spatial sampling of sdo/aia
resolution_trace_var = .5**2 #spatial sampling of rhessi/trace
resolution_vbi_var = 0.017**2 #spatial sampling of DKIST/VBI in H-alpha filter
bbsogst_var = 0.034**2 #spatial sampling of BBSO/GST at H-alpha

if blur == 1:
   
    # subtract the vbi resolution from that of the instrument in question to get 
    # the width of the PSF we want to convolve.  Covolution of two gaussians is a 
    # Gaussian of variance equal to the sum of the *variances* of the two Gaussians
    # , so assume observations are "already" convolved with th DKIST PSF

    # var = std**2
    # var_total**2 = var_dkist**2 + var_aia**2
    pixels_psf_sig = round((np.sqrt(bbsogst_var-resolution_vbi_var))/halpha_samp)
    convolved = gaussian_filter(np.asarray(frame),pixels_psf_sig)
    frame = convolved

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
ax.pcolormesh(frame,cmap='magma')
ax.set_aspect('equal')
ax.invert_yaxis()

cc = plt.ginput(2,timeout = 120)

x0, y0, x1, y1 = cc[0][0], cc[0][1], cc[1][0], cc[1][1]

length = int(np.hypot(x1-x0,y1-y0))

x, y = np.linspace(y0, y1, length), np.linspace(x0, x1, length)

# Find the intensities nearest to the (x,y) coordinates above
# This essentially finds the intensity profile along the cut.
zi = frame[x.astype(int), y.astype(int)]

# Use skimage to find the intensity profile along the line.
# skimage.measure.profile_line returns a 1D array of intensity values 
# in a directly line from (x0,y0) to (x1,y1), of length equal to the 
# ceil of the computed length of the line (in units of pixels)
profile = skimage.measure.profile_line(frame,[y0,x0],[y1,x1])

# Convert the length of the skimage output to arcsec
xdirection = np.arange(len(profile))*spatial_samp

# Plot intensity profile in separate window
fig,ax=plt.subplots(dpi=300)
ax.plot(profile,'-')
plt.show()

bb = plt.ginput(10,timeout = 120)
        
# Extract the start and end of the fitting window
for i in range(0,len(bb),2):
    st=int(bb[i][0])
    end=int(bb[i+1][0])+1
    
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
        #ax.legend(loc=2)
        
        
        # Plot the frame in one panel with the selected line, and the intensity
        # profile in the second panel along the selected line/
        fig, axes = plt.subplots(nrows=2,dpi=200)
        axes[0].imshow(frame,cmap='magma')
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
        try:
            popt,pcov = scipy.optimize.curve_fit(double_gaussian,\
                                                 xdirection[st:end+1],\
                                                     profile[st:end+1],p0=p0,
                                                     maxfev=200000)
            
        except RuntimeError:
            width1s.append(np.nan) 
            width2s.append(np.nan) 
            widtherr1s.append(np.nan)
            widtherr2s.append(np.nan)
            r2s.append(np.nan)
            amp1s.append(np.nan)
            amp2s.append(np.nan)
            print('RunTime Error!')
            
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
        #ax.legend()
        
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
                 startx,starty,endx,endy,r2s,amp1s,amp2s,note,y0,y1,x0,x1)
    elif gauss2 == 0:
        np.savez(filenamesave,widths,widtherrs,startx,starty,endx,endy,r2s,amps,
                 note,y0,y1,x0,x1)    
        
muted = tc.tol_cset('muted')
plt.close('all')
