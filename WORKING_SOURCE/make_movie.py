#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 09:57:36 2025

@author: coletamburri
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from skimage.exposure import match_histograms


# def save_numpy_array_as_movie(filepath, array, fps=30):
#     """
#     Saves a 3D numpy array as a movie file.

#     Args:
#         filepath (str): Path to save the movie file (e.g., 'output.avi').
#         array (numpy.ndarray): 3D numpy array representing the video frames 
#                                  (shape: (num_frames, height, width, channels)).
#         fps (int, optional): Frames per second. Defaults to 30.
#     """
#     if array.ndim != 3:
#         raise ValueError("Input array must be 4-dimensional (num_frames, height, width, channels)")

#     fourcc = cv2.VideoWriter_fourcc(*'H264')
#     height, width = array.shape[1], array.shape[2]
#     out = cv2.VideoWriter(filepath, fourcc, fps, (width, height))

#     for frame in array:
#         # Ensure the frame is in uint8 format
#         frame = (frame * 255).astype(np.uint8) if frame.max() <= 1.0 else frame.astype(np.uint8)
#         out.write(frame)  # Write the frame to the video
#     out.release()
    
infile = '/Users/coletamburri/Desktop/August_2024_DKIST_Flares/VBI_X_class/runningdifference_Xclassbrightening.fits'
#test = '/Users/coletamburri/Desktop/VBI_Destretching/postdestretch_histomatch_dataCube.fits'

outfile = '/Users/coletamburri/Desktop/August_2024_DKIST_Flares/VBI_X_class/runningdiff-movie_full.mp4'
#outfiletest = '/Users/coletamburri/Desktop/test.mp4'
fits_data = fits.open(infile)[0].data[1:]

#histogram matching first
gbref = fits_data[0]
timegb,yy,xx=fits_data.shape

dataCubeTrackedhist=[]
for idata in range(timegb):
    image=fits_data[idata,:,:]
    matched = match_histograms(image, gbref)
    dataCubeTrackedhist.append(matched)

fig, ax = plt.subplots(dpi=300)
im = ax.imshow(fits_data[0])  # Display the first frame initially

def animate(i):
    im.set_array(fits_data[i])
    return (im,)
# 722,1147
# 3303,2552
ani = animation.FuncAnimation(fig, animate, frames=len(fits_data), interval=200, blit=True)

ani.save(outfile, writer='ffmpeg')