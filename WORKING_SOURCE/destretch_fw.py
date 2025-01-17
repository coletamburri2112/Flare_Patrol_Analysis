import sys
import gc
import numpy as np
from scipy import ndimage

'''
Author:   Friedrich Woeger
          DKIST/National Solar Observatory
          3665 Discovery Dr.
          Boulder, CO 80303
          USA
'''

def fracWin2d(tileSize, fraction):
  apodSize = int(tileSize*fraction + 0.5)

  x = np.arange(apodSize)/(apodSize)*(np.pi/2.0)
  y = np.sin(x)**2
  z = np.ones(tileSize)
  z[0:apodSize] = y
  z[tileSize-apodSize:tileSize] = np.flip(y)

  out = np.outer(z,z)

  return out

def globalTracking(imageCube, rMean=3, subTile=None):
  ##### global correlation tracking of cube on running mean
  # Track on a tile of size 256**2 in the middle 
  if subTile is None:
    subTile = [imageCube.shape[1]//2-128, imageCube.shape[2]//2-128, 256]
  window = fracWin2d(subTile[2], 0.4375)
  window/= np.mean(window)
  
  # loop through images
  for i in range(0, imageCube.shape[0]):
    # set reference - TRAILING running average
    if (rMean < imageCube.shape[0]):
      if (i-rMean >= 0):
        ref = np.mean(imageCube[i-rMean:i, subTile[0]:subTile[0]+subTile[2], subTile[1]:subTile[1]+subTile[2]], axis=0)
      elif (i > 0):
        ref = np.mean(imageCube[0:i, subTile[0]:subTile[0]+subTile[2], subTile[1]:subTile[1]+subTile[2]], axis=0)
      else:
        ref = imageCube[i, subTile[0]:subTile[0]+subTile[2], subTile[1]:subTile[1]+subTile[2]]
    else:
      ref = np.mean(imageCube[:, subTile[0]:subTile[0]+subTile[2], subTile[1]:subTile[1]+subTile[2]], axis=0)
    # current image to be worked on
    img = imageCube[i, subTile[0]:subTile[0]+subTile[2], subTile[1]:subTile[1]+subTile[2]]
    # compute normalization stats
    mRef = np.mean(ref)
    sRef = np.std(ref)
    mImg = np.mean(img)
    sImg  = np.std(img)

    # correlation, each FT scales with (# of pixel in image)
    refFT = np.fft.rfft2( (ref-mRef)/sRef * window )
    imgFT = np.fft.rfft2( (img-mImg)/sImg * window )
    # inverse FT scales with 1 / (# of pixel in image), but only once so we need to revert that
    xCorr = np.fft.fftshift( np.fft.irfft2( np.conj(imgFT)*refFT ) ) / (ref.shape[0]*ref.shape[1])
    # find integer shift
    idx = np.argmax(xCorr)
    in0 = idx // xCorr.shape[0] - xCorr.shape[0]//2
    in1 = idx %  xCorr.shape[0] - xCorr.shape[1]//2
    # print(i, in0, in1)

    if ((np.abs(in0) > xCorr.shape[0]//2) | \
        (np.abs(in1) > xCorr.shape[1]//2)):
      continue
    else:
      imageCube[i, :, :] = np.roll(imageCube[i, :, :], (in0, in1), axis=(0,1))
    
  return imageCube

# imageShape: 2D shape of image
# tileSize:   requested size of the tiles, required to be divisible by two
# overlap:    boolean describing whether tiles should overlap by half of tileSize
def segmentCoords(imageShape, tileSize, overlap=False):
  if (tileSize % 2 != 0):
    sys.exit("tileSize needs to be even.")
  # The following will (temporarily) store the number of non-overlapping tiles
  # in the numTiles array. We will calculate any offset before generating
  # the number of overlapping tiles in the numTiles array, if requested.
  # Note: the computation is done for y and x simultaneously
  numTiles  = np.asarray(imageShape)//tileSize
  # Calculate the edge offset. edgeOffsets will contain the total offset of both sides.
  edgeOffsets   = np.asarray(imageShape) - numTiles*tileSize
  # NOW create the array containing the number of overlapping tiles (=2n-1)
  if (overlap):
    numTiles += numTiles - 1
    tileOffset = tileSize//2    # each tile is offset from the next by half its size
  else:
    tileOffset = tileSize       # tiles are butting each other

  # output array of upper edge coords of each tile
  outCoords= np.empty((numTiles[0]*numTiles[1], 2), dtype='int16')
  # compute coords
  cnt = 0
  for i in range(numTiles[0]):
    for j in range(numTiles[1]):
      outCoords[cnt, :] = [edgeOffsets[0]//2 + i*tileOffset,
                           edgeOffsets[1]//2 + j*tileOffset]
      cnt += 1

  return numTiles, outCoords

# image:      2D image array
# tileSize:   requested size of the tiles, required to be divisible by two
# numTiles:   2 element array describing number of tiles in X and Y (output of segmentCoords)
# tileCoords: array containing the edge coords of each tile (output of segmentCoords)
#             yDim = xDim (square array) is assumed
def segmentImage(image, tileSize, numTiles, tileCoords):
  if (tileSize % 2 != 0):
    sys.exit("tileSize needs to be even.")
  outTiles = np.empty((numTiles[0]*numTiles[1], tileSize, tileSize), dtype='float32')
  for i in range(numTiles[0]*numTiles[1]):
    outTiles[i, :, :] = image[tileCoords[i, 0]:tileCoords[i, 0]+tileSize,
                              tileCoords[i, 1]:tileCoords[i, 1]+tileSize]
  return outTiles

# compute the per tile cross correlation between refCube tiles and tileCube tiles
# refCube:  tile cube for reference tiles - shape (nTiles, yDim, xDim)
# tileCube: tile cube for tiles to be correlated - shape (nTiles, yDim, xDim)
#           yDim = xDim (square array) is assumed
def computePerTileXCorrelation(refCube, tileCube):
  # window to reduce edge effects
  # this will bias the correlation coefficient
  refWin = fracWin2d(refCube.shape[1], 0.125)
  refWin/= np.mean(refWin)
  tileWin = fracWin2d(tileCube.shape[1], 0.4375)
  tileWin/= np.mean(tileWin)

  # compute normalization stats
  mRef = np.mean(refCube, axis=(1, 2), keepdims=True)
  sRef = np.std(refCube,  axis=(1, 2), keepdims=True)
  mTc  = np.mean(tileCube,axis=(1, 2), keepdims=True)
  sTc  = np.std(tileCube, axis=(1, 2), keepdims=True)

  # correlation, each FT scales with (# of pixel in image)
  refFT = np.fft.rfftn( (refCube-mRef)/sRef * refWin[np.newaxis, ...],  axes=(1, 2) )
  tcFT  = np.fft.rfftn( (tileCube-mTc)/sTc  * tileWin[np.newaxis, ...], axes=(1, 2) )
  # inverse FT scales with 1 / (# of pixel in image), but only once so we need to revert that
  xCube = np.fft.fftshift( np.fft.irfftn( np.conj(tcFT)*refFT, axes=(1, 2) ), axes = (1,2) ) / (refCube.shape[1]*refCube.shape[2])

  return xCube

# xCorr: 3D cross-correlation data array with shape (nTiles, yDim, xDim)
#        yDim = xDim (square array) is assumed
def interpolatePeak(xCorr):
  # output array initialized with zeros
  result = np.zeros((xCorr.shape[0], 2))
  # start interpolation loop
  for i in range(xCorr.shape[0]):
    ### find index of max value of current tile
    idx = np.argmax(xCorr[i, :, :])
    in0 = idx // xCorr.shape[1]
    in1 = idx %  xCorr.shape[1]
    ### bail if to large
    if ((np.abs(in0 - xCorr.shape[1]//2) > xCorr.shape[1]//2) | \
        (np.abs(in1 - xCorr.shape[2]//2) > xCorr.shape[2]//2)):
      result[i, :] = [0.0, 0.0]
    ### interpolation
    else:
      # y-dir
      if ((in0 > 0) & (in0 < xCorr.shape[1]-1)):
        den = xCorr[i, in0-1, in1] - 2.0*xCorr[i, in0, in1] + xCorr[i, in0+1, in1]
        s0 = 0.5*(xCorr[i, in0-1, in1] - xCorr[i, in0+1, in1])/den
        if ((not np.isfinite(s0)) | (s0 >= 1.0)):
          s0 = 0.0
      else:
        s0 = 0.0
      # x-dir
      if ((in1 > 0) & (in1 < xCorr.shape[2]-1)):
        den = xCorr[i, in0, in1-1] - 2.0*xCorr[i, in0, in1] + xCorr[i, in0, in1+1]
        s1  = 0.5*(xCorr[i, in0, in1-1] - xCorr[i, in0, in1+1])/den
        if ((not np.isfinite(s1)) | (s1 >= 1.0)):
          s1 = 0.0
      else:
        s1 = 0.0
      # now compute the 'correct' integer shift value
      in0 = in0 - xCorr.shape[1]//2
      in1 = in1 - xCorr.shape[2]//2
    result[i, :] = [in0+s0, in1+s1]

  return result

def regXcorrResults(posFrame, numTiles, filter=3):
  # remap results list to 2D array
  t00 = np.reshape(posFrame[:,0], numTiles) # y-shifts on 2D coarse grid
  t01 = np.reshape(posFrame[:,1], numTiles) # x-shifts on 2D coarse grid
  # median filter
  t10 = ndimage.median_filter(t00, size=filter)
  t11 = ndimage.median_filter(t01, size=filter)
  # gaussian filter
  #t10 = ndimage.gaussian_filter(t00, filter)
  #t11 = ndimage.gaussian_filter(t01, filter)

  return np.asarray([[t10,t11]])


def computeDistortionGrid(numTiles, tileCoords, tileSize, xCorShifts, imageShape):
  # remap results list to 2D array
  t00 = np.reshape(xCorShifts[:,0], numTiles) # y-shifts on 2D coarse grid
  t01 = np.reshape(xCorShifts[:,1], numTiles) # x-shifts on 2D coarse grid

  # size to interpolate to
  y0 = np.min(np.asarray(tileCoords)[:,0])
  y1 = np.max(np.asarray(tileCoords)[:,0])
  dy = 2./(tileSize) # overlapping tiles
  x0 = np.min(np.asarray(tileCoords)[:,1])
  x1 = np.max(np.asarray(tileCoords)[:,1])
  dx = 2./(tileSize) # overlapping tiles
  # interpolate upwards
  #   because of the way mgrid works with the dx, we will have tileSize/2-1
  #   too many values, has to be taken into account when embedding
  map = np.mgrid[0:numTiles[0]:dy, 0:numTiles[1]:dx]
  t10 = ndimage.map_coordinates(t00, map, order=2, mode='nearest')
  t11 = ndimage.map_coordinates(t01, map, order=2, mode='nearest')
  #   embed in full arrays, so that things are centered
  t20 = np.zeros(imageShape)
  t20[int(y0+tileSize//2):int(y1+tileSize//2), int(x0+tileSize//2):int(x1+tileSize//2)] = \
      t10[:-int(tileSize//2), :-int(tileSize//2)]
  t21 = np.zeros(imageShape)
  t21[int(y0+tileSize//2):int(y1+tileSize//2), int(x0+tileSize//2):int(x1+tileSize//2)] = \
      t11[:-int(tileSize//2), :-int(tileSize//2)]

  output = np.asarray([t20,t21])

  return output

def warpIm(image, distGrid):
  imgShape = image.shape
  # destretch
  grid = np.mgrid[:imgShape[0], :imgShape[1]]
  distortionGrid = grid - distGrid
  destretched = ndimage.map_coordinates(image, distortionGrid, order=3, mode='nearest')
  return destretched

def destretchSeq(imageSequence, tileSizeInput, rMean=3, globalTrack=None):

  dataCube = globalTracking(imageSequence, rMean=rMean, subTile=globalTrack) 

  for tileSize in tileSizeInput:
    numTiles, tileCoords = segmentCoords(dataCube.shape[1:], tileSize, overlap=True)

    ########## compute actual destretch setting
    # resultant array for all image
    # start with it being a copy of the original input
    # loop through all images
    for i in range(0, dataCube.shape[0]):
      # set reference - TRAILING running average
      if (rMean < dataCube.shape[0]):
        if (i-rMean >= 0):
          ref = np.mean(dataCube[i-rMean:i, :, :], axis=0)
        elif (i > 0):
          ref = np.mean(dataCube[0:i, :, :], axis=0)
        else:
          ref = dataCube[i, :, :]
      else:
        ref = np.mean(dataCube, axis=0)
      # use reference to compute shifts for individual image
      refCube  = segmentImage(ref, tileSize, numTiles, tileCoords)
      tCube    = segmentImage(dataCube[i, :, :], tileSize, numTiles, tileCoords)
      xCube    = computePerTileXCorrelation(refCube, tCube)
      posCube  = interpolatePeak(xCube)
      posCubeR = regXcorrResults(posCube, numTiles)
      curGrid  = computeDistortionGrid(numTiles, tileCoords, tileSize, posCubeR, dataCube.shape[1:])

      # apply for usage in the TRAILING running average and the next iteration
      dataCube[i, :, :] = warpIm(dataCube[i, :, :], curGrid)
      print(f"Image {i+1} of {dataCube.shape[0]} ...  ", end="\r")
    #####

  return dataCube
