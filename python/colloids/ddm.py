#!/usr/bin/env python
# -*- coding: utf-8 -*-

#    Copyright 2015 Mathieu Leocmach
#
#    This file is part of Colloids.
#
#    Colloids is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Colloids is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Colloids.  If not, see <http://www.gnu.org/licenses/>.
#

# Implementation of the basics of Differential Dynamics Microscopy
# Cerbino, R. & Trappe, V. Differential dynamic microscopy: Probing wave vector dependent dynamics with a microscope. Phys. Rev. Lett. 100, 1–4 (2008).

import numpy as np
from matplotlib.image import imread
import numexpr

def readImages(path, N, t0=1):
    """Read N successive images using the given and stock them in a numpy array 
    NxHxW (image shape is HxW) of type np.uint8.
    The path is a string to be formatted like 'directory/images_{:06d}.tif'
    so that path.format(124) -> 'directory/images_000124.tif'.
    """
    #get the images shape while checking that the last image do exist
    H, W = imread(path.format(N-1+t0)).shape[:-1] #imread makes 4 channels out of one
    images = np.zeros([N, H, W], np.uint8)
    for t in range(N):
        images[t] = imread(path.format(t+t0))[:,:,0]
    return images
    
    
def spectreDiff(im0, im1):
    """Compute the squared modulus of the 2D Fourier Transform of the difference between im0 and im1"""
    return numexpr.evaluate(
            'real(abs(f))**2',
            {'f': np.fft.fft2(im1-im0.astype(float))}
            )
    
class RadialAverager(object):
    """Radial average of a 2D array centred on (0,0), like the result of fft2d."""
    def __init__(self, shape):
        assert len(shape)==2
        self.dists = np.sqrt(np.fft.fftfreq(shape[0])[:,None]**2 +  np.fft.fftfreq(shape[1])[None,:]**2)
        self.bins = np.fft.fftfreq(max(shape))[:max(shape)/2]
        self.hd = np.histogram(self.dists, self.bins)[0]
    
    def __call__(self, im):
        assert im.shape == self.dists.shape
        hw = np.histogram(self.dists, self.bins, weights=im)[0]
        return hw/self.hd

def radialAverage(im):
    """Radial average of a 2D array centred on (0,0), like the result of fft2d."""
    dists = np.sqrt(np.fft.fftfreq(im.shape[0])[:,None]**2 +  np.fft.fftfreq(im.shape[1])[None,:]**2)
    bins = np.fft.fftfreq(max(im.shape))[:max(im.shape)/2]
    hd = np.histogram(dists, bins)[0]
    hw = np.histogram(dists, bins, weights=im)[0]
    return hw/hd
    
def timeAveraged(images, dt, navmax=100):
    """Does at most navmax spectreDiff on regularly spaced couples of images. 
    Separation within couple is dt."""
    result = np.zeros(images.shape[1:])
    step = max([(len(images)-dt)/navmax, 1])
    couples = np.arange(0, len(images)-dt, step)
    #print step, len(couples)
    for t in couples:
        result += spectreDiff(images[t], images[t+dt])
    return result / len(couples)
    
def logSpaced(L, num=50):
    """Generate an array of log spaced integers smaller than L"""
    return np.unique(np.logspace(
        start=0, stop=np.log(L)/np.log(2), 
        num=num, base=2, endpoint=False
        ).astype(int))
    
def ddm(images, navmax=100, num=50):
    """Does timeAveraged and radialAverage for log-spaced time intervals.
    Returns (intervals, ddm)"""
    dts = logSpaced(len(images), num)
    ra = RadialAverager(images.shape[1:])
    D = np.zeros((len(dts), len(ra.hd)))
    for i, dt in enumerate(dts):
        D[i] = ra(timeAveraged(images, dt, navmax))
    return dts, D
