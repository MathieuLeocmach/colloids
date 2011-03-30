#    Copyright 2009 Mathieu Leocmach
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
from __future__ import with_statement #for python 2.5, useless in 2.6
import numpy as np
import os.path, subprocess, shlex
from colloids import lif, vtk
from scipy.ndimage.filters import gaussian_filter, sobel
from scipy.ndimage.morphology import grey_erosion, grey_dilation


def bestParams(inputPath, outputPath, radMins=np.arange(2.5, 5, 0.5), radMaxs=np.arange(6, 32, 4), t=0, serie=None):
    if serie==None:
        serie = lif.Reader(inputPath).chooseSerieIndex()
    s = lif.Reader(inputPath).getSeries()[serie]
    out_head = outputPath+'_thr0_radMin%(m)g_radMax%(M)g_'
    if s.getNbFrames>1:
        dt = s.getTimeLapse()
        dt_min, dt_sec = divmod(dt, 60)
        dt_sec, dt_msec = divmod(dt_sec, 1)
        if dt_min > 0:
            out_head += '%dmin_'%dt_min
            dt_msec=0
        if dt_sec > 1:
            out_head += '%ds'%dt_sec
            if dt_msec > 0.01:
                out_head += '%d'%np.floor(dt_msec*100)
            out_head +='_'
    out_noext = out_head+'t%(t)0'+('%d'%len('%d'%s.getNbFrames()))+'d'
    for radMin in radMins:
	for radMax in radMaxs:
            if not os.path.exists(
                (out_noext+'.intensity')%{'m':radMin, 'M':radMax, 't':t}
                ):
                p = subprocess.Popen(
                    shlex.split(
                        (
                            ('/home/mathieu/test/bin/tracker -L -i %(lif)s -o %(out)s --serie %(serie)d'%{
                                'lif':inputPath,
                                'out':out_head,
                                'serie':serie,
                                }) +
                            ' --radiusMin %(m)g --radiusMax %(M)g --threshold 0 --onlyTimeStep %(t)d --exportIntensities --fftPlanning 32'
                        )%{
                            'm':radMin,
                            'M':radMax, 't':t
                            }
                        ),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                    )
                out, err = p.communicate()
                if not p.returncode==0:
                    print 'At radMin=%(m)g and radMax=%(M)g'%{'m':radMin, 'M':radMax}
                    print out
                    print err
                    return
    hists = np.asarray([[
	np.histogram(
		np.loadtxt(
                    (out_noext+'.intensity')%{'m':radMin, 'M':radMax, 't':t}
                    ),
		bins=range(300)
		)[0]
	for radMax in radMaxs] for radMin in radMins])
    return hists

def getBandPassMask(serie, radMin, radMax):
    maskshape = serie.getFrameShape()
    maskshape[-1] = maskshape[-1]/2 + 1
    mask = np.zeros(maskshape, dtype=bool)
    freqs = np.repeat(np.expand_dims(serie.getFrameShape(),1),2, axis=1)
    freqs[:,0] /= float(2*radMin)
    freqs[:,1] /= float(2*radMax)
    freqs[-1] *= serie.getZXratio()
    xs = (np.repeat(np.arange(maskshape[0]/2, dtype=float)[:, np.newaxis],2, axis=1)/freqs[0])**2
    ys = (np.repeat(np.arange(maskshape[1]/2, dtype=float)[:, np.newaxis],2, axis=1)/freqs[1])**2
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            z = freqs[-1] * np.where(x+y<1, np.sqrt(1 - x - y), 0)
            mask[i, j, z[1]:z[0]] = True
            mask[-i, j, z[1]:z[0]] = True
            mask[i, -j, z[1]:z[0]] = True
            mask[-i, -j, z[1]:z[0]] = True
    return mask

def saveIntensityVTK(fname):
	v=vtk.Polydata()
	v.points = np.loadtxt(fname, skiprows=2)
	v.scalars = [('Intensity', np.loadtxt(os.path.splitext(fname)[0]+'.intensity'))]
	n = fname.split('_t')
	n[-2]+='_test'
	n[-1] = n[-1][:-3]+'vtk'
	v.save('_t'.join(n))

def diff_of_gaussians(im, k=1.6, n=3):
    """
Construct an octave of the scale space by the difference of gaussians [1].

 im input image. Any dimension is correct. Only floating type of sufficient precision yeild coherent output (see scipy.ndimage.filters.gaussian_filter).
 k minimum bluring.
 n number of intervals to divide the octave into.

Using this scale space, blobs can be detected between k*sqrt(2)*2**(1.0/n) and k*sqrt(2)*2

[1] David G Lowe, International Journal Of Computer Vision 60, 91-110 (2004)."""
    return np.diff([
        gaussian_filter(im, s) for s in k*(2**(1.0/n))**np.arange(n+3)
        ], axis=0)

def pixel_centers_2Dscale(im, k=1.6, n=3):
    """Blob finder : find the local maxima in an octave of the scale space"""
    assert im.ndim==2, "work only with 2D images"
    DoG2D = diff_of_gaussians(np.asarray(im, float), k, n)
    centers_scale = np.bitwise_and(
	    DoG2D[1:-1]==grey_erosion(DoG2D, [3]*3)[1:-1],
	    DoG2D[1:-1]<0)
    centers_scale[:,:,0] = 0
    centers_scale[:,:,-1] = 0
    centers_scale[:,0] = 0
    centers_scale[:,-1] = 0
    return centers_scale

def pretty_2Dscale_track(im, k=1.6, n=3):
    """Pretty display of the blobs localized in the input image

The red channel is the original image.
The green channel displays a bright pixel at the position of each blob.
The crosses in the blue channel have the size of the detected blob."""
    centers_scale = pixel_centers_2Dscale(im, k=1.6, n=3)
    rads = np.asarray(2*np.sqrt(2)*k*(2**(1.0/n))**np.arange(1, n+1), int)/2
    return np.dstack((
        im, 255 * centers_scale.max(0),
        255*np.max([
            np.bitwise_or(
                grey_dilation(d, [1,2*r+1]),
                grey_dilation(d, [2*r+1,1]))
            for r,d in zip(rads, centers_scale)], 0)
        ))

def local_disp(c, DoG):
    """Interpolate the scale-space to find the extremum with subpixel resolution"""
    m = np.maximum(np.zeros_like(c), c-2)
    M = np.minimum(DoG.shape, c+3)
    ngb = DoG[m[0]:M[0], m[1]:M[1], m[2]:M[2]]
    grad = np.asarray([sobel(ngb, axis=a) for a in range(ngb.ndim)])
    hess = np.asarray([
            [sobel(ga, axis=a) for ga in grad] for a in range(ngb.ndim)])
    s, y, x = c-m
    dc = -4**(ngb.ndim-1)*np.dot(
            np.linalg.inv(hess[:,:,s,y,x]),
            grad[:,s,y,x])
    value = DoG[s,y,x]+0.5*np.dot(grad[:,s,y,x],dc)
    return dc, value

def centers_2Dscale(im, k=1.6, n=3):
    """Blob finder : find the local maxima in the scale space"""
    assert im.ndim==2, "work only with 2D images"
    DoG2D = diff_of_gaussians(np.asarray(im, float), k, n)
    centers_scale = np.bitwise_and(
	    DoG2D==grey_erosion(DoG2D, [3]*3),
	    DoG2D<0)
    #remove maxima on the borders
    centers_scale[:,:,0] = 0
    centers_scale[:,:,-1] = 0
    centers_scale[:,0] = 0
    centers_scale[:,-1] = 0
    centers_scale[0] = 0
    centers_scale[-1] = 0
    #from array to coordinates
    centers = np.transpose(np.where(centers_scale))
    #subpixel resolution (first try)
    dcenval = [local_disp(c, DoG2D) for c in centers]
    dcenters = np.asarray([d[0] for d in dcenval])
    vals = np.asarray([d[1] for d in dcenval])
    #if the displacement is larger than 0.5 in any direction,
    #the center is shifted to the neighbouring pixel
    for p in range(len(dcenters)):
	if np.absolute(dcenters[p]).max()>0.5:
		nc = (centers[p]+(dcenters[p]>0.5))-(dcenters[p]<-0.5)
		ndc, nv = local_disp(nc, DoG0)
		#remove the center if it is moving out of its new pixel (unstable)
		if np.absolute(ndc).max()>0.5:
			centers[p] = -1
			continue
		centers[p] = nc
		dcenters[p] = ndc
		vals[p] = nv
    
    return (centers+dcenters)[np.where(np.bitwise_and(
        centers[:,0]>-1,
        vals<0))[0]]


