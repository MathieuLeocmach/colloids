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
import os.path, subprocess, shlex, string, re, time
from colloids import lif, vtk
from scipy.ndimage.filters import gaussian_filter, gaussian_filter1d, sobel, uniform_filter
from scipy.ndimage.morphology import grey_erosion, grey_dilation, binary_dilation
from scipy.ndimage import measurements
from scipy.sparse.linalg import splu, spsolve
from scipy import sparse
from scipy import weave
from scipy.weave import converters
import numexpr
import PIL.Image
import unittest

coefprime = np.array([1,-8, 0, 8, -1])
coefsec = np.array([-1, 16, -30, 16, -1])

def readTIFF16(path, bigendian=True):
    """Read 16bits TIFF"""
    im = PIL.Image.open(path)
    out = np.fromstring(
        im.tostring(), 
        np.uint8
        ).reshape(tuple(list(im.size)+[2]))
    if bigendian:
        return (np.array(out[:,:,0], np.uint16) << 8) + out[:,:,1]
    else:
        return (np.array(out[:,:,1], np.uint16) << 8) + out[:,:,0]

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
    centers_scale = np.bitwise_and(np.bitwise_and(
	    DoG2D[1:-1]==grey_erosion(DoG2D, [7]*3)[1:-1],
	    DoG2D[1:-1]<0), grey_dilation(DoG2D, [3]*3)[1:-1]<0)
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

class LocalDispl:
    def __init__(self, space):
        self.space = space
        self.ngb = np.zeros([5]*self.space.ndim, float)
        self.grad = np.zeros([self.space.ndim]+[3]*self.space.ndim, float)
        self.hess =  np.zeros([self.space.ndim]*2, float)

    def __call__(self, c):
        m = np.maximum(np.zeros_like(c), c-2)
        M = np.minimum(self.space.shape, c+3)
        if np.min(m==c-2) and np.min(M==c+3):
            #not on a border, thus quick calculation
            self.ngb[:] = self.space[tuple([slice(u,U) for u, U in zip(m,M)])]
            for a in range(self.space.ndim):
                self.grad[a] = sobel(self.ngb, axis=a)[tuple([slice(1,-1)]*self.space.ndim)]
            for a in range(self.space.ndim):
                for b in range(a, self.space.ndim):
                    self.hess[a,b] = sobel(self.grad[a], axis=b)[tuple([1]*self.space.ndim)]
                    self.hess[b,a] = self.hess[a,b]
            gr = self.grad[tuple([Ellipsis]+[1]*self.space.ndim)]
            dc = -4**(self.space.ndim-1)*np.dot(np.linalg.inv(self.hess), gr)
            value = self.ngb[tuple([2]*self.space.ndim)]+0.5*np.dot(gr,dc)
            return dc, value
        else:
            #on a border, use errorproof but slow scipy code
            return local_disp(c, self.space)
            
        
        
def local_disp(c, DoG):
    """Interpolate the scale-space to find the extremum with subpixel resolution"""
    m = np.maximum(np.zeros_like(c), c-2)
    M = np.minimum(DoG.shape, c+3)
    ngb = DoG[tuple([slice(u,U) for u, U in zip(m,M)])]
    grad = np.asarray([sobel(ngb, axis=a) for a in range(ngb.ndim)])
    hess = np.asarray([
            [sobel(ga, axis=a) for ga in grad] for a in range(ngb.ndim)])
    sl = tuple([Ellipsis]+(c-m).tolist())
    dc = -4**(ngb.ndim-1)*np.dot(
            np.linalg.inv(hess[sl]),
            grad[sl])
    value = ngb[sl]+0.5*np.dot(grad[sl],dc)
    return dc, value

def elongated(c, DoG, elong_max=10):
    """Tell if a point of interest describes an elongated structure"""
    assert DoG.ndim==3, """Works only for 2D images"""
    m = np.maximum(np.zeros_like(c), c-2)[1:]
    M = np.minimum(DoG.shape, c+3)[1:]
    ngb = DoG[tuple([c[0]]+[slice(u,U) for u, U in zip(m,M)])]
    grad = np.asarray([sobel(ngb, axis=a) for a in range(ngb.ndim)])
    hess = np.asarray([
            [sobel(ga, axis=a) for ga in grad] for a in range(ngb.ndim)])
    sl = tuple([Ellipsis]+(c[1:]-m).tolist())
    det = np.linalg.det(hess[sl])
    if det<0:
        #the point is not a maximum
        return True
    return np.trace(hess[sl])**2/det > (elong_max+1.0)**2/elong_max
    

def find_blob(im, k=1.6, n=3):
    """Blob finder : find the local maxima in an octave of the scale space"""
    DoG = diff_of_gaussians(np.asarray(im, float), k, n)
    centers_scale = np.bitwise_and(
	    DoG==grey_erosion(DoG, [3]*DoG.ndim),
	    DoG<0)
    #remove maxima on the borders
    for a in range(centers_scale.ndim):
        centers_scale[tuple([slice(None)]*(centers_scale.ndim-1-a)+[0])]=0
        centers_scale[tuple([slice(None)]*(centers_scale.ndim-1-a)+[-1])]=0
    #from array to coordinates
    centers = np.transpose(np.where(centers_scale))
    if len(centers)==0:
        return np.zeros([0, DoG.ndim+1])
    #subpixel resolution (first try)
    dcenval = [local_disp(c, DoG) for c in centers]
    dcenters = np.asarray([d[0] for d in dcenval])
    vals = np.asarray([d[1] for d in dcenval])
    #if the displacement is larger than 0.5 in any direction,
    #the center is shifted to the neighbouring pixel
    for p in range(len(dcenters)):
        if np.absolute(dcenters[p]).max()>0.5:
            nc = (centers[p]+(dcenters[p]>0.5))-(dcenters[p]<-0.5)
            ndc, nv = local_disp(nc, DoG)
            #remove the center if it is moving out of its new pixel (unstable)
            if np.absolute(ndc).max()>0.5:
                centers[p] = -1
                continue
            centers[p] = nc
            dcenters[p] = ndc
            vals[p] = nv
    #filter out the unstable or dim centers
    good = np.bitwise_and(centers[:,0]>-1, vals<0)
    if not good.max():
        return np.zeros([0, DoG.ndim+1])
    centers = centers[good]
    dcenters = dcenters[good]
    #remove elongated features (2D image only)
##    if im.ndim == 2:
##        good = [not elongated(c, DoG, elong_max) for c in centers]
##        centers = centers[good]
##        dcenters = dcenters[good]
    return np.column_stack((centers+dcenters, vals[good]))

def save_2Dcenters(file_pattern, frame):
    for z, im in enumerate(frame):
	centers = np.vstack((
		find_blob(np.repeat(np.repeat(im, 2, 0), 2,1))*[1, 0.5, 0.5,1],
		find_blob(im)+[3,0,0,0]
	))[:,[2,1,0,3]]
	centers[:,2] = 1.6/np.sqrt(2)*2**(centers[:,2]/3)
	np.savetxt(
		file_pattern%(z, 'dat'),
		np.vstack((
			[1, len(centers), 1],
			[256, 256, 1.6*np.sqrt(2)*2**(7.0/3)],
			centers[:,:-1]
			)), fmt='%g'
		)
	np.savetxt(file_pattern%(z, 'intensity'), centers[:,-1], fmt='%g')
    

def load_clusters(trajfile):
    """Load the piles of 2D centers as formed by linker"""
    clusters = []
    path = os.path.split(trajfile)[0]
    with open(trajfile) as f:
        #parsing header
        ZXratio = float(f.next()[:-1].split("\t")[1])
        pattern, ext = os.path.splitext(f.next()[:-1].split("\t")[0])
        token, = f.next()[:-1].split("\t")
        offset,size = map(int, f.next()[:-1].split("\t"))
        m = re.match('(.*)'+token+'([0-9]*)',pattern)
        head = m.group(1)
        digits = len(m.group(2))
        recomposed = os.path.join(path, head+token+'%0'+str(digits)+'d%s')
        #load coordinates in the (x, y, r) space for each z
        slices = [np.hstack((
            np.loadtxt(recomposed%(z,ext), skiprows=2),
            -np.loadtxt(recomposed%(z,'.intensity'))[:, np.newaxis]
            ))
            for z in range(size)]
        #parse cluster description
        for line in f:
            z0 = int(line[:-1])
            pos = map(int, string.split(f.next()[:-1],'\t'))
            clusters.append(np.asarray([
                [s[p,0], s[p,1], z*ZXratio,
                 s[p,-2],
                 s[p,-1]] for z, s, p in zip(
                    range(z0, z0+len(pos)),
                    slices[z0:],
                    pos)
                ]))
    return clusters
    
def clusters2particles(clusters, k=1.6, n=3, noDuplicate=True, outputFinders=False):
    particles = []
    #allocate the memory for each possible size of MultiscaleBlobFinder
    N = max(map(len, clusters))
    finders = [MultiscaleBlobFinder([l], n, 3) for l in range(5, N+1)]
    for cl in clusters:
        #a cluster that appears in less than 3 slices is not a real particle
        if len(cl)<3:
            continue
        #No blob can be extracted out of too short signal
        #we just take the middle of the cluster
        if len(cl)<5:
            m = np.argmax(cl[:,-2])
            if m>0 and m+1<len(cl):
                blobs = [np.asarray([[
                    measurements.center_of_mass(cl[m-1:m+2,-2])[0]+m-1,
                    cl[0.5*len(cl),-2], 0.0
                    ]])]
            else:
                blobs = [np.asarray([[
                    measurements.center_of_mass(cl[:,-2])[0],
                    cl[0.5*len(cl),-2], 0.0
                    ]])]
        else:
            #xy gradient along z
            gradxy = np.sqrt(sobel(cl[:,0], axis=0)**2 + sobel(cl[:,1], axis=0)**2)
            #try to split the cluster at the location of strong xy(z) gradient
            gradblobs = finders[len(cl)-5](gradxy, k)
            #keep only strong position change
            if len(gradblobs):
                gradblobs = gradblobs[gradblobs[:,-1]<-0.05]
            if len(gradblobs):
                #split into subclusters
                clusters += np.array_split(cl, np.sort(
                    np.rint(gradblobs[:,0]).astype(int)
                    ))
                #stop the treatment of the actual cluster.
                #subclusters will be treated further in the loop
                continue
            #get the blobs for each signal, remove signals without blob
            blobs = filter(len, [
                finders[len(cl)-5](u, k) for u in [cl[:,3], cl[:,4]]
                ])
            if len(blobs)==0:
                #no blob in any signal, we just take the position of the maximum radius
                m = np.argmax(cl[:,-2])
                if m>0 and m+1<len(cl):
                    blobs = [np.asarray([[
                        measurements.center_of_mass(cl[m-1:m+2,-2])[0]+m-1,
                        cl[0.5*len(cl),-2], 0.0
                        ]])]
                else:
                    blobs = [np.asarray([[
                        measurements.center_of_mass(cl[:,-2])[0],
                        cl[0.5*len(cl),-2], 0.0
                        ]])]
        #sort the blob of each signal by intensity (negative)
        blobs = np.vstack([bs[np.argsort(bs[:,-1])] for bs in blobs])
        #Remove overlapping centers than may appear at different scales
        #in the different signals.
        #The most intense blob in apparent radius is the best,
        #then, the second intense in apparent radius, etc.
        #then, the blobs in intensity
        if noDuplicate:
            #The blobs in a stack of 2D centers are few (<30),
            #a simple O(N**2) algorithm is enough
            out = []
            for i in blobs:
                for j in out:
                    if (i[0]-j[0])**2 < (i[1]+j[1])**2:
                        break
                else:
                    out.append(i)
            blobs = np.vstack(out)
        #Interpolate the x,y,r,intensity values at the position of each blob
        grad = gaussian_filter1d(cl, k/2, axis=0, order=1)
        for z, s, v in blobs:
            #smoothed = (gaussian_filter1d(cl, 1.6*2**(s/3-1), axis=0)-cl.mean(0))*np.sqrt(2)+cl.mean(0)
            #grad = gaussian_filter1d(cl, 1.6*2**(s/3-1), axis=0, order=1)
            zi = np.rint(z)
            if zi==len(cl):
                zi = len(cl)-1
            if zi<0:
                zi = 0
            dz = z-zi
            particles.append(cl[zi]+grad[zi]*dz)
    if outputFinders:
        return np.asarray(particles), finders
    else:
        return np.asarray(particles)


def label_clusters(centers):
    """Group 2D centers by proximity.
centers is a (N,3) array of the coordinates"""
    dom = np.zeros(centers.max(0)+1, bool)
    dom[tuple(np.asarray(centers.T, int))] = True
    lab, nlab = measurements.label(dom, np.ones([3]*3))
    clusters = [[] for l in range(nlab+1)]
    for p, (x,y,z) in enumerate(centers):
        clusters[lab[x, y, z]].append(p)
    return clusters

def split_clusters(clusters, centers):
    for i, k in enumerate(clusters):
        if len(k)==0:
            continue
        cl = centers[k]
        M = cl[cl[:,-1].argmax()]
        core = np.sum((cl[:,:-1]-M[:-1])**2, -1)<(cl[:,-1]+M[-1])**2
        if core.sum()==len(core):
            continue
        out = [p for p,c in zip(k, core) if not c]
        clusters.append(out)
        clusters[i] = [p for p,c in zip(k, core) if c]
    return clusters
    
def draw_spheres(shape, pos, radii):
    im = np.zeros(shape, np.uint8)
    if np.isscalar(radii):
        radii = radii * np.ones(len(pos))
    assert len(pos)==len(radii)
    assert np.min(pos.max(0)<shape[::-1]) and np.min([0,0,0]<=pos.min(0)), "points out of bounds"
    code = """
    #pragma omp parallel for
    for(int p=0; p<Npos[0]; ++p)
    {
        const double rsq = pow(radii(p), 2);
        for(int i = std::max(0, (int)(pos(p,2)-radii(p))); i<std::min(Nim[0], (int)(pos(p,2)+radii(p))+1); ++i)
        {
            const double di = pow(pos(p,2)-i, 2);
            for(int j = std::max(0, (int)(pos(p,1)-radii(p))); j<std::min(Nim[1], int(pos(p,1)+radii(p))+1); ++j)
            {
                const double dj = pow(pos(p,1)-j, 2);
                for(int k = std::max(0, (int)(pos(p,0)-radii(p))); k<std::min(Nim[2], (int)(pos(p,0)+radii(p))+1); ++k)
                {
                    const double dsq = pow(pos(p,0)-k, 2) + dj + di;
                    if(dsq<=rsq)
                        im(i,j,k) = 255;
                }
            }
        }
    }
    """
    weave.inline(
        code,['pos', 'radii', 'im'],
        type_converters = converters.blitz,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return im;
    
def draw_rods(pos, radii, shape=None, im=None):
    if im is None:
        im = np.zeros(shape, bool)
    if shape is None:
        shape = im.shape
    assert len(pos)==len(radii)
    assert np.min(pos.max(0)<shape[::-1]) and np.min([0,0,0]<=pos.min(0)), "points out of bounds"
    code = """
    #pragma omp parallel for
    for(int p=0; p<Npos[0]; ++p)
    {
        const int i = pos(p,0), j = pos(p,1);
        im(blitz::Range(
            std::max(0, (int)(pos(p,2)-radii(p))),
            std::min(Nim[0], (int)(pos(p,2)+radii(p))+1)-1
            ), j, i) = 1;
    }
    """
    weave.inline(
        code,['pos', 'radii', 'im'],
        type_converters = converters.blitz,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return im;
    
    
def get_deconv_kernel(im, k=1.6, pxZ = 1.0, pxX=1.0):
    """Compute the deconvolution kernel from a priori isotropic image. 
    Returned kernel is in Fourier space."""
    assert im.ndim == 3
    imbl = gaussian_filter(im, k)
    sblx = (np.abs(np.fft.rfft(imbl, axis=2))**2).mean(0).mean(0)
    sblz = (np.abs(np.fft.rfft(imbl, axis=0))**2).mean(1).mean(1)
    f2 = np.interp(
        np.fft.fftfreq(2*len(sblz), pxZ)[:len(sblz)], 
        np.fft.fftfreq(2*len(sblx), pxX)[:len(sblx)], sblx
        )/sblz
    return np.sqrt(f2)
    
def deconvolve(im, kernel):
    """Deconvolve the input image. Suppose no noise (already blurred input)."""
    return np.fft.irfft(np.fft.rfft(im, axis=0) * kernel[:,None,None], axis=0, n=im.shape[0])

def radius2scale(R, k=1.6, n=3.0, dim=3):
    """Converts a radius (in pixels) to a scale index (logarithmic scale corresponding to the inner working of a MultiscaleTracker)"""
    return n / np.log(2) * np.log(
        R/k * np.sqrt(n*(2**(2.0/n) - 1)/(2 * dim* np.log(2)))
        ) -1

def scale2radius(x, k=1.6, n=3.0, dim=3):
    """Converts a scale index (logarithmic scale corresponding to the inner working of a MultiscaleTracker) to a radius (in pixels)"""
    return k * 2**((x+1)/float(n))*np.sqrt(
        2 * dim * np.log(2) / float(n) / (2**(2.0/n)-1)
        )
        
def radius2sigma(R, n=3.0, dim=3):
    """Converts a radius (in pixels) to a scale (logarithmic pixel scale)"""
    return R / np.sqrt(2*dim* np.log(2) / n/(1 - 2**(-2.0/n)))
    
def sigma2radius(sigma, n=3.0, dim=3):
    """Converts a scale (logarithmic pixel scale) to a radius (in pixels)"""
    return sigma * np.sqrt(2*dim* np.log(2) / n/(1 - 2**(-2.0/n)))
        
support_functions = """
#include <boost/math/special_functions/erf.hpp>
    inline double halfG(const double &d, const double &R, const double &sigma)
    {
        return exp(-pow(R+d, 2)/(2*pow(sigma,2)))*sqrt(2/M_PI)*sigma/d + boost::math::erf((R+d)/sigma/sqrt(2));
    }
    inline double G(const double &d, const double &R, const double &sigma)
    {
        return halfG(d,R,sigma) + halfG(-d, R, sigma);
    }
    inline double G(const double &R, const double &sigma)
    {
        const double x = R/sigma/sqrt(2);
        return boost::math::erf(x) - x*exp(-pow(x, 2))*2/sqrt(M_PI);
    }
    inline double DoG(const double &d, const double &R, const double &sigma, const double &alpha)
    {
        return G(d,R, alpha*sigma) - G(d, R, sigma);
    }
    inline double halfG_dsigma(const double &d, const double &R, const double &s2)
    {
        return (R*R+d*R+s2)*exp(-pow(R+d, 2)/(2*s2))/sqrt(2*M_PI)/d/s2;
    }
    inline double G_dsigma(const double &d, const double &R, const double &s2)
    {
        return halfG_dsigma(d,R,s2) + halfG_dsigma(-d, R, s2);
    }
    inline double DoG_dsigma(const double &d, const double &R, const double &sigma, const double &alpha)
    {
        return alpha*G_dsigma(d,R, pow(alpha*sigma, 2)) - G_dsigma(d, R, sigma*sigma);
    }
    inline double G_dsigma(const double &R, const double &sigma)
    {
        return -pow(R,3)/pow(sigma,4)*sqrt(2/M_PI)*exp(-pow(R,2)/2/pow(sigma,2));
    }
    inline double halfG_dsigma_dR(const double &d, const double &R, const double &s2)
    {
        return -R*(pow(R+d, 2)-s2)*exp(-pow(R+d, 2)/(2*s2))/sqrt(2*M_PI)/d/(s2*s2);
    }
    inline double G_dsigma_dR(const double &d, const double &R, const double &s2)
    {
        return halfG_dsigma_dR(d,R,s2) + halfG_dsigma_dR(-d, R, s2);
    }
    inline double DoG_dsigma_dR(const double &d, const double &R, const double &sigma, const double &alpha)
    {
        return alpha*G_dsigma_dR(d,R, pow(alpha*sigma, 2)) - G_dsigma_dR(d, R, sigma*sigma);
    }
    inline double G_dsigma_dR(const double &R, const double &sigma)
    {
        return pow(R,2)*(pow(R,2)-3*pow(sigma,2))/pow(sigma,6)*sqrt(2/M_PI)*exp(-pow(R,2)/2/pow(sigma,2));
    }
    """

def global_rescale_weave(sigma0, bonds, dists, R0=None, n=3):
    """Takes into account the overlapping of the blurred spot of neighbouring particles to compute the radii of all particles. Suppose all particles equally bright.
    
    parameters
    ----------
    sigma0 :  array((N))
        The radii output by MultiscaleTracker via scale2radius.
    bonds : array((M,2), int)
        First output of particles.get_bonds
    dists : array((M), float)
        Second output of particles.get_bonds
    R0 : array((N))
        Previous iteration's radii
    n : int
        Same as in MultiscaleTracker
    """
    assert len(bonds)==len(dists) 
    alpha = 2**(1.0/n)
    if R0==None:
        R0 = sigma2radius(sigma0, n=float(n))
    v0 = np.zeros([len(sigma0)])
    tr = np.zeros([len(sigma0)])
    jacob = np.zeros([len(bonds),2])
    code = """
    //using blitz++ expression
    v0 = -alpha*pow(R0,3)/pow(alpha*sigma0,4)*sqrt(2/M_PI)*exp(-pow(R0, 2)/2/pow(alpha*sigma0,2)) + pow(R0,3)/pow(sigma0,4)*sqrt(2/M_PI)*exp(-pow(R0, 2)/2/pow(sigma0,2));
    tr = alpha*pow(R0,2)*(pow(R0,2)-3*pow(alpha*sigma0,2))/pow(alpha*sigma0,6)*sqrt(2/M_PI)*exp(-pow(R0,2)/2/pow(alpha*sigma0,2)) - pow(R0,2)*(pow(R0,2)-3*pow(sigma0,2))/pow(sigma0,6)*sqrt(2/M_PI)*exp(-pow(R0,2)/2/pow(sigma0,2));
    #pragma omp parallel for
    for(int b=0; b<Nbonds[0]; ++b)
    {
        int i = bonds(b,0), j = bonds(b,1);
        jacob(b,0) = DoG_dsigma_dR(dists(b), R0(j), sigma0(i), alpha);
        jacob(b,1) = DoG_dsigma_dR(dists(b), R0(i), sigma0(j), alpha);
        v0(i) += DoG_dsigma(dists(b), R0(j), sigma0(i), alpha);
        v0(j) += DoG_dsigma(dists(b), R0(i), sigma0(j), alpha);
    }
    """
    weave.inline(
        code,['bonds', 'dists', 'sigma0', 'R0', 'alpha', 'v0', 'jacob', 'tr'],
        type_converters =converters.blitz,
        support_code = support_functions,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    jacob0 = sparse.lil_matrix(tuple([len(sigma0)]*2))
    for (i,j), (a,b) in zip(bonds, jacob):
        jacob0[i,j] = a
        jacob0[j,i] = b
    for i,a in enumerate(tr):
        jacob0[i,i] = a
    return R0 + spsolve(jacob0.tocsc(), -v0)
    
def global_rescale_intensity(sigma0, bonds, dists, intensities, R0=None, n=3):
    """Takes into account the overlapping of the blurred spot of neighbouring particles to compute the radii of all particles. The brightness of the particles is taken into account.
    
    parameters
    ----------
    sigma0 :  array((N))
        The radii output by MultiscaleTracker via scale2radius.
    bonds : array((M,2), int)
        First output of particles.get_bonds
    dists : array((M), float)
        Second output of particles.get_bonds
    intensities : array((N), float)
        Output of solve_intensities
    R0 : array((N), float)
        Previous iteration's radii
    n : int
        Same as in MultiscaleTracker
    """
    assert len(bonds)==len(dists) 
    alpha = 2**(1.0/n)
    if R0==None:
        R0 = sigma2radius(sigma0, n=float(n))
    v0 = np.zeros([len(sigma0)])
    tr = np.zeros([len(sigma0)])
    jacob = np.zeros([len(bonds),2])
    code = """
    //using blitz++ expression
    //tr = intensities * (alpha*pow(R0,2)*(pow(R0,2)-3*pow(alpha*sigma0,2))/pow(alpha*sigma0,6)*sqrt(2/M_PI)*exp(-pow(R0,2)/2/pow(alpha*sigma0,2)) - pow(R0,2)*(pow(R0,2)-3*pow(sigma0,2))/pow(sigma0,6)*sqrt(2/M_PI)*exp(-pow(R0,2)/2/pow(sigma0,2)));
    #pragma omp parallel for
    for(int i=0; i<Nv0[0]; ++i)
    {
        v0(i) = intensities(i) * (alpha*G_dsigma(R0(i), alpha*sigma0(i)) - G_dsigma(R0(i), sigma0(i)));
        tr(i) = intensities(i) * (alpha*G_dsigma_dR(R0(i), alpha*sigma0(i)) - G_dsigma_dR(R0(i), sigma0(i)));
    }
    #pragma omp parallel for
    for(int b=0; b<Nbonds[0]; ++b)
    {
        int i = bonds(b,0), j = bonds(b,1);
        const double d = std::max(dists(b), R0(i)+R0(j));
        jacob(b,0) = intensities(j) * DoG_dsigma_dR(d, R0(j), sigma0(i), alpha);
        jacob(b,1) = intensities(i) * DoG_dsigma_dR(d, R0(i), sigma0(j), alpha);
        v0(i) += intensities(j) * DoG_dsigma(d, R0(j), sigma0(i), alpha);
        v0(j) += intensities(i) * DoG_dsigma(d, R0(i), sigma0(j), alpha);
    }
    """
    weave.inline(
        code,['bonds', 'dists', 'sigma0', 'intensities', 'R0', 'alpha', 'v0', 'jacob', 'tr'],
        type_converters =converters.blitz,
        support_code = support_functions,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    jacob0 = sparse.lil_matrix(tuple([len(sigma0)]*2))
    for (i,j), (a,b) in zip(bonds, jacob):
        jacob0[i,j] = a
        jacob0[j,i] = b
    for i,a in enumerate(tr):
        jacob0[i,i] = a
    return R0 + spsolve(jacob0.tocsc(), -v0)
    
def solve_intensities(sigma0, bonds, dists, intensities, R0=None, n=3):
    """Takes into account the overlapping of the blurred spot of neighbouring particles to compute the brightness of all particles.
    
    parameters
    ----------
    sigma0 :  array((N))
        The radii output by MultiscaleTracker via scale2radius.
    bonds : array((M,2), int)
        First output of particles.get_bonds
    dists : array((M), float)
        Second output of particles.get_bonds
    intensities : array((N), float)
        Value of the Difference of Gaussian at the place and scale of each particle, e.g. the last column of the output of MultiscaleTracker
    R0 : array((N), float)
        Previous iteration's radii
    n : int
        Same as in MultiscaleTracker
    """
    assert len(bonds)==len(dists) 
    alpha = 2**(1.0/n)
    if R0==None:
        R0 = sigma2radius(sigma0, n=float(n))
    tr = np.zeros([len(sigma0)])
    ofd = np.zeros([len(bonds),2])
    code = """
    #pragma omp parallel for
    for(int i=0; i<Ntr[0]; ++i)
        tr(i) = G(R0(i), alpha*sigma0(i)) - G(R0(i), sigma0(i));
    #pragma omp parallel for
    for(int b=0; b<Nbonds[0]; ++b)
    {
        int i = bonds(b,0), j = bonds(b,1);
        const double d = dists(b);
        ofd(b,0) = DoG(d, R0(j), sigma0(i), alpha);
        ofd(b,1) = DoG(d, R0(i), sigma0(j), alpha);
    }
    """
    weave.inline(
        code,['bonds', 'dists', 'sigma0', 'R0', 'alpha', 'ofd', 'tr'],
        type_converters =converters.blitz,
        support_code = support_functions,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    mat = sparse.lil_matrix(tuple([len(sigma0)]*2))
    for (i,j), (a,b) in zip(bonds, ofd):
        mat[i,j] = a
        mat[j,i] = b
    for i,a in enumerate(tr):
        mat[i,i] = a
    return spsolve(mat.tocsc(), intensities)
    
    
halfG_dsigma = lambda d, R, sigma : (R**2+d*R+sigma**2)*np.exp(-(R+d)**2/(2*sigma**2))/np.sqrt(2*np.pi)/d/sigma**2

def G_dsigma(d, R, sigma):
    if d==0 or (R+d**2==R):
        return -(R**3)/sigma**4*np.sqrt(2/np.pi)*np.exp(-R**2/2/sigma**2)
    else:
        return halfG_dsigma(d,R,sigma) + halfG_dsigma(-d, R, sigma)

DoG_dsigma = lambda d, R, sigma, alpha : alpha*G_dsigma(d,R,alpha*sigma) - G_dsigma(d, R, sigma)

halfG_dsigma_dR = lambda d, R, sigma : -R*((R+d)**2-sigma**2)*np.exp(-(R+d)**2/(2*sigma**2))/np.sqrt(2*np.pi)/d/sigma**4
    
def G_dsigma_dR(d, R, sigma):
    if d==0 or (R+d**2==R):
        return R**2*(R**2-3*sigma**2)/sigma**6*np.sqrt(2/np.pi)*np.exp(-R**2/2/sigma**2)
    else:                                                                
        return halfG_dsigma_dR(d,R,sigma) + halfG_dsigma_dR(-d, R, sigma)

DoG_dsigma_dR = lambda d, R, sigma, alpha : alpha*G_dsigma_dR(d,R,alpha*sigma) - G_dsigma_dR(d, R, sigma)

def global_rescale(coords, sigma0, R0=None, bonds=None, n=3):
    alpha = 2**(1.0/n)
    if R0==None:
        R0 = sigma2radius(sigma0, n=float(n))
        v0 = np.zeros([len(coords)])
    else:
        v0 = DoG_dsigma(0, R0, sigma0, alpha)
    jacob0 = sparse.lil_matrix(tuple([len(coords)]*2))
    for i,(r,s) in enumerate(zip(R0, sigma0)):
        jacob0[i,i] = DoG_dsigma_dR(0, r, s, alpha)
    if bonds==None:
        bonds=[]
        for i,p in enumerate(coords):
            for j,dsq in enumerate(np.sum(numexpr.evaluate(
                    '(q-p)**2',
                    {'p':p, 'q':coords[i+1:]}
                    ), axis=-1)):
                if dsq<4*(R0[i]+R0[i+1+j])**2:
                    bonds.append([i, i+1+j, np.sqrt(dsq)])
    if len(bonds[0])==2:
        bonds = [[i, j, np.sqrt(np.sum((coords[i]-coords[j])**2))] for i,j in bonds]
    for i, j, d in bonds:
        jacob0[i, j] = DoG_dsigma_dR(d, R0[j], sigma0[i], alpha)
        jacob0[j, i] = DoG_dsigma_dR(d, R0[i], sigma0[j], alpha)
        v0[i] += DoG_dsigma(d, R0[j], sigma0[i], alpha)
        v0[j] += DoG_dsigma(d, R0[i], sigma0[j], alpha)
    return R0 + spsolve(jacob0.tocsc(), -v0)
    
class CrockerGrierFinder:
    """A single scale blob finder using Crocker & Grier algorithm"""
    def __init__(self, shape=(256,256), dtype=np.float32):
        """Allocate memory once"""
        self.blurred = np.empty(shape, dtype)
        self.background = np.empty(shape, dtype)
        self.dilated = np.empty_like(self.blurred)
        self.binary = np.empty(self.blurred.shape, bool)
        
    def fill(self, image, k=1.6, uniform_size=None, background=None):
        """All the image processing when accepting a new image."""
        assert self.blurred.shape == image.shape, """Wrong image size:
%s instead of %s"""%(image.shape, self.blurred.shape)
        #fill the first layer by the input
        self.blurred[:] = image
        #Gaussian filter
        gaussian_filter(self.blurred, k, output=self.blurred)
        #background removal
        if background is None:
            if uniform_size is None:
                uniform_size = int(10*k)
            if uniform_size>0:
                uniform_filter(self.blurred, uniform_size, output=self.background)
                self.blurred -= self.background
        else:
            self.blurred -= background
        #Dilation
        grey_dilation(self.blurred, [3]*self.blurred.ndim, output=self.dilated)
        
    def initialize_binary(self, maxedge=1.1, threshold=None):
        if threshold is None:
            self.binary[:] = self.blurred == self.dilated
        else:
            self.binary = numexpr.evaluate(
                '(b==d) & (b>thr)',
                {'b':self.blurred, 'd':self.dilated, 'thr':threshold}
                )
        #eliminate particles on the edges of image
        for a in range(self.binary.ndim):
            self.binary[tuple([slice(None)]*(self.binary.ndim-1-a)+[slice(0,2)])]=False
            self.binary[tuple([slice(None)]*(self.binary.ndim-1-a)+[slice(-2, None)])]=False
        #eliminate blobs that are edges
        if self.blurred.ndim==2 and maxedge>0 :
            for p in np.transpose(np.where(self.binary)):
                #xy neighbourhood
                ngb = self.blurred[tuple([slice(u-1, u+2) for u in p])]
                #compute the XYhessian matrix coefficients
                hess = [
                    ngb[0, 1] - 2*ngb[1, 1] + ngb[-1,1],
                    ngb[1, 0] - 2*ngb[1, 1] + ngb[1,-1],
                    ngb[0,0] + ngb[-1,-1] - ngb[0,-1] - ngb[-1,0]
                    ]
                #determinant of the Hessian, for the coefficient see
                #H Bay, a Ess, T Tuytelaars, and L Vangool,
                #Computer Vision and Image Understanding 110, 346-359 (2008)
                detH = hess[0]*hess[1] - hess[2]**2
                ratio = (hess[0]+hess[1])**2/(4.0*hess[0]*hess[1])
                if detH<0 or ratio>maxedge:
                    self.binary[tuple(p.tolist())] = False
                    
    def no_subpix(self):
        """extracts centers positions and values from binary without subpixel resolution"""
        nb_centers = self.binary.sum()
        if nb_centers==0 or self.binary.min():
            return np.zeros([0, self.blurred.ndim+1])
        #original positions of the centers
        c0 = np.transpose(np.where(self.binary))
        vals = self.blurred[self.binary]
        return np.column_stack((vals, c0))
        
    def subpix(self):
        """Extract and refine to subpixel resolution the positions and size of the blobs"""
        nb_centers = self.binary.sum()
        if nb_centers==0 or self.binary.min():
            return np.zeros([0, self.blurred.ndim+1])
        centers = np.empty([nb_centers, self.blurred.ndim+1])
        #original positions of the centers
        c0 = np.transpose(np.where(self.binary))
        if self.binary.ndim==2:
            im = self.blurred
            code = """
            #pragma omp parallel for
            for(int p=0; p<Nc0[0]; ++p)
            {
                const int i = c0(p,0), j = c0(p,1);
                centers(p,0) = blitz::sum(im(blitz::Range(i-2,i+2), blitz::Range(j-2,j+2)))/25.;
                centers(p,1) = i - blitz::sum(im(blitz::Range(i-2,i+2), j) * coefprime) / blitz::sum(im(blitz::Range(i-2,i+2), j) * coefsec);
                centers(p,2) = j - blitz::sum(im(i, blitz::Range(j-2,j+2)) * coefprime) / blitz::sum(im(i, blitz::Range(j-2,j+2)) * coefsec);
            }
            """
            weave.inline(
                code,['c0', 'centers', 'coefprime', 'coefsec', 'im'],
                type_converters =converters.blitz,
                extra_compile_args =['-O3 -fopenmp -mtune=native'],
                extra_link_args=['-lgomp'],
                verbose=2)
            return centers
        for i, p in enumerate(c0):
            #neighbourhood
            ngb = self.blurred[tuple([slice(u-2, u+3) for u in p])]
            for dim in range(ngb.ndim):
                a = ngb[tuple([1]*(ngb.ndim-1-dim)+[slice(None)]+[1]*dim)]
                centers[i,dim+1] = p[dim] - np.dot(coefprime,a)/np.dot(coefsec,a)
            centers[i,0] = ngb.mean()
        return centers
        
    def __call__(self, image, k=1.6, maxedge=1.1, threshold=None, uniform_size=None, background=None):
        """Locate bright blobs in an image with subpixel resolution.
Returns an array of (x, y, intensity)"""
        self.fill(image, k, uniform_size, background)
        self.initialize_binary(maxedge, threshold)
        centers = self.subpix()[:,::-1]
        return centers

class OctaveBlobFinder:
    """Locator of bright blobs in an image of fixed shape. Works on a single octave."""
    def __init__(self, shape=(256,256), nbLayers=3, dtype=np.float32):
        """Allocate memory once"""
        self.layersG = np.empty([nbLayers+3]+list(shape), dtype)
        self.layers = np.empty([nbLayers+2]+list(shape), dtype)
        self.eroded = np.empty_like(self.layers)
        self.binary = np.empty(self.layers.shape, bool)
        self.time_fill = 0.0
        self.time_subpix = 0.0
        self.ncalls = 0
        self.noutputs = 0
        self.sizes = np.empty(nbLayers)

    def get_iterative_radii(self, k):
        nbLayers = len(self.layersG)-3
        #target blurring radii
        sigmas = k*2**(np.arange(nbLayers+3)/float(nbLayers))
        #corresponding blob sizes and iterative blurring radii
        return np.rint(sigmas*np.sqrt(2)).astype(int), np.sqrt(np.diff(sigmas**2))
        
    def fill(self, image, k=1.6):
        """All the image processing when accepting a new image."""
        t0 = time.clock()
        #total 220 ms
        assert self.layersG[0].shape == image.shape, """Wrong image size:
%s instead of %s"""%(image.shape, self.layersG[0].shape)
        #fill the first layer by the input (already blurred by k)
        self.layersG[0] = image
        self.sizes, sigmas_iter = self.get_iterative_radii(k)
        #Gaussian filters
        for l, layer in enumerate(self.layersG[:-1]):
            gaussian_filter(layer, sigmas_iter[l], output=self.layersG[l+1])
        #Difference of gaussians
        for l in range(len(self.layers)):
            self.layers[l] = self.layersG[l+1] - self.layersG[l]
        #Erosion 86.2 ms
        grey_erosion(self.layers, [3]*self.layers.ndim, output=self.eroded)
        #scale space minima, whose neighbourhood are all negative 10 ms
        self.time_fill += time.clock()-t0

    def initialize_binary(self, maxedge=1.1, first_layer=False, maxDoG=None):
        """Convert the DoG layers into the binary image, True at center position.
        
        Centers are local minima of the DoG with negative value 
            if maxDog is None, the DoG value should be further from 0 than machine precision. 
            else, the DoG value must be lower than maxDog
        Centers at the edge of the image are excluded. 
        On 2D images, if maxedge is positive, elongated blobs are excluded if the ratio of the eignevalues of the Hessian matrix is larger than maxedge.
        Optionally, the local spatial minima in the first DoG layer can be considered as centers.
        """
        if maxDoG is None:
            maxDoG = 0
        #local minima in the DoG on both space and scale are obtained from erosion
        self.binary = numexpr.evaluate(
            '(l==e) & (l<maxDoG) & (l**2+1.0>1.0)',
            {
                'l': self.layers,
                'e': self.eroded,
                'maxDoG': maxDoG
                }
            )
        #If the first DoG layer is taken into account, 
        #its centers are translated in the layer above.
        #Necessary for subpixel
        if first_layer:
            #if a local maximum in the Gaussian layer 1 is 
            #within 1px of a potential (but doomed) center in the DoG layer 0
            #add it to binary layer 1
            self.binary[0] = binary_dilation(
                self.binary[0], 
                np.ones([3]*(self.layers.ndim-1))
                )
            self.binary[0] &= grey_dilation(
                self.layersG[1], 
                [3]*(self.layers.ndim-1)
                ) == self.layersG[1]
            self.binary[1] |= self.binary[0]
        #centers in the first and last layers are discarded
        #do not remove these lines or you break subpixel
        self.binary[0] = False
        self.binary[-1] = False
        #eliminate particles on the edges of image
        for r, bi in zip(self.sizes[1:-1], self.binary[1:-1]):
            for a in range(bi.ndim):
                bi[tuple([slice(None)]*(bi.ndim-1-a)+[slice(0,r)])]=False
                bi[tuple([slice(None)]*(bi.ndim-1-a)+[slice(-r, None)])]=False
        #eliminate blobs that are edges
        if self.layers.ndim==3 and maxedge>0 :
            for r, bi, layer in zip(self.sizes[1:-1], self.binary[1:-1], self.layers[1:-1]):
                for p in np.transpose(np.where(bi)):
                    #xy neighbourhood
                    ngb = layer[tuple([slice(u-1, u+2) for u in p])]
                    #compute the XYhessian matrix coefficients
                    hess = [
                        ngb[0, 1] - 2*ngb[1, 1] + ngb[-1,1],
                        ngb[1, 0] - 2*ngb[1, 1] + ngb[1,-1],
                        ngb[0,0] + ngb[-1,-1] - ngb[0,-1] - ngb[-1,0]
                        ]
                    #determinant of the Hessian, for the coefficient see
                    #H Bay, a Ess, T Tuytelaars, and L Vangool,
                    #Computer Vision and Image Understanding 110, 346-359 (2008)
                    detH = hess[0]*hess[1] - hess[2]**2
                    ratio = (hess[0]+hess[1])**2/(4.0*hess[0]*hess[1])
                    if detH<0 or ratio>maxedge:
                        bi[tuple(p.tolist())] = False

    def no_subpix(self):
        """extracts centers positions and values from binary without subpixel resolution"""
        nb_centers = self.binary.sum()
        if nb_centers==0 or self.binary.min():
            return np.zeros([0, self.layers.ndim+1])
        #original positions of the centers
        c0 = np.transpose(np.where(self.binary))
        vals = self.layers[self.binary]
        return np.column_stack((vals, c0))

    def subpix(self, method=1):
        """Extract and refine to subpixel resolution the positions and size of the blobs"""
        nb_centers = self.binary.sum()
        if nb_centers==0 or self.binary.min():
            return np.zeros([0, self.layers.ndim+1])
        centers = np.empty([nb_centers, self.layers.ndim+1])
        #original positions of the centers
        c0 = np.transpose(np.where(self.binary))
        if method==0:
            for i, p in enumerate(c0):
                #neighbourhood
                ngb = self.layers[tuple([slice(u-1, u+2) for u in p])]
                grad = np.asarray([
                    (n[-1]-n[0])/2.0
                    for n in [
                        ngb[tuple(
                            slice(None) if a==u else 1
                            for u in range(ngb.ndim)
                            )]
                        for a in range(ngb.ndim)]
                    ])
                hess = np.empty([ngb.ndim]*2)
                for a in range(ngb.ndim):
                    n = ngb[tuple([1]*(ngb.ndim-1-a)+[slice(None)]+[1]*a)]
                    hess[a,a] = n[-1]+n[0]-2*n[1]
                for a in range(ngb.ndim-1):
                    for b in range(a+1,ngb.ndim):
                        n = ngb[tuple(
                            slice(None) if u==a or u==b else 1
                            for u in range(ngb.ndim)
                            )]
                        hess[a,b] = (n[0,0] + n[-1,-1] - n[0,-1] - n[-1,0])/4.0
                        hess[b,a] = hess[a,b]
                dx = - np.dot(np.linalg.inv(hess), grad)
                centers[i,1:] = p + dx
                centers[i,0] = ngb[tuple([1]*ngb.ndim)]+0.5*np.dot(dx,grad)
        else:
            for i, p in enumerate(c0):
                #neighbourhood, three pixels in the scale axis,
                #but according to scale in space
                r = self.sizes[p[0]]
                rv = [1]+[r]*(self.layers.ndim-1)
                ngb = np.copy(self.layers[tuple(
                    [slice(p[0]-1,p[0]+2)]+[
                        slice(u-r, u+r+1) for u in p[1:]
                        ]
                    )])
                #label only the negative pixels
                labels = measurements.label(ngb<0)[0]
                lab = labels[tuple(rv)]
                #value
                centers[i,0] = measurements.mean(ngb, labels, [lab])
                #pedestal removal
                ped = measurements.maximum(ngb, labels, [lab])
                if ped!=self.layers[tuple(p.tolist())]: #except if only one pixel or uniform value
                    ngb -= ped
                #center of mass
                centers[i,1:] = (np.asanyarray(measurements.center_of_mass(
                    ngb, labels, [lab]
                    ))-rv)+p
                #the subscale resolution is calculated using only 3 pixels
                n = ngb[tuple(
                    [slice(None)]+[r]*(self.layers.ndim-1)
                    )].ravel()
                denom = (n[2] - 2 * n[1] + n[0])
                if (abs(denom)+1.0)**2 > 1.0:
                    centers[i,1] = p[0] - (n[2] - n[0]) / 2.0 / denom
                else: centers[i,1] = p[0]
        return centers
        
    def __call__(self, image, k=1.6, maxedge=1.1, first_layer=False, maxDoG=None):
        """Locate bright blobs in an image with subpixel resolution.
Returns an array of (x, y, r, -intensity in scale space)"""
        self.ncalls += 1
        self.fill(image, k)
        self.initialize_binary(maxedge, first_layer, maxDoG)
        t0 = time.clock()
        centers = self.subpix()[:,::-1]
        self.time_subpix += time.clock() - t0
        #convert scale to size
        n = (len(self.layers)-2)
        centers[:,-2] = scale2radius(centers[:,-2], k, n, self.layers.ndim-1)
        self.noutputs += len(centers)
        return centers
        
        
class MultiscaleBlobFinder:
    """Locator of bright blobs in an image of fixed shape. Works on more than one octave, starting at octave -1."""
    def __init__(self, shape=(256,256), nbLayers=3, nbOctaves=3, dtype=np.float32, Octave0=True):
        """Allocate memory for each octave"""
        shapes = np.vstack([np.ceil([s*2.0**(Octave0-o) for s in shape]) for o in range(nbOctaves)])
        self.preblurred = np.empty(shapes[0], dtype)
        self.octaves = [
            OctaveBlobFinder(s, nbLayers, dtype)
            for s in shapes if s.min()>8
            ] #shortens the list of octaves if no blob can be detected in that small window
        if not Octave0:
            self.octaves.insert(0, OctaveBlobFinder([0]*len(shape), nbLayers, dtype))
        self.Octave0 = Octave0
        self.time = 0.0
        self.ncalls = 0
        
    def __call__(self, image, k=1.6, Octave0=True, removeOverlap=True, maxedge=1.1, deconvKernel=None, first_layer=False, maxDoG=None):
        """Locate blobs in each octave and regroup the results"""
        if not self.Octave0:
            Octave0 = False
        self.ncalls += 1
        t0 = time.clock()
        if len(self.octaves)==0:
            return np.zeros([0, image.ndim+2])
        #upscale the image for octave -1
        #halfbl = gaussian_filter(np.array(im, , k/2.0)
        if Octave0:
            im2 = np.copy(image)
            for a in range(image.ndim):
                im2 = np.repeat(im2, 2, a)
            #preblur octave -1
            gaussian_filter(im2, k, output=self.preblurred)
            #locate blobs in octave -1
            centers = [self.octaves[0](self.preblurred, k, maxedge, maxDoG=maxDoG)]
        else:
            centers = []
        if len(self.octaves)>1:
            if deconvKernel is not None:
                assert len(deconvKernel) == image.shape[0]/2+1
                assert image.ndim == 3
                #deconvolve the Z direction by a precalculated kernel
                #To avoid noise amplification, the blurred image is deconvolved, not the raw one
                deconv = deconvolve(gaussian_filter(image.astype(float), k), deconvKernel)
                #remove negative values
                centers += [self.octaves[1](np.maximum(0, deconv), maxedge=maxedge, first_layer=first_layer, maxDoG=maxDoG)]
            else:
                centers += [self.octaves[1](gaussian_filter(image, k), maxedge=maxedge, first_layer=first_layer, maxDoG=maxDoG)]
        #subsample the -3 layerG of the previous octave
        #which is two times more blurred that layer 0
        #and use it as the base of new octave
        for o, oc in enumerate(self.octaves[2:]):
            centers += [oc(
                self.octaves[o+1].layersG[-3][
                    tuple([slice(None, None, 2)]*image.ndim)],
                k, maxedge, maxDoG=maxDoG
                )]
        #merge the results and scale the coordinates and sizes
        centers = np.vstack([
            c * ([2**(o-Octave0)]*(1+image.ndim)+[1])
            for o, c in enumerate(centers)
            ])
        if len(centers)<2:
            return centers
        if not removeOverlap:
            return centers
        #remove overlaping objects (keep the most intense)
        #scales in dim*N^2, thus expensive if many centers and in high dimensions
        #Using a spatial index may be faster
        out = []
        for i in centers[np.argsort(centers[:,-1])]:
            for j in out:
                if np.sum((i[:-2]-j[:-2])**2) < (i[-2]+j[-2])**2:
                    break
            else:
                out.append(i)
        self.time += time.clock() - t0
	return np.vstack(out)

def treatFrame(serie, t, file_pattern, finder=None ):
    if finder is None:
        finder = MultiscaleBlobFinder(serie.get2DShape())
    stack = serie.getFrame(T=t)
    for z, im in enumerate(stack):
        centers = finder(im)
        np.savetxt(
            file_pattern%(t, z, 'dat'),
            np.vstack((
                    [1, len(centers), 1],
                    serie.get2DShape()+[1.6*np.sqrt(2)*2**(7.0/3)],
                    centers[:,:-1]
                    )), fmt='%g'
            )
        np.savetxt(file_pattern%(t, z, 'intensity'), centers[:,-1], fmt='%g')
    pro = subprocess.Popen([
        '/home/mathieu/test/bin/linker',
        file_pattern%(t,0,'dat'),
        '_z', '1', '%g'%serie.getZXratio(), '%d'%len(stack)
        ], stdout=subprocess.PIPE)
    trajfile = pro.communicate()[0].split()[-1]
    clusters = load_clusters(trajfile)
    particles = clusters2particles(clusters)
    np.save(os.path.splitext(trajfile)[0], particles)

def localize2D3D(serie, file_pattern, cleanup=True):
    finder = MultiscaleBlobFinder(serie.get2DShape())
    for t in range(serie.getNbFrames()):
        treatFrame(serie, t, file_pattern, finder)
        if cleanup:
            for z in range(serie.getFrameShape()[-1]):
                os.remove(file_pattern%(t, z, 'dat'))
                os.remove(file_pattern%(t, z, 'intensity'))
              
        
  
                
if __name__ == '__main__':
    import os.path
    #tests
    class TestCrockerGrier(unittest.TestCase):

        def test_init(self):
            im = np.load(os.path.join(os.path.dirname(__file__), 'dillute_raw.npy'))
            finder = CrockerGrierFinder(im.shape)
            centers = finder(im, threshold=0.1)
            self.assertEqual(len(centers), 5)
    unittest.main()
