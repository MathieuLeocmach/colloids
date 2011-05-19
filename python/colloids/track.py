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
from scipy.ndimage.filters import gaussian_filter, gaussian_filter1d, sobel
from scipy.ndimage.morphology import grey_erosion, grey_dilation, binary_dilation
from scipy.ndimage import measurements


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
        if len(cl)<5:
            continue
        #get the blobs for each signal, remove signals without blob
        blobs = filter(len, [finders[len(cl)-5](u, k) for u in [
            -np.sqrt(sobel(cl[:,0], axis=0)**2 + sobel(cl[:,1], axis=0)**2),
            cl[:,3], cl[:,4]
            ]])
        if len(blobs)==0: #no blob in any signal
            continue
        #sort the blob of each signal by intensity (negative)
        blobs = np.vstack([bs[np.argsort(bs[:,-1])] for bs in blobs])
        #Remove overlapping centers than may appear at different scales
        #in the different signals.
        #The most intense blob in the gradient of position is the best,
        #then, the second intense in the gradient of position, etc.
        #then, the blobs in apparent radius
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

class OctaveBlobFinder:
    """Locator of bright blobs in an image of fixed shape. Works on a single octave."""
    def __init__(self, shape=(256,256), nbLayers=3):
        """Allocate memory once"""
        self.layersG = np.empty([nbLayers+3]+list(shape), float)
        self.layers = np.empty([nbLayers+2]+list(shape), float)
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
        self.layers[:] = np.diff(self.layersG, axis=0) #5.99 ms
        #Erosion 86.2 ms
        grey_erosion(self.layers, [3]*self.layers.ndim, output=self.eroded)
        #scale space minima, whose neighbourhood are all negative 10 ms
        self.time_fill += time.clock()-t0

    def initialize_binary(self):
        self.binary = self.layers==self.eroded
        self.binary[0] = False
        self.binary[-1] = False
        #eliminate particles on the edges of image
        for r, bi in zip(self.sizes[1:-1], self.binary[1:-1]):
            for a in range(bi.ndim):
                bi[tuple([slice(None)]*(bi.ndim-1-a)+[slice(0,r)])]=False
                bi[tuple([slice(None)]*(bi.ndim-1-a)+[slice(-r, None)])]=False
        #eliminate blobs that are edges
        if self.layers.ndim==3:
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
                    if detH<0 or ratio>1.1:
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
                    [slice(None)]+[slice(r-1,r+2)]*(self.layers.ndim-1)
                    )]
                n -= n.max()
                ds = measurements.center_of_mass(n)[0]-1
                if np.abs(ds)>0.33:
                    if ds<0:
                        ds = (n[1]-n[0]).sum()/(n[:-1]-n[:-1].max()).sum()
                    else:
                        ds = (n[2]-n[1]).sum()/(n[1:]-n[1:].max()).sum()
                #the scale axis is logarythmic
                centers[i,1] = np.exp(ds)*p[0]
        return centers
        
    def __call__(self, image, k=1.6):
        """Locate bright blobs in an image with subpixel resolution.
Returns an array of (x, y, r, -intensity in scale space)"""
        self.ncalls += 1
        self.fill(image, k)
        self.initialize_binary()
        t0 = time.clock()
        centers = self.subpix()[:,::-1]
        self.time_subpix += time.clock() - t0
        #convert scale to size
        centers[:,-2] = 0.5*(1+k)*(k*np.sqrt(2)*2**(centers[:,-2]/(len(self.layers)-2))-1.0/(len(self.layers)-2))
        self.noutputs += len(centers)
        return centers
        
        
class MultiscaleBlobFinder:
    """Locator of bright blobs in an image of fixed shape. Works on more than one octave, starting at octave -1."""
    def __init__(self, shape=(256,256), nbLayers=3, nbOctaves=3):
        """Allocate memory for each octave"""
        shapes = np.vstack([np.ceil([s*2.0**(1-o) for s in shape]) for o in range(nbOctaves)])
        self.preblurred = np.empty(shapes[0])
        self.octaves = [
            OctaveBlobFinder(s, nbLayers)
            for s in shapes if s.min()>8
            ] #shortens the list of octaves if no blob can be detected in that small window
        self.time = 0.0
        self.ncalls = 0
        
    def __call__(self, image, k=1.6):
        """Locate blobs in each octave and regroup the results"""
        self.ncalls += 1
        t0 = time.clock()
        if len(self.octaves)==0:
            return np.zeros([0, image.ndim+2])
        #upscale the image for octave -1
        im2 = np.copy(image)
        for a in range(image.ndim):
            im2 = np.repeat(im2, 2, a)
        #preblur octave -1
        gaussian_filter(im2, k, output=self.preblurred)
        #locate blobs in octave -1
        centers = [self.octaves[0](self.preblurred, k)]
        #subsample the -3 layerG of the previous octave
        #which is two times more blurred that layer 0
        #and use it as the base of new octave
        for o, oc in enumerate(self.octaves[1:]):
            centers += [oc(
                self.octaves[o].layersG[-3][
                    tuple([slice(None, None, 2)]*image.ndim)],
                k,
                )]
        #merge the results and scale the coordinates and sizes
        centers = np.vstack([
            c * ([2**(o-1)]*(1+image.ndim)+[1])
            for o, c in enumerate(centers)
            ])
        self.time += time.clock() - t0
	return centers

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
                    [256, 256, 1.6*np.sqrt(2)*2**(7.0/3)],
                    centers[:,:-1]
                    )), fmt='%g'
            )
        np.savetxt(file_pattern%(t, z, 'intensity'), centers[:,-1], fmt='%g')
    pro = subprocess.Popen([
        '/home/mathieu/test/bin/linker',
        file_pattern%(t,0,'dat'),
        '_z', '5', '%g'%serie.getZXratio(), '%d'%len(stack)
        ], stdout=subprocess.PIPE)
    trajfile = pro.communicate()[0].split()[-1]
    trajfile = os.path.join(os.path.split(file_pattern)[0], trajfile)
    clusters = load_clusters(trajfile)
    particles = clusters2particles(clusters)
    np.save(os.path.splitext(trajfile)[0], particles)

def localize2D3D(serie, file_pattern, cleanup=True):
    finder = MultiscaleBlobFinder(serie.get2DShape())
    for t in range(serie.getNbFrames()):
        treatFrame(serie, t, file_pattern, finder)
        if cleanup:
            for z in range(len(stack)):
                os.remove(file_pattern%(t, z, 'dat'))
                os.remove(file_pattern%(t, z, 'intensity'))
