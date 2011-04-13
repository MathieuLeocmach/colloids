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
import os.path, subprocess, shlex, string, re
from colloids import lif, vtk
from scipy.ndimage.filters import gaussian_filter, gaussian_filter1d, sobel
from scipy.ndimage.morphology import grey_erosion, grey_dilation
from scipy.ndimage import measurements
from rtree import Rtree


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

def multiple_centers(centers):
    """Remove overlapping centers (2D)"""
    #sort the centers by decreasing radius
    s = s[np.argsort(s[:,0])][::-1]
    idx = Rtree()
    for p,(r, y, x) in enumerate(s):
        bb = (x-r, y-r, x+r, y+r)
        inter = list(idx.intersection(bb))
        #don't add centers that overlap
        if len(inter)>0:
            qs = s[inter]
            if np.sum(np.sum((qs[:,1:]-s[p][1:])**2)<(qs[:,0]+r)**2)>0:
                continue
        idx.add(p, bb)
    return s[list(idx.intersection(idx.bounds))]

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

def cluster2radgrad(cluster, k=1.6):
    """Enhance the cusps between stacked particles in the radius(z) curve
by substracting the magnitude of the gradient of x(z) and y(z) ;
and adding the intensity"""
    grad = np.sqrt(np.sum(gaussian_filter1d(
        cluster[:,:2], k, axis=0, order=1
        )**2, axis=-1))
    smoothI = gaussian_filter1d(cluster[:,4], k)
    pos_maxratio = np.argmax(cluster[:,3]/smoothI)
    return cluster[:,3] - grad + smoothI/smoothI[pos_maxratio]

def filter_blobs1d(blobs, k=1.6, n=3):
    """Remove overlapping centers than may appear at different scales.
Keeps the center with the strongest (negative) signal."""
    if len(blobs)==0:
        return blobs
    #sort by intensity
    inb = blobs[np.argsort(blobs[:,-1])]
    out = []
    for i in inb:
            for j in out:
                    if (i[1]-j[1])**2 < (k*(2**(i[0]/n)+2**(j[0]/3)))**2/2:
                            break
            else:
                    out.append(i)
    return np.asarray(out)
    
def clusters2particles(clusters, k=1.6, n=3, noDuplicate=True):
    particles = []
    for cl in clusters:
        blobs = np.vstack((
            find_blob(cluster2radgrad(np.repeat(cl,2,0)), k, n)*[1, 0.5, 1],
            find_blob(cluster2radgrad(cl), k, n)+[n, 0, 0],
            (find_blob(cluster2radgrad(cl[::2]), k, n)*[1,2,1])+[2*n, 0, 0]
            ))
        if noDuplicate:
            blobs = filter_blobs1d(blobs)
        grad = gaussian_filter1d(cl, k/2, axis=0, order=1)
        for s, z, v in blobs:
            #smoothed = (gaussian_filter1d(cl, 1.6*2**(s/3-1), axis=0)-cl.mean(0))*np.sqrt(2)+cl.mean(0)
            #grad = gaussian_filter1d(cl, 1.6*2**(s/3-1), axis=0, order=1)
            zi = np.rint(z)
            if zi==len(cl):
                zi = len(cl)-1
            if zi<0:
                zi = 0
            dz = z-zi
            particles.append(cl[zi]+grad[zi]*dz)
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
