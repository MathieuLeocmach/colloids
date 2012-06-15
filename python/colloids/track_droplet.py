import numpy as np
from colloids import track, particles
from scipy.misc import imsave, imshow, imread
import sys, itertools, re
from os.path import split,join,isfile, splitext
import Image

def readTIFF16(path):
    """Read 16bits TIFF"""
    im = Image.open(path)
    out = np.fromstring(
        im.tostring(), 
        np.uint8
        ).reshape(tuple(list(im.size)+[2]))
    return (np.array(out[:,:,0], np.uint16)<<8)+out[:,:,1]

#read input file name
if len(sys.argv) >1:
    filename = sys.argv[1]
else:
    while(True):
        filename = raw_input("filename --> ")
        if isfile(filename):
            break
        else:
            print '"%s" is not an existing file' % filename

#how is made the file series name template
#the last groupe of digits afet a _ is supposed to be the time step
m = re.match('(.*)_([0-9]*)(.*)', filename)
assert m is not None, '"%s" do not match pattern' % filename
pattern = m.group(1)+'_%'+('0%dd'%len(m.group(2)))+(''.join(m.groups()[2:]))
path, name = split(m.group(1))
outpattern = join(path,'treated_'+name)+'_%'+('0%dd'%len(m.group(2)))+(''.join(m.groups()[2:]))

#read the first image
im = imread(filename)
if im.ndim ==0:
    im = readTIFF16(filename)
cent = np.zeros_like(im)
finder = track.MultiscaleBlobFinder(im.shape, nbOctaves=100)



for t in itertools.count():
    if not isfile(pattern%t):
        break
    im = imread(pattern%t)
    if im.ndim ==0:
        im = readTIFF16(pattern%t)
    #localisation
    centers = finder(im, Octave0=False, removeOverlap=True)
    centers = centers[centers[:,-1]<-1]
    scales = track.radius2scale(centers[:,-2], dim=2)
    centers = centers[scales>0]
    scales = scales[scales>0]
    #output result image
    cent.fill(0)
    for (ym,xm), (yM, xM) in zip(np.maximum(0, centers[:,:2]-centers[:,2][:,None]), np.minimum(centers[:,:2]+centers[:,2][:,None]+1, im.shape[::-1])):
        cent[xm:xM,ym:yM] += 50/256.*im.max()
    imsave(outpattern%t, np.dstack((im, cent, cent)))
    #look for Gaussian intensities
    for i, s in enumerate(scales):
        o = int(s/3)
        centers[i,-1] = finder.octaves[o].layersG[s-3*o, centers[i,1]/2**o, centers[i,0]/2**o]
    #Trajectory linking
    if t==0:
        linker = particles.Linker(len(centers))
        pos1 = centers
    else:
        pos0 = pos1
        pos1 = centers
        pairs, distances = particles.get_links_size(pos0[:,:-2], pos0[:,-2], pos1[:,:-2], pos1[:,-2], maxdist=2)
        linker.addFrame(len(pos1), pairs, distances)
    #output
    np.savetxt(
        splitext(pattern%t)[0]+'.csv', 
        np.column_stack((linker.pos2tr[-1], centers)),
        fmt='%g'
        )
        
        
    
