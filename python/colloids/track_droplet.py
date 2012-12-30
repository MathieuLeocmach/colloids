import numpy as np
from colloids import track, particles, phaseCorrelation
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
    return (np.array(out[:,:,0], np.uint16) << 8) + out[:,:,1]

if __name__ == "__main__":
    #read input file name
    if len(sys.argv) >1:
        filename = sys.argv[1]
    else:
        while(True):
            filename = raw_input("filename --> ")
            try:
                if isfile(filename) or isfile(filename%0):
                    break
                else:
                    print '"%s" is not an existing file' % filename
            except TypeError:
                print '"%s" is not a valid pattern' % filename

    if isfile(filename):
        #how is made the file series name template
        #the last groupe of digits afet a _ is supposed to be the time step
        m = re.match('(.*)_([0-9]*)(.*)', filename)
        assert m is not None, '"%s" do not match pattern' % filename
        pattern = m.group(1)+'_%'+('0%dd'%len(m.group(2)))+(''.join(m.groups()[2:]))
        path, name = split(m.group(1))
        outpattern = join(path,'treated_'+name)+'_%'+('0%dd'%len(m.group(2)))+(''.join(m.groups()[2:]))
    else:
        pattern = filename
        path, name = split(filename)
        outpattern = join(path,'treated_'+name)
        filename = pattern%0

    #read the first image
    im = imread(filename)
    if im.ndim ==0:
        im = readTIFF16(filename)
    cent = np.zeros_like(im)
    #estimate the best blurring radius
    finder = track.MultiscaleBlobFinder(im.shape, nbOctaves=100)
    centers = finder(im, Octave0=False, removeOverlap=False)
    sigma = 1.6*2**((1+np.argmax([oc.binary[1:-1].sum(-1).sum(-1) for oc in finder.octaves[1:]]))/3.)
    #create a monoscale Crocker & Grier finder
    finder = track.CrockerGrierFinder(im.shape)
    #create drift file
    fdrift = open(splitext(pattern%0)+'.shifts', 'w')

    for t in itertools.count():
        if not isfile(pattern%t):
            break
        im = imread(pattern%t)
        if im.ndim ==0:
            im = readTIFF16(pattern%t)
        #localisation
        centers = finder(im, sigma)
        #output result image
        cent.fill(0)
        for (ym,xm), (yM, xM) in zip(np.maximum(0, centers[:,:2]-sigma), np.minimum(centers[:,:2]+sigma+1, im.shape[::-1])):
            cent[xm:xM,ym:yM] += 50/256.*im.max()
        imsave(outpattern%t, np.dstack((im, cent, cent)))
        #Trajectory linking
        if t==0:
            linker = particles.Linker(len(centers))
            pos1 = centers
            sp1 = np.fft.fftn(
                im*np.hanning(im.shape[0])[:,None]*np.hanning(im.shape[1])
                )
        else:
            pos0 = pos1
            pos1 = centers
            #look for overall shift between successive pictures
            sp0 = sp1
            sp1 = np.fft.fftn(
                im*np.hanning(im.shape[0])[:,None]*np.hanning(im.shape[1])
                )
            R = sp0 * np.conjugate(sp1)
            shift = phaseCorrelation.getDispl(np.abs(np.fft.ifftn(R/np.abs(R))))
            fdrift.write('%d %d\n'%(shift[1], shift[0]))
            pairs, distances = particles.get_links(pos0[:,:2], np.ones(len(pos0)), pos1[:,:2]+shift[::-1], np.ones(len(pos1)), maxdist=6.*sigma)
            linker.addFrame(len(pos1), pairs, distances)
        #output
        np.savetxt(
            splitext(pattern%t)[0]+'.csv', 
            np.column_stack((linker.pos2tr[-1], centers)),
            fmt='%g'
            )
        
        
    
