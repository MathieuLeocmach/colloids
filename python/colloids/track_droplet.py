import numpy as np
from colloids import track, particles, phaseCorrelation
from scipy.misc import imsave, imshow, imread
import sys, itertools, re
from os.path import split,join,isfile, splitext
import Image
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline
from colloids.colors import colorscale
import matplotlib, pylab as plt
from matplotlib.colors import ListedColormap

matplotlib.cm.register_cmap("bluetored", ListedColormap(colorscale(np.linspace(0,1,256,False)[None,:])[0]))

def readTIFF16(path, bigendian=True):
    """Read 16bits TIFF"""
    im = Image.open(path)
    out = np.fromstring(
        im.tostring(), 
        np.uint8
        ).reshape(tuple(list(im.size)+[2]))
    if bigendian:
        return (np.array(out[:,:,0], np.uint16) << 8) + out[:,:,1]
    else:
        return (np.array(out[:,:,1], np.uint16) << 8) + out[:,:,0]
    
def nth_period(a, n=2, guess=150, width=50):
    period = guess
    for i  in range(n):
        period = (np.argmin(a[period*(i+1)-width:period*(i+1)+width])+period*(i+1)-width)/float(i+1)
    return period

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
        divpattern = join(path,'bgdiv_'+name)+'_%'+('0%dd'%len(m.group(2)))+(''.join(m.groups()[2:]))
    else:
        pattern = filename
        path, name = split(filename)
        outpattern = join(path,'treated_'+name)
        divpattern = splitext(join(path,'bgdiv_'+name))[0]+'.jpg'
        filename = pattern%0

    #read the first image
    im = imread(filename)
    if im.ndim ==0:
        im = readTIFF16(filename)
    cent = np.zeros_like(im)
    if len(sys.argv) >2:
        sigma = float(sys.argv[2])
    else:
        #estimate the best blurring radius
        finder = track.MultiscaleBlobFinder(im.shape, nbOctaves=100)
        centers = finder(im, k=1.6, Octave0=False, removeOverlap=False)
        sigma = 1.6*2**(np.argmax([
            oc.binary[1:-1].sum(-1).sum(-1) 
            for oc in finder.octaves[1:]
            ])/3.0)
    #create a monoscale Crocker & Grier finder
    finder = track.CrockerGrierFinder(im.shape)
    #background
    background = track.gaussian_filter(np.array(im, float), 10*sigma)
    background[background==0] = 1.0
    #drift or flow
    overlap = 2**int(np.log(20*sigma)/np.log(2))
    window_size = 2*overlap
    #having 3 subwindows is not enought
    use_subgrid_phasecor = 5*overlap < min(im.shape)
    if not use_subgrid_phasecor:
        window_size = im.shape[0]
    xgrid, ygrid = phaseCorrelation.get_coordinates(im.shape, window_size, overlap)
    #create drift file
    #fdrift = open(splitext(pattern%0)[0]+'.drift', 'w')
    window = np.hanning(window_size)[:,None]*np.hanning(window_size)
    #drifts = []

    for t in itertools.count():
        if not isfile(pattern%t):
            break
        im = imread(pattern%t)
        if im.ndim ==0:
            im = readTIFF16(pattern%t)
        #localisation with background removal
        centers = finder(im/background, sigma, uniform_size=0)
        #output result image
        cent.fill(0)
        for (ym,xm), (yM, xM) in zip(np.maximum(0, centers[:,:2]-sigma), np.minimum(centers[:,:2]+sigma+1, im.shape[::-1])):
            cent[xm:xM,ym:yM] += 50/256.*im.max()
        imsave(outpattern%t, np.dstack((im, cent, cent)))
        plt.imsave(divpattern%t, finder.blurred, cmap='bluetored', vmin=1, vmax=2)
        #Trajectory linking
        if t==0:
            linker = particles.Linker(len(centers))
            pos1 = centers
            if use_subgrid_phasecor:
                sp1 = np.fft.fftn(
                    phaseCorrelation.sub_windows(im/background, window_size, overlap) * window,
                    axes=[1,2]
                    )
            else:
                sp1 = np.fft.fftn(im/background*window)
        else:
            pos0 = pos1
            pos1 = centers
            #look for drift between successive pictures
            sp0 = sp1
            if use_subgrid_phasecor:
                #a phase correlation per subwindow
                sp1 = np.fft.fftn(
                    phaseCorrelation.sub_windows(im/background, window_size, overlap) * window,
                    axes=[1,2]
                    )
                R = sp0 * np.conjugate(sp1)
                u,v, confidence = np.transpose([
                    phaseCorrelation.getDisplConfidence(r) 
                    for r in np.abs(np.fft.ifftn(R/np.abs(R), axes=[1,2]))
                    ])
                #confidence-weighted average on the neighbouring windows
                weights = gaussian_filter(confidence.reshape(xgrid.shape), 1.0)-1
                u = gaussian_filter((u*(confidence-1)).reshape(xgrid.shape), 1.0)/weights
                v = gaussian_filter((v*(confidence-1)).reshape(xgrid.shape), 1.0)/weights
                #interpolate using 2D splines
                uspl = RectBivariateSpline(xgrid[0], ygrid[::-1,0], u, s=0)
                vspl = RectBivariateSpline(xgrid[0], ygrid[::-1,0], v, s=0)
                dx = vspl.ev(pos1[:,0], pos1[:,1])
                dy = uspl.ev(pos1[:,0], pos1[:,1])
            else:
                #only global phase correlation
                sp1 = np.fft.fftn(im/background*window)
                R = sp0 * np.conjugate(sp1)
                dx, dy, confidence = phaseCorrelation.getDisplConfidence(
                    np.abs(np.fft.ifftn(R/np.abs(R)))
                    )
                dx = dx*np.ones(len(pos1))
                dy = dy*np.ones(len(pos1))
            #shift = phaseCorrelation.getDispl(np.abs(np.fft.ifftn(R/np.abs(R))))
            #fdrift.write('%d %d\n'%(shift[1], shift[0]))
            #drifts.append(shift)
            pairs, distances = particles.get_links(pos0[:,:2], np.ones(len(pos0)), pos1[:,:2]+np.column_stack((dx, dy)), np.ones(len(pos1)), maxdist=6.*sigma)
            linker.addFrame(len(pos1), pairs, distances)
            if t>1:
                #try to link trajectories terminated two time steps ago with newly created trajectories (no match with one step ago)
                
                starting_pos = np.where(linker.pos2tr[-1]>linker.nbtrajs[-2])[0]
                pairs, distances = particles.get_links(
                    terminated[:,1:3], 
                    np.ones(len(terminated)), 
                    (pos1[:,:2]+np.column_stack((dx, dy)))[starting_pos], 
                    np.ones(len(starting_pos)), 
                    maxdist=6.*sigma
                    )
                if len(pairs)>0:
                    #translate pairs into indices the linker understands
                    pairs[:,0] = terminated_tr[pairs[:,0]]
                    pairs[:,1] = starting_pos[pairs[:,1]]
                    #update the linker
                    has_grown = linker.update(pairs, distances)
                    #add intermediate positions at the end of the previous time step
                    newly = pos1[has_grown[linker.pos2tr[-1]]]
                    newly[:,:2] += np.column_stack((dx, dy))[has_grown[linker.pos2tr[-1]]]
                    intermediate = 0.5 * (terminated[has_grown[terminated_tr]] + newly)
                    np.savetxt(
                        open(splitext(pattern%(t-1))[0]+'.csv', 'a'),
                        np.column_stack((np.where(has_grown)[0], intermediate)),
                        fmt='%g')
                    
                
            #remember the trajectories that end at previous time step (not continuing at present time step)
            terminated_pos_tr = [
                (tr[-1], i)
                for i, start, tr in zip(range(linker.nbtrajs[-2]), linker.trajstart, linker.tr2pos)
                if start+len(tr) == t]
            if len(terminated_pos_tr)==0:
                terminated_pos = np.zeros(0,int) 
                terminated_tr = np.zeros(0,int)
                terminated = np.zeros([0,3])
            else:
                terminated_pos, terminated_tr = np.transpose(terminated_pos_tr)
                terminated = pos0[terminated_pos]
                if use_subgrid_phasecor:
                    dx = vspl.ev(terminated[:,0], terminated[:,1])
                    dy = uspl.ev(terminated[:,0], terminated[:,1])
                    terminated[:,:2] -= np.column_stack((dx, dy))
                else:
                    terminated[:,:2] -= [dx[0], dy[0]]
        #remember inner image of the traker to measure intensities afterward
        blurred0 = np.copy(finder.blurred)
        #output
        np.savetxt(
            splitext(pattern%t)[0]+'.csv', 
            np.column_stack((linker.pos2tr[-1], centers)),
            fmt='%g'
            )
    #fdrift.close()
    #Focus on the trajectories spanning the whole time and half of it
    isspanning = np.array(map(len, linker.tr2pos)) == len(linker.pos2tr)
    halfspanning = np.intersect1d(linker.pos2tr[0], linker.pos2tr[len(linker.pos2tr)/2])
    ishalfspanning = np.zeros(len(linker.tr2pos), bool)
    ishalfspanning[halfspanning] = True
    intensities = np.zeros((len(linker.pos2tr), isspanning.sum()))
    halfinten = np.zeros((len(linker.pos2tr)/2, len(halfspanning)))
    for t, pos2tr in enumerate(linker.pos2tr):
        name = splitext(pattern%t)[0]+'.csv'
        alli = np.loadtxt(name, usecols=[3])
        intensities[t] = alli[isspanning[pos2tr]][np.argsort(pos2tr[isspanning[pos2tr]])]
        if t<len(halfinten):
            halfinten[t] = alli[ishalfspanning[pos2tr]][np.argsort(pos2tr[ishalfspanning[pos2tr]])]
            #halfinten[t] = alli[np.argsort(pos2tr[ishalfspanning[pos2tr]])]
    np.savetxt(join(path,'intensities.csv'), intensities, fmt='%g')
    np.savetxt(join(path,'half_intensities.csv'), halfinten, fmt='%g')
    #compute periods and export them with the positions in the first frame
    first = np.loadtxt(splitext(pattern%0)[0]+'.csv', usecols=[1,2])
    periods = []
    for n in itertools.count():
        try:
            periods.append([nth_period(i,1+n) for i in intensities.T])
        except ValueError:
            break
    np.savetxt(join(path,'pos_periods.csv'), np.column_stack([first[isspanning[linker.pos2tr[0]]]]+periods), fmt='%g')
    periods = []
    for n in itertools.count():
        try:
            periods.append([nth_period(i,1+n) for i in halfinten.T])
        except ValueError:
            break
    np.savetxt(join(path,'pos_halfperiods.csv'), np.column_stack([first[ishalfspanning[linker.pos2tr[0]]]]+periods), fmt='%g')
    
    
