import numpy as np
from colloids import track
from scipy.misc import imsave, imshow, imread
import sys, itertools, re
from os.path import splitext,isfile

if len(sys.argv) >1:
    filename = sys.argv[1]
else:
    while(True):
        filename = raw_input("filename --> ")
        if isfile(filename):
            break
        else:
            print '"%s" is not an existing file' % filename

m = re.match('(.*)_([0-9]*)(.*)', filename)
assert m is not None, '"%s" do not match pattern' % filename
pattern = m.group(1)+'_%'+('0%dd'%len(m.group(2)))+(''.join(m.groups()[2:]))

im = imread(filename)
cent = np.zeros_like(im)
finder = track.MultiscaleBlobFinder(im.shape, nbOctaves=100)
f = open(m.group(1)+'.nb', 'w')
for t in itertools.count():
    if not isfile(pattern%t):
        break
    im = imread(pattern%t)
    centers = finder(im, Octave0=False, removeOverlap=False)
    centers = centers[centers[:,-1]<-1]
    f.write('%d\n'%len(centers))
    cent.fill(0)
    for x, y, r, i in centers:
        cent[y-r:y+r+1,x-r:x+r+1] += 50
    imsave('treated_'+pattern%t, np.dstack((im, cent, cent)))
f.close()
    