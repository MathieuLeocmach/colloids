import numpy as np
from scipy.misc import imshow, imread
import sys, itertools, re
from os.path import split,join,isfile

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
path, name = split(m.group(1))
outpattern = join(path,'spectrum_'+name)+'_%'+('0%dd'%len(m.group(2)))+'.dat'

im = imread(filename)
dists = np.zeros(im.shape)
for x, dx in enumerate((np.arange(im.shape[0])-im.shape[0]/2)**2):
	for y, dy in enumerate((np.arange(im.shape[1])-im.shape[1]/2)**2):
		dists[x,y] = dx+dy
dists = np.sqrt(dists)

for t in itertools.count():
    if not isfile(pattern%t):
        break
    im = imread(pattern%t)
    sp = np.fft.fftshift(np.fft.fft2(im))
    Sq = np.histogram(dists.ravel(), bins=range(sp.shape[0]/2), weights=np.abs(sp).ravel())[0]
    nb = np.histogram(dists.ravel(), bins=range(sp.shape[0]/2))[0]
    Sq[nb>0] /= nb[nb>0]
    np.savetxt(outpattern%t, Sq)

