import numpy as np
from scipy.misc import imshow, imread
import sys, itertools, re
from os.path import split,join,isfile
import numexpr

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
shape = im.shape
wavenumbers = map(np.fft.fftfreq, shape)
dists = numexpr.evaluate('sqrt(qx**2+qy**2)', {
    'qx':wavenumbers[0][:,None],
    'qy':wavenumbers[1][None,:shape[-1]/2+1]
    })
nbq, qs = np.histogram(dists.ravel(), wavenumbers[0][:len(wavenumbers[0])/2])
ws = map(np.hamming, shape)

for t in itertools.count():
    if not isfile(pattern%t):
        break
    im = imread(pattern%t)
    im = im - im.mean()
    for d, w in enumerate(ws):
        im *= w[tuple([None]*d + [slice(None)] + [None]*(im.ndim-d-1))]
    spectrum = numexpr.evaluate(
        'abs(f).real**2', 
        {'f': np.fft.rfftn(im)}
        )
    Sq = np.histogram(dists.ravel(), qs, weights=spectrum.ravel())[0]
    Sq[nbq>0] /= nbq[nbq>0]
    np.savetxt(outpattern%t, Sq)

