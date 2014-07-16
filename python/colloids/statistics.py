import numpy as np
import Gnuplot
import numexpr
from scipy import weave

g=Gnuplot.Gnuplot()
g('set pm3d map')

def plothist(data, title=None, bins=50, normed=False):
    h, b = np.histogram(data, bins=bins, normed=normed)
    return Gnuplot.Data(
        b[:-1], h,
        title=title,
        with_='steps'
        )

def columnData(data, xvals=None, yvals=None, **keywords):
    """Prepare the arguments needed for an usual Gnuplot grid data in order to have a column data (no interpolation)"""
    if xvals is None:
        xvals = np.arange(data.shape[0])
    if yvals is None:
        yvals = np.arange(data.shape[1])
    data = np.repeat(np.repeat(data, 2, axis=0), 2, axis=1)
    xvals2 = np.repeat(xvals, 2)
    xvals2[1:-1:2] = xvals[1:]
    xvals2[-1] = xvals[-1]+2*(xvals[-1]-xvals[-2])
    yvals2 = np.repeat(yvals, 2)
    yvals2[1:-1:2] = yvals[1:]
    yvals2[-1] = yvals[-1]+2*(yvals[-1]-yvals[-2])
    keywords['binary'] = False
    return Gnuplot.GridData(data, xvals2, yvals2, **keywords)

def plot2dhist(datax, datay, title=None, bins=50, normed=False, cbmax=None, cbmin=None, logscale=False):
    h, bx, by = np.histogram2d(datax, datay, bins=bins, normed=normed)
    if logscale:
        h[np.where(h>0)] = np.log10(h[np.where(h>0)])
    h2 = np.repeat(np.repeat(h, 2, axis=0), 2, axis=1)
    bx2 = np.repeat(bx[:-1], 2)
    bx2[1::2] = bx[:-1]+(bx[1]-bx[0])
    by2 = np.repeat(by[:-1], 2)
    by2[1::2] = by[:-1]+(by[1]-by[0])
    g('set palette gray negative')
    g.splot(Gnuplot.GridData(h2, bx2, by2, binary=0))

def digitized2dhist(datax, datay, bins=50, zbins=5, logscale=True):
    h, bx, by = np.histogram2d(datax, datay, bins=bins)
    #test iterability
    try:
        [e for e in zbins]
        return -np.rot90(np.digitize(h.ravel(), zbins).reshape(h.shape))
    except TypeError:
        if logscale:
            return -np.rot90(np.digitize(
                h.ravel(),
                np.logspace(0, np.log10(h.max()), zbins)
                ).reshape(h.shape))
        else:
            return np.digitize(
                h.ravel(),
                np.linspace(1, h.max(), zbins)
                ).reshape(h.shape)

def meanMap1D(x, values, bins=50):
    number, b = np.histogram(x, bins)
    total, b = np.histogram(x, bins=b, weights=values)
    selection = np.where(number>0)
    return b[selection], total[selection]/number[selection], number[selection]

def meanMap(x, y, values, bins=50):
    number, bx, by = np.histogram2d(x, y, bins)
    total, bx, by = np.histogram2d(x, y, bins=[bx, by], weights=values)
    return np.where(number==0, -1, total/number), bx[:-1], by[:-1]

def plotMeanMap(x, y, values, bins=50, cbmax=None, cbmin=None):
    """plot the mean map in order to have real step functions"""
    h, bx, by = meanMap(x, y, values, bins)
    h2 = np.repeat(np.repeat(h, 2, axis=0), 2, axis=1)
    bx2 = np.repeat(bx, 2)
    bx2[1::2] = bx+(bx[1]-bx[0])
    by2 = np.repeat(by, 2)
    by2[1::2] = by+(by[1]-by[0])
    if not cbmax:
            cbmax = h.max()
    if cbmin:
        low = cbmin
    else:
        low = np.extract(h>-1, h).min();
    g(
        'set palette defined ('
        +'%g "white", %g "black", %g "purple", %g "red", %g "yellow")'
        % (
            low-(cbmax-low)/100,
            low,
            (2*cbmax+3*low)/5,
            (3*cbmax+2*low)/5,
            cbmax
            )
        )
    g('set cbrange [%g:%g]' % (low-(cbmax-low)/100, cbmax))
    g.splot(Gnuplot.GridData(h2, bx2, by2, binary=0))

def timecorrel(data):
    """time correlation of an array of size [time, particles"""
    k = np.copy(data)
    if k.dtype.kind != 'c':
        k -= data.mean()
    p = numexpr.evaluate(
        'k * complex(k0.real, -k0.imag)',
        {'k':k, 'k0':k[0]}
        ).mean(axis=1)
    n = numexpr.evaluate('abs(k).real**2').mean()
    return numexpr.evaluate('p/n')

def time_correlation(data, av=10):
    """read the particle-wise scalar from a time serie of files and compute the time correlation"""
    if av==0:
        c = np.zeros_like(data.sum(axis=1))
        for t0 in range(len(data)-1):
            c[:len(data)-t0] += timecorrel(data[t0:])
        for dt in range(len(c)):
            c[dt] /= len(c)-dt
        return c
    else:
        return np.mean([timecorrel(data[t0:len(data)+t0-av]) for t0 in range(av)], axis=0)

def correlate_qlm(qlm, average=True):
    qlm_cor = np.real(qlm*np.conj(qlm[0]))[:,:,[0]+2*range(1,qlm.shape[-1])].sum(axis=-1)
    if average:
            return (qlm_cor/qlm_cor[0]).mean(axis=-1)
    else:
            return qlm_cor/qlm_cor[0]

def time_correlation_qlm(qlm, av=10):
    if av==0:
        c = np.zeros(len(qlm), float)
        for t0 in range(len(qlm)-1):
            c[:len(qlm)-t0] += correlate_qlm(qlm[t0:])
        for dt in range(len(c)):
            c[dt] /= len(c)-dt
        return c
    else:
        return np.mean([correlate_qlm(qlm[t0:len(qlm)+t0-av]) for t0 in range(av)], axis=0)


def space_correlation(positions, values, bins):
    """spatial correlation of a scalar observable"""
    v = values - values.mean()
    total = np.zeros(len(bins)-1, values.dtype)
    rdf = np.zeros_like(total)
    for p, u in zip(positions, np.conjugate(v)):
	distsq = np.sum((positions-p)**2, axis=-1)
	h = np.histogram(distsq, bins=bins, weights=v*u)[0]
	n = np.histogram(distsq, bins=bins)[0]
	rdf += n
	nonzero = n>0
	total[nonzero] += h[nonzero]/n[nonzero]
    return total/(v*np.conj(v)).sum(), rdf/bins[1:]**2


def pca(data):
    cov = np.cov(data.T)
    evalues, evect = np.linalg.eig(cov)
    s = (-np.absolute(evalues)).argsort()
    evect = evect[s]
    evalues = evalues[s]
    nored = np.dot(evect.T, data.T)
    g = Gnuplot.Gnuplot()
    d = Gnuplot.Data(nored[0], nored[1])
    g.plot(d)
    red = [np.dot(evect[:,:i].T, data.T) for i in range(1,len(evect)+1)]
    redback = np.asarray(
            [np.dot(evect[:,:i+1], r).T for i,r in enumerate(red)]
            )
    raw_input('ok ?')
    return redback

def dca(data, groups):
    Ws = np.asarray([np.cov(data[g].T) for g in groups])
    W = np.average(Ws, axis=0, weights=map(len, groups))
    mus = np.asarray([data[g].mean(axis=0) for g in groups])
    Bs = np.asarray([np.outer(m-mu, m-mu) for m in mus])
    B = np.average(Bs, axis=0, weights=map(len, groups))
    V = B+W
    evalues, evect = np.linalg.eig(np.dot(np.linalg.inv(V), B))
    evalues=np.real_if_close(evalues)
    s = (-np.absolute(evalues)).argsort()
    evalues = evalues[s]
    evect = evect[s]
    nored = np.dot(evect.T, cloudSC.T)
    d = [Gnuplot.Data(nored[0][g], nored[1][g]) for g in groups]
    g = Gnuplot.Gnuplot()
    eval('g.plot('+', '.join(['d[%d]'%i for i in range(len(d))])+')')
    raw_input('ok?')
    return evalues, evect
    
class StructureFactor2D:
    def __init__(self, shape):
        self.im = np.zeros(shape)
        #mask of wavenumbers
        self.wavenumbers = map(np.fft.fftfreq, shape)
        self.dists = numexpr.evaluate('sqrt(qx**2+qy**2)', {
            'qx':self.wavenumbers[0][:,None],
            'qy':self.wavenumbers[1][None,:]
            })
        #bin the wavenumbers
        self.nbq, self.qs = np.histogram(self.dists.ravel(), self.wavenumbers[0][:len(self.wavenumbers[0])/2])
        self.S = np.zeros(self.nbq.shape)
        self.ws = map(np.hamming, shape)
        self.Ns = []
        
    def __call__(self, pos, accumulate=True):
        if len(pos)==0:
            self.Ns.append(0)
            return
        #draw a white pixel at the position of each particle
        self.im.fill(0)
        for x, y, z in pos:
            self.im[x,y,z] = 1
        #remove offset
        self.im -= self.im.mean()
        #windowing
        for d, w in enumerate(self.ws):
            self.im *= w[tuple([None]*d + [slice(None)] + [None]*(2-d))]
        #do the (half)Fourier transform
        spectrum = np.abs(np.fft.rfftn(self.im)[:,:,0])**2
        #radial average (sum)
        S = np.histogram(self.dists.ravel(), self.qs, weights=spectrum.ravel())[0]/spectrum.mean()
        if accumulate:
            self.S += S*len(pos)
            self.Ns.append(len(pos))
        else:
            return S
        
    def get_S(self):
        #radial average (division)
        S = np.copy(self.S)
        S[self.nbq>0] /= self.nbq[self.nbq>0]
        S[0] = np.var(self.Ns) * len(self.Ns)
        return S
        
class StructureFactor3D:
    def __init__(self, shape):
        self.im = np.zeros(shape)
        #mask of wavenumbers
        self.wavenumbers = map(np.fft.fftfreq, shape)
        self.dists = numexpr.evaluate('sqrt(qx**2+qy**2+qz**2)', {
            'qx':self.wavenumbers[0][:,None,None],
            'qy':self.wavenumbers[1][None,:,None],
            'qz':self.wavenumbers[1][None,None,:shape[-1]/2+1]
            })
        self.spectrum = np.zeros_like(self.dists)
        #bin the wavenumbers
        self.nbq, self.qs = np.histogram(self.dists.ravel(), self.wavenumbers[0][:len(self.wavenumbers[0])/2])
        self.S = np.zeros(self.nbq.shape)
        self.ws = map(np.hamming, shape)
        self.Ns = []
        
    def __call__(self, pos, periodic=False, accumulate=True):
        if len(pos)==0:
            self.Ns.append(0)
            return
        #draw a white pixel at the position of each particle
        self.im.fill(0)
        for x, y, z in pos:
            self.im[x,y,z] = 1
        #remove offset
        self.im -= self.im.mean()
        if not periodic:
            #windowing
            for d, w in enumerate(self.ws):
                self.im *= w[tuple([None]*d + [slice(None)] + [None]*(2-d))]
        #do the (half)Fourier transform
        self.spectrum = numexpr.evaluate(
            'abs(f).real**2', 
            {'f': np.fft.rfftn(self.im)}
            )

        #radial average (sum)
        S = np.histogram(self.dists.ravel(), self.qs, weights=self.spectrum.ravel())[0]/self.spectrum.mean()
        if accumulate:
            self.S += S*len(pos)
            self.Ns.append(len(pos))
        else:
            return S
        
    def get_S(self):
        #radial average (division)
        S = np.copy(self.S)
        S[self.nbq>0] /= self.nbq[self.nbq>0]
        S[0] = np.var(self.Ns) * len(self.Ns)
        return S
        
        
class ImageStructureFactor:
    """A class to compute rapially averaged structure factor of a 2D image"""
    def __init__(self, shape):
        assert len(shape) == 2, "only 2D images implemented"
        L = max(shape)
        self.qs = np.fft.fftfreq(L)[:L/2]
        self.dists = np.sqrt(np.fft.fftfreq(shape[0])[:,None]**2 + np.fft.fftfreq(shape[1])**2)
        self.dcount = np.histogram(self.dists.ravel(), bins=self.qs)[0]
        self.has_window = False
        
    def set_window(self,w='hanning'):
        if hasattr(np, w):
            self.window = getattr(np,w)(self.dists.shape[0])[:,None] * getattr(np,w)(self.dists.shape[1])[:,None]
            self.has_window = True
        elif w == False:
            self.has_window = False
        elif isinstance(w, np.ndarray):
            assert w.shape == self.dists.shape
            self.window = w
            self.has_window = True
        
    def __call__(self, im):
        if self.has_window:
            spectrum = np.abs(np.fft.fft2(im*self.window))**2
        else:
            spectrum = np.abs(np.fft.fft2(im))**2
        return np.histogram(self.dists.ravel(), bins=self.qs, weights=spectrum.ravel())[0] / self.dcount
        
        
