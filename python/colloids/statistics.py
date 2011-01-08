import numpy as np
import Gnuplot

g=Gnuplot.Gnuplot()
g('set pm3d map')

def plothist(data, title=None, bins=50, normed=False):
    h, b = np.histogram(data, bins=bins, normed=normed)
    return Gnuplot.Data(
        b[:-1], h,
        title=title,
        with_='steps'
        )

def plot2dhist(datax, datay, title=None, bins=50, normed=False, cbmax=None, cbmin=None, logscale=False):
    h, bx, by = np.histogram2d(datax, datay, bins=bins, normed=normed)
    h[np.where(h==0)] = -1
    if logscale:
        h[np.where(h>0)] = np.log10(h[np.where(h>0)])
    h2 = np.repeat(np.repeat(h, 2, axis=0), 2, axis=1)
    bx2 = np.repeat(bx[:-1], 2)
    bx2[1::2] = bx[:-1]+(bx[1]-bx[0])
    by2 = np.repeat(by[:-1], 2)
    by2[1::2] = by[:-1]+(by[1]-by[0])
    if not cbmax:
            cbmax = h.max()
    if cbmin:
        low = cbmin
    else:
        low = np.extract(h>-1, h).min();
    g('set palette defined (%g "white", %g "black", %g "purple", %g "red", %g "yellow")' % (low-(cbmax-low)/100, low, (cbmax+2*low)/3, (2*cbmax+low)/3, cbmax))
    g('set cbrange [%g:%g]' % (low-(cbmax-low)/100, cbmax))
    g.splot(Gnuplot.GridData(h2, bx2, by2, binary=0))

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
    k = data
    if k.dtype.kind != 'c':
        k -= data.mean()
    return (k*np.conj(k[0])).mean(axis=1)/(np.abs(k)**2).mean()

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
