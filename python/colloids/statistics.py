import numpy as np
import Gnuplot

def plothist(data, title=None, bins=50, normed=False):
    h, b = np.histogram(data, bins=bins, normed=normed)
    return Gnuplot.Data(
        b[:-1], h,
        title=title,
        with_='steps'
        )

def plot2dhist(datax, datay, title=None, bins=50, normed=False):
    h, bx, by = np.histogram2d(datax, datay, bins=bins, normed=normed)
    return Gnuplot.GridData(
        h, bx[:-1], by[:-1],
        with_='lines',
        title=title)

def meanMap(x, y, values, bins=50):
    number, bx, by = np.histogram2d(x, y, bins)
    total, bx, by = np.histogram2d(x, y, bins=[bx, by], weights=values)
    return np.where(number==0, -1, total/number), bx[:-1], by[:-1]

def plotMeanMap(x, y, values, bins=50, cbmax=None):
    """plot the mean map in order to have real step functions"""
    h, bx, by = meanMap(x, y, values, bins)
    h2 = np.repeat(np.repeat(h, 2, axis=0), 2, axis=1)
    bx2 = np.repeat(bx, 2)
    bx2[1::2] = bx+(bx[1]-bx[0])
    by2 = np.repeat(by, 2)
    by2[1::2] = by+(by[1]-by[0])
    if not cbmax:
            cbmax = h.max()
    low = np.extract(h>-1, h).min();
    g('set palette defined (%g "white", %g "black", %g "purple", %g "red", %g "yellow")' % (low-(cbmax-low)/100, low, (cbmax+2*low)/3, (2*cbmax+low)/3, cbmax))
    g('set cbrange [%g:%g]' % (low-(cbmax-low)/100, cbmax))
    g.splot(Gnuplot.GridData(h2, bx2, by2, binary=0))

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
