#
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

# Implementing
# http://en.wikipedia.org/wiki/Phase_correlation

import lif
import Image
import numpy as np
from numpy.fft import rfft2, irfft2

def inplaceWindowFunction(a):
    """apply a Hamming window function to each dimension of the input array"""
    for i,v in enumerate(a.shape):
        u = np.ones((len(a.shape)))
        u[i] = v
        a *= np.hamming(v).reshape(u)
    return a

def windowFunction(a):
    """return a copy of the input with a Hamming window function applied to each of it's dimension"""
    return inplaceWindowFunction(np.copy(a))

def phaseCorrel(a,b):
    """phase correlation calculation"""
    R = rfft2(a)*np.conj(rfft2(b))
    R /= np.absolute(R)
    return irfft2(R,a.shape)

def showNormalized(a):
    """normalizing array to display as image"""
    m = np.amin(a)
    M = np.amax(a)
    n = 255*(a-m)/(M-m)
    im = Image.fromarray(np.uint8(n))
    im.show()

def getDispl(a):
    """Get the periodic displacement of the maximum of an array with respect to (0,0)"""
    l = a.argmax()
    displ = [l / a.shape[0], l % a.shape[0]]
    for i,v in enumerate(displ):
        if v > a.shape[i]/2:
            displ[i] -= a.shape[i]

    return displ

#aquire initial data
reader = lif.LifReader("D:\Users\ishizaka\Documents\dataMathieu\Tsuru11dm_phi=52.53_J36.lif")
serie = reader.getSeries()[3]
size2D = tuple(serie.getNumberOfElements()[:2])
a0 = np.fromstring(serie.get2DString(T=0),np.ubyte).reshape(size2D[1],size2D[0])
a1 = np.fromstring(serie.get2DString(T=29),np.ubyte).reshape(size2D[1],size2D[0])

#in fact, windowing the inputs produces strong noise
#r = phaseCorrel(windowFunction(a0),windowFunction(a1))
r = phaseCorrel(a0,a1)

print getDispl(r)

showNormalized(r)
