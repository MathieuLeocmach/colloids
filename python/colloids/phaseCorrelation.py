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
from numpy.fft import rfft2, irfft2,fft2, ifft2
import numpy.lib.stride_tricks
from numpy import ma

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
    R = fft2(a)*np.conj(fft2(b))
    R /= np.absolute(R)
    return ifft2(R,a.shape)

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
    displ = np.array(np.unravel_index(l, a.shape))
    for i,v in enumerate(displ):
        if v > a.shape[i]/2:
            displ[i] -= a.shape[i]

    return displ
    
def getDisplConfidence(a):
    """Get the periodic displacement of the maximum of an array with respect to (0,0), and the confidence, expressed as the rqtio of intensity of the first by the second maximum"""
    b = np.fft.fftshift(a)
    l = b.argmax()
    first_peak = a.max()
    displ = np.array(np.unravel_index(l, a.shape))
    masked = b.view(ma.MaskedArray)
    lb = np.maximum(displ-2, 0)
    ub = np.minimum(displ+3, a.shape)
    masked[lb[0]:ub[0], lb[1]:ub[1]] = ma.masked
    second_peak = masked.max()
    displ -= np.array(a.shape)/2

    return list(displ) + [first_peak/second_peak]
    
def get_phaseshift(a,b):
    assert a.shape == b.shape
    window = np.hanning(a.shape[0])[:,None]*np.hanning(a.shape[1])
    R = np.fft.fftn(a)*np.conjugate(np.fft.fftn(b*window))
    return getDispl(np.abs(np.fft.ifftn(R/np.abs(R))))
    
def sub_windows(a, window_size, overlap ):
    """
    This is a nice numpy trick. The concept of numpy strides should be
    clear to understand this code.
    Basically, we have a 2d array and we want to perform cross-correlation
    over the interrogation windows. An approach could be to loop over the array
    but loops are expensive in python. So we create from the array a new array
    with three dimension, of size (n_windows, window_size, window_size), in which
    each slice, (along the first axis) is an interrogation window.
    """
    sz = a.itemsize
    shape = a.shape
    
    strides = (sz*shape[1]*(window_size-overlap), sz*(window_size-overlap), sz*shape[1], sz)
    shape = ( int((shape[0] - window_size)/(window_size-overlap))+1, int((shape[1] - window_size)/(window_size-overlap))+1 , window_size, window_size)
    
    return numpy.lib.stride_tricks.as_strided( a, strides=strides, shape=shape ).reshape(-1, window_size, window_size)
    
def get_field_shape ( image_size, window_size, overlap ):
    """Compute the shape of the resulting flow field.
Given the image size, the interrogation window size and
the overlap size, it is possible to calculate the number
of rows and columns of the resulting flow field.
Parameters
----------
image_size: two elements tuple
a two dimensional tuple for the pixel size of the image
first element is number of rows, second element is
the number of columns.
window_size: int
the size of the interrogation window.
overlap: int
the number of pixel by which two adjacent interrogation
windows overlap.
Returns
-------
field_shape : two elements tuple
the shape of the resulting flow field
"""
    
    return ( (image_size[0] - window_size)//(window_size-overlap)+1,
             (image_size[1] - window_size)//(window_size-overlap)+1 )
             
def get_coordinates( image_size, window_size, overlap ):
    """Compute the x, y coordinates of the centers of the interrogation windows.
Parameters
----------
image_size: two elements tuple
a two dimensional tuple for the pixel size of the image
first element is number of rows, second element is
the number of columns.
window_size: int
the size of the interrogation windows.
overlap: int
the number of pixel by which two adjacent interrogation
windows overlap.
Returns
-------
x : 2d np.ndarray
a two dimensional array containing the x coordinates of the
interrogation window centers, in pixels.
y : 2d np.ndarray
a two dimensional array containing the y coordinates of the
interrogation window centers, in pixels.
"""

    # get shape of the resulting flow field
    field_shape = get_field_shape( image_size, window_size, overlap )

    # compute grid coordinates of the interrogation window centers
    x = np.arange( field_shape[1] )*(window_size-overlap) + (window_size-1)/2.0
    y = np.arange( field_shape[0] )*(window_size-overlap) + (window_size-1)/2.0
    
    return np.meshgrid(x,y[::-1])
    
if __name__ == "__main__":
    #aquire initial data
    reader = lif.Reader("D:\Users\ishizaka\Documents\dataMathieu\Tsuru11dm_phi=52.53_J36.lif")
    serie = reader.getSeries()[3]
    size2D = tuple(serie.getNumberOfElements()[:2])
    a0 = np.fromstring(serie.get2DString(T=0),np.ubyte).reshape(size2D[1],size2D[0])
    a1 = np.fromstring(serie.get2DString(T=29),np.ubyte).reshape(size2D[1],size2D[0])

    #in fact, windowing the inputs produces strong noise
    #r = phaseCorrel(windowFunction(a0),windowFunction(a1))
    r = phaseCorrel(a0,a1)

    print getDispl(r)

    showNormalized(r)
