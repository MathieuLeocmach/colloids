#
#    Copyright 2013 Mathieu Leocmach
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

import numpy as np
import scipy.constants as const
from scipy.special import gamma
import warnings

def logderive(x,f,width):
    """A special purpose routine for finding the first and second logarithmic derivatives of slowly varying, but noisy data.  It returns a smoother version of 'f' and its first and second log derivative."""
    lx = np.log(x)
    ly = np.log(f)
    f2 = np.zeros_like(f)
    df = np.zeros_like(f)
    ddf = np.zeros_like(f)
    
    for i, u in enumerate(lx):
        #Data are unequally spaced, so we cannot use a simple gaussian_filter
        w = np.exp(-(lx-u)**2 / (2*width**2))
        #truncate the Gaussian, run faster
        ww = w > 0.03
        #polynomial fit weighted by the Gaussian.
        p = np.poly1d(np.polyfit(lx[ww], ly[ww], 2, w=w[ww]))
        f2[i] = np.exp(p(u))
        df[i] = p.deriv(1)(u)
        ddf[i] = p.deriv(2)(u)
    return f2, df, ddf



def msd2G(dt, msd, a, T, dim=3, clip=0.03, width=0.7):
    """One particle microrheology. Returns omega (1/s), G (Pa)
    
    Keyword arguments:
    dt -- in seconds
    msd -- in microns^2
    a -- radius in microns
    T -- in Kelvin
    dim -- the dimensionality of the data
    clip -- a fraction of G(s) below which G'(w) and G"(w) are meaningless (see NOTES)
    width -- the width of the gaussian that is used for the polyfit (see Notes below)
    
    NOTES:
        It needs more than a 7-8 points per decade of time/frequency
        and really hates noise and long-wavelength ripples in the data.

        G'(w) and G"(w) are clipped at 0.03x G(s) as they are almost
        certainly meaningless below that point unless the data is
        *extremely* clean.  Set 'clip' to less than 0.03 to see more.
        See Tom Mason's paper: PRL *79*, 3284, (1997) for details.

        set the width to something bigger for noisy data, aok if the data
        is very weakly curved. recommended starting value: 0.7
        
        This function is the translation of the IDL implementation by John C. Crocker in 1999
    """
    #convert the inputs in meters - Kilograms - Seconds
    am = a * 1e-6
    omega = 1./ dt
    msdm = msd * 1e-12
    C = dim * constants.k * T / (3 * np.pi * am)
    
    #use 2nd order local formula for G(s)-- good to 1% of G(s)
    m, d, dd = logderive(dt, msdm, width)
    Gs = C / (m * gamma(1+d) * (1 + dd/2))
    
    #use 2nd order local formula for G'(w),G"(w)-- good to 2% of G(w)
    g, da, dda = logderive(omega, Gs, width)
    G  = g / (1+dda) * (np.exp(0.5j * np.pi * da) - (np.pi/2 - 1) * (1-da) * dda)
    
    if (np.abs(dd).max() > 0.15) or (np.abs(dda).max() > 0.15):   #original value: 0.15
        warnings.warn('High curvature in data, moduli may be unreliable!')
    
    #clip off the suspicious (i.e. unreliable) data
    w = G.real < Gs*clip
    if n.sum() > 0: 
        G.real[w]=0
    w = G.imag < Gs*clip
    if w.sum() > 0: 
        G.imag[w]=0
        
    return omega, G
