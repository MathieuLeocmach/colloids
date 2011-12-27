#
#    Copyright 2011 Mathieu Leocmach
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
import rtree.index
from scipy.special import sph_harm
from scipy import weave
from scipy.weave import converters
        
        
def weave_qlm(pos, ngbs, inside, l=6):
    qlm = np.zeros([len(pos), l+1], np.complex128)
    support = """
    #include <boost/math/special_functions/spherical_harmonic.hpp>
    """
    code = """
    #pragma omp parallel for
    for(int i=0; i<Nngbs[0]; ++i)
    {
        if(!inside(i))
            continue;
        double cart[3] = {0, 0, 0};
        double sph[3] = {0, 0, 0};
        int nb = 0;
        for(int j=0; j<Nngbs[1]; ++j)
        {
            int q = ngbs(i,j);
            if(q<0 || q>=Npos[1])
                continue;
            nb++;
            for(int d=0; d<3; ++d)
                cart[d] = pos(i,d) - pos(q, d);
            sph[0] = sqrt(cart[0]*cart[0] + cart[1]*cart[1] + cart[2]*cart[2]);
            if(abs(cart[2])==sph[0] || sph[0]*sph[0]+1.0 == 1.0)
            {
                sph[1] = 0;
                sph[2] = 0;
            }
            else
            {
                sph[1] = acos(cart[2]/sph[0]);
                sph[2] = atan2(cart[1], cart[0]);
                if(sph[2]<0)
		            sph[2] += 2.0*M_PI;
            }
		    for(int m=0; m<Nqlm[1]; ++m)
		        qlm(i,m) += boost::math::spherical_harmonic(Nqlm[1]-1, m, sph[1], sph[2]);
        }
        for(int m=0; m<Nqlm[1]; ++m)
            qlm(i,m) /= nb;
    }
    """
    weave.inline(
        code,['pos', 'ngbs', 'inside', 'qlm'],
        type_converters =converters.blitz,
        support_code = support,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return qlm
    


def cart2sph(cartesian):
    """Convert Cartesian coordinates [[x,y,z],] to spherical coordinates [[r,phi,theta],]
phi is cologitudinal and theta azimutal"""
    spherical = np.zeros_like(cartesian)
    #distances
    spherical[:,0] = np.sum(cartesian**2, axis=-1)
    sel = np.asarray([
        r!=z or r+1.0>1.0
        for z, r in zip(cartesian[:,-1]**2, spherical[:,0])
        ])
    spherical[:,0] = np.sqrt(spherical[:,0])
    #colatitudinal phi [0, pi[
    spherical[sel,1] = np.arccos(cartesian[sel,-1]/spherical[sel,0])
    #azimutal (longitudinal) theta [0, 2pi[
    spherical[sel,2] = np.arctan2(cartesian[sel,1], cartesian[sel,0])
    spherical[spherical[:,2]<0, 2] += 2*np.pi
    return spherical
    

def ql(qlm):
    l = qlm.shape[1]-1
    s = 2*np.sum(np.abs(qlm[:,1:])**2, -1) + np.abs(qlm[:,0])**2
    return np.sqrt(4*np.pi/(2*l+1)*s)
def wl(qlm):
    l = qlm.shape[1]-1
    w = np.zeros(qlm.shape[0], qlm.dtype)
    for m1 in range(-l, l+1):
        for m2 in range(-l, l+1):
            m3 = -m1-m2
            if -l<=m3 and m3<=l:
                w+= get_w3j(l, [m1, m2, m3]) * get_qlm(qlm, m1) * get_qlm(qlm, m2) * get_qlm(qlm, m3)
    return w
    
def get_qlm(qlms, m):
    if m>=0:
        return qlms[:,m]
    if (-m)%2 == 0:
        return np.conj(qlms[:,-m])
    return -np.conj(qlms[:,-m])
		
_w3j = [
    [1],
    np.sqrt([2/35., 1/70., 2/35., 3/35.])*[-1,1,1,-1],
    np.sqrt([
        2/1001., 1/2002., 11/182., 5/1001.,
        7/286., 5/143., 14/143., 35/143., 5/143.,
        ])*[3, -3, -1/3.0, 2, 1, -1/3.0, 1/3.0, -1/3.0, 1],
    np.sqrt([
        1/46189., 1/46189.,
        11/4199., 105/46189.,
        1/46189., 21/92378.,
        1/46189., 35/46189., 14/46189.,
        11/4199., 21/4199., 7/4199.,
        11/4199., 77/8398., 70/4199., 21/4199.
        ])*[-20, 10, 1, -2, -43/2.0, 3, 4, 2.5, -6, 2.5, -1.5, 1, 1, -1, 1, -2],
    np.sqrt([
        10/96577., 5/193154.,
        1/965770., 14/96577.,
        1/965770., 66/482885.,
        5/193154., 3/96577., 77/482885.,
        65/14858., 5/7429., 42/37145.,
        65/14858., 0.0, 3/7429., 66/37145.,
        13/74290., 78/7429., 26/37145., 33/37145.,
        26/37145., 13/37145., 273/37145., 429/37145., 11/7429.,
        ])*[
            7, -7, -37, 6, 73, -3,
            -5, -8, 6, -1, 3, -1,
            1, 0, -3, 2, 7, -1, 3, -1,
            1, -3, 1, -1, 3],
    np.sqrt([
        7/33393355., 7/33393355.,
        7/33393355., 462/6678671.,
        7/33393355., 1001/6678671.,
        1/233753485., 77/6678671., 6006/6678671.,
        1/233753485., 55/46750697., 1155/13357342.,
        1/233753485., 2926/1757545., 33/46750697., 3003/6678671.,
        119/1964315., 22/2750041., 1914/474145., 429/5500082.,
        17/13750205., 561/2750041., 77/392863., 143/27500410., 2002/392863.,
        323/723695., 1309/20677., 374/144739., 143/144739., 1001/206770.,
        323/723695., 7106/723695., 561/723695., 2431/723695., 2002/103385., 1001/103385.
        ])*[
            -126, 63, 196/3.0, -7, -259/2.0, 7/3.0,
            1097/3.0, 59/6.0, -2,
            4021/6.0, -113/2.0, 3,
            -914, 1/3.0, 48, -3,
            -7/3.0, 65/3.0, -1, 3,
            214/3.0, -3, -2/3.0, 71/3.0, -1,
            3, -1/3.0, 5/3.0, -2, 1/3.0,
            2/3.0, -1/3.0, 2, -4/3.0, 2/3.0, -1]
    ]
_w3j_m1_offset = np.array([0,1,2,4,6,9,12,16,20,25,30], int)

def get_w3j(l, ms):
    sm = np.sort(np.abs(ms))
    return _w3j[l/2][_w3j_m1_offset[sm[-1]]+sm[0]]
