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
import numexpr
        
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
            if(q<0 || q>=Npos[0])
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
        if(nb>0)
            for(int m=0; m<Nqlm[1]; ++m)
                qlm(i,m) /= (double)nb;
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
    
def boo_product(qlm1, qlm2):
    if qlm1.ndim==2 and qlm2.ndim==2:
        prod = np.empty([len(qlm1), len(qlm2)])
        code ="""
        #pragma omp parallel for
        for(int p=0; p<Nqlm1[0]; ++p)
            for(int q=0; q<Nqlm2[0]; ++q)
            {
                prod(p,q) = real(qlm1(p,0)*conj(qlm2(q,0)));
                for(int m=1; m<Nqlm1[1]; ++m)
                    prod(p,q) += 2.0*real(qlm1(p,m)*conj(qlm2(q,m)));
                prod(p,q) *= 4.0*M_PI/(2.0*(Nqlm1[1]-1)+1);
            }
        """
        weave.inline(
            code,['qlm1', 'qlm2', 'prod'],
            type_converters =converters.blitz,
            extra_compile_args =['-O3 -fopenmp'],
            extra_link_args=['-lgomp'],
            verbose=2, compiler='gcc')
        return prod
    else:
        n = np.atleast_2d(numexpr.evaluate(
            """real(complex(real(a), -imag(a)) * b)""",
            {'a':qlm1, 'b':qlm2}
            ))
        p = numexpr.evaluate(
            """4*pi/(2*l+1)*(2*na + nb)""",
            {
                'na': n[:,1:].sum(-1),
                'nb': n[:,0],
                'l': n.shape[1]-1,
                'pi': np.pi
                })
        return p

def ql(qlm):
    q = np.empty(len(qlm))
    code = """
    #pragma omp parallel for
    for(int i=0; i<Nqlm[0]; ++i)
    {
        q(i) = std::norm(qlm(i,0));
        for(int m=1; m<Nqlm[1]; ++m)
            q(i) += 2.0 * std::norm(qlm(i,m));
        q(i) = sqrt(4*M_PI / (2*(Nqlm[1]-1)+1) * q(i));
    }
    """
    weave.inline(
        code,['qlm', 'q'],
        type_converters =converters.blitz,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return q
    
def wl(qlm):
    support = """
    inline std::complex<double> get_qlm(const blitz::Array<std::complex<double>,2> &qlms, const int &i, const int &m)
    {
        if(m>=0)
            return qlms(i,m);
        if((-m)%2 == 0)
            return std::conj(qlms(i,-m));
        return -std::conj(qlms(i,-m));
    }
    inline const double& getW3j(const blitz::Array<double, 1> &w3j, const int &m1, const int &m2)
    {
        const int w3j_m1_offset[11] = {0,1,2,4,6,9,12,16,20,25,30};
        std::list<int> ms;
        ms.push_back(abs(m1));
        ms.push_back(abs(m2));
        ms.push_back(abs(m1+m2));
        ms.sort();
        const int off = w3j_m1_offset[ms.back()]+ms.front();
        return w3j(off);
    }
    """
    w = np.zeros(qlm.shape[0])
    code = """
    const int l = Nqlm[1]-1;
    #pragma omp parallel for
    for(int i=0; i<Nqlm[0]; ++i)
        for(int m1=-l; m1<=l; ++m1)
        {
            std::complex<double> qlm1 = get_qlm(qlm, i, m1);
            for(int m2=-l; m2<=l; ++m2)
            {
                int m3 = -m1-m2;
                if( m3<-l || m3>l)
                    continue;
                std::complex<double> qlm2 = get_qlm(qlm, i, m2), 
                    qlm3 = get_qlm(qlm, i, m3);
                w(i) += std::real(getW3j(w3j, m1, m2) * qlm1 * qlm2 * qlm3);
            }
        }
    """
    w3j = _w3j[(qlm.shape[1]-1)/2]
    weave.inline(
        code,['qlm', 'w', 'w3j', '_w3j_m1_offset'],
        type_converters =converters.blitz,
        support_code = support,
        headers = ['<list>'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return w
    
def get_qlm(qlms, m):
    if m>=0:
        return qlms[:,m]
    if (-m)%2 == 0:
        return np.conj(qlms[:,-m])
    return -np.conj(qlms[:,-m])

def gG_l(pos, qlms, Qlms, is_center, Nbins, maxdist):
    """
    Spatial correlation of the qlms and the Qlms
    """
    assert len(pos) == len(qlms)
    assert len(qlms) == len(Qlms)
    assert len(is_center) == len(pos)
    maxsq = float(maxdist**2)
    hQ = np.zeros(Nbins)
    hq = np.zeros(Nbins)
    g = np.zeros(Nbins, int)
    code = """
    #pragma omp parallel for
    for(int i=0; i<Npos[0]; ++i)
    {
        if(!is_center(i)) 
            continue;
        for(int j=0; j<Npos[0]; ++j)
        {
            if(i==j) continue;
            double disq = 0.0;
            for(int dim=0; dim<3;++dim)
                disq += pow(pos(i,dim)-pos(j,dim), 2);
            if(disq>=(double)maxsq)
                continue;
            const int r = sqrt(disq/(double)maxsq)*Nbins;
            double pq = real(qlms(i,0)*conj(qlms(j,0)));
            for(int m=1; m<Nqlms[1]; ++m)
                pq += 2.0*real(qlms(i,m)*conj(qlms(j,m)));
            pq *= 4.0*M_PI/(2.0*(Nqlms[1]-1)+1);
            double pQ = real(Qlms(i,0)*conj(Qlms(j,0)));
            for(int m=1; m<NQlms[1]; ++m)
                pQ += 2.0*real(Qlms(i,m)*conj(Qlms(j,m)));
            pQ *= 4.0*M_PI/(2.0*(NQlms[1]-1)+1);
            #pragma omp critical
            {
                ++g(r);
                hq(r) += pq;
                hQ(r) += pQ;
            }
        }
    }
    """
    weave.inline(
        code,['qlms', 'Qlms', 'pos', 'maxsq', 'Nbins', 'hQ', 'hq', 'g','is_center'],
        type_converters =converters.blitz,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return hq, hQ, g

def periodic_gG_l(pos, L, qlms, Qlms, Nbins):
    """
    Spatial correlation of the qlms and the Qlms in a periodic box of size L
    """
    assert len(pos) == len(qlms)
    assert len(qlms) == len(Qlms)
    assert len(is_center) == len(pos)
    maxdist = L/2.0
    maxsq = float(maxdist**2)
    hQ = np.zeros(Nbins)
    hq = np.zeros(Nbins)
    g = np.zeros(Nbins, int)
    code = """
    #pragma omp parallel for
    for(int i=0; i<Npos[0]; ++i)
    {
        if(!is_center(i)) 
            continue;
        for(int j=0; j<Npos[0]; ++j)
        {
            if(i==j) continue;
            double disq = 0.0;
            for(int dim=0; dim<3;++dim)
                disq += pow(fmod(pos(i,dim)-pos(j,dim), L), 2);
            if(disq>=(double)maxsq)
                continue;
            const int r = sqrt(disq/(double)maxsq)*Nbins;
            double pq = real(qlms(i,0)*conj(qlms(j,0)));
            for(int m=1; m<Nqlms[1]; ++m)
                pq += 2.0*real(qlms(i,m)*conj(qlms(j,m)));
            pq *= 4.0*M_PI/(2.0*(Nqlms[1]-1)+1);
            double pQ = real(Qlms(i,0)*conj(Qlms(j,0)));
            for(int m=1; m<NQlms[1]; ++m)
                pQ += 2.0*real(Qlms(i,m)*conj(Qlms(j,m)));
            pQ *= 4.0*M_PI/(2.0*(NQlms[1]-1)+1);
            #pragma omp critical
            {
                ++g(r);
                hq(r) += pq;
                hQ(r) += pQ;
            }
        }
    }
    """
    weave.inline(
        code,['qlms', 'Qlms', 'pos', 'maxsq', 'Nbins', 'hQ', 'hq', 'g','is_center', 'L'],
        type_converters =converters.blitz,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return hq, hQ, g
            
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
