#depreciated functions using scipy.weave

import numpy as np
from scipy.special import sph_harm
try:
    from scipy import weave
    from scipy.weave import converters
except ImportError:
    try:
        import weave
        from weave import converters
    except ImportError:
        pass
import numexpr
from colloids import periodic

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
    
def periodic_qlm(pos, ngbs, period, l=6, weights=None):
    assert len(pos) == len(ngbs)
    if weights is None:
        weights = np.ones(ngbs.shape, float)/ngbs.shape[1]
    assert weights.shape == ngbs.shape
    qlm = np.zeros([len(pos), l+1], np.complex128)
    code = """
    #pragma omp parallel for
    for(int i=0; i<Nngbs[0]; ++i)
    {
        double cart[3] = {0, 0, 0};
        double sph[3] = {0, 0, 0};
        //int nb = 0;
        double total_weight = 0.0;
        for(int j=0; j<Nngbs[1]; ++j)
        {
            int q = ngbs(i,j);
            if(q<0 || q>=Npos[0])
                continue;
            //nb++;
            double weight = weights(i,j);
            total_weight += weight;
            for(int d=0; d<3; ++d)
                cart[d] = periodic_dist(pos(i,d), pos(q, d), period);
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
                qlm(i,m) += weight * boost::math::spherical_harmonic(Nqlm[1]-1, m, sph[1], sph[2]);
        }
        if(total_weight>0)
            for(int m=0; m<Nqlm[1]; ++m)
                qlm(i,m) /= total_weight;
    }
    """
    weave.inline(
        code,['pos', 'ngbs', 'period', 'weights', 'qlm'],
        type_converters =converters.blitz,
        headers = ['<boost/math/special_functions/spherical_harmonic.hpp>'],
        support_code = periodic.dist_code,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return qlm
    
def steinhardt_g_l(pos, bonds, is_center, Nbins, maxdist, l=6):
    """
    Spatial correlation of the bond's spherical harmonics
    """
    assert len(is_center) == len(pos)
    maxsq = float(maxdist**2)
    hq = np.zeros(Nbins)
    g = np.zeros(Nbins, int)
    qlms = np.zeros([len(bonds), l+1], np.complex128)
    bpos = np.zeros([len(bonds), 3])
    code = """
    //position and spherical harmonics of each bond
    #pragma omp parallel for
    for(int b=0; b<Nbonds[0]; ++b)
    {
        int i = bonds(b,0), j = bonds(b,1);
        blitz::Array<double,1> cart(pos(i, blitz::Range::all()) - pos(j, blitz::Range::all()));
        bpos(b, blitz::Range::all()) = pos(j, blitz::Range::all()) + 0.5 * cart;
        double sph[3] = {0, 0, 0};
        sph[0] = sqrt(blitz::sum(blitz::pow(cart, 2)));
        if(abs(cart(2))==sph[0] || sph[0]*sph[0]+1.0 == 1.0)
        {
            sph[1] = 0;
            sph[2] = 0;
        }
        else
        {
            sph[1] = acos(cart(2)/sph[0]);
            sph[2] = atan2(cart(1), cart(0));
            if(sph[2]<0)
                sph[2] += 2.0*M_PI;
        }
        for(int m=0; m<Nqlms[1]; ++m)
            qlms(b,m) = boost::math::spherical_harmonic(Nqlms[1]-1, m, sph[1], sph[2]);
    }
    #pragma omp parallel for
    for(int b=0; b<Nbonds[0]; ++b)
    {
        int i = bonds(b,0), j = bonds(b,1);
        if(!is_center(i) || !is_center(j))
            continue;
        for(int c=0; c<Nbonds[0]; ++c)
        {
            const double disq = blitz::sum(blitz::pow(bpos(b, blitz::Range::all()) - bpos(c, blitz::Range::all()),2));
            if(disq>=(double)maxsq)
                continue;
            const int r = sqrt(disq/(double)maxsq)*Nbins;
            double pq = real(qlms(b,0)*conj(qlms(c,0)));
            for(int m=1; m<Nqlms[1]; ++m)
                pq += 2.0*real(qlms(b,m)*conj(qlms(c,m)));
            pq *= 4.0*M_PI/(2.0*(Nqlms[1]-1)+1);
            #pragma omp critical
            {
                ++g(r);
                hq(r) += pq;
            }
        }
    }
    """
    weave.inline(
        code,['qlms', 'pos', 'bonds', 'bpos', 'maxsq', 'Nbins', 'hq', 'g','is_center'],
        type_converters =converters.blitz,
        headers=['<boost/math/special_functions/spherical_harmonic.hpp>'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return hq, g
