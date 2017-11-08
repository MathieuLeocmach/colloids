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
