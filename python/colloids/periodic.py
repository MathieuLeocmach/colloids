import numpy as np
from colloids import vtk, statistics
import os, subprocess, shlex
try:
    from scipy import weave
    from scipy.weave import converters
except ImportError:
    try:
        import weave
        from weave import converters
    except ImportError:
        pass

def periodify(u, v, periods=None):
    """Given two arrays of points in a d-dimentional space with periodic boundary conditions, find the shortest vector between each pair"""
    assert u.shape == v.shape
    diff = np.array(v, float) - u
    if periods is None:
        return diff
    assert len(periods) == u.shape[-1]
    #ensures the largest coordinate is smaller than half a period
    half = 0.5*np.array(periods)
    pers = np.tile(periods, [len(diff),1])
    toolarge = diff > half
    while np.any(toolarge):
        diff[toolarge] -= pers[toolarge]
        toolarge = diff > half
    #same for negative coordinates
    toosmall = diff < -half
    while np.any(toosmall):
        diff[toosmall] += pers[toosmall]
        toosmall = diff < -half
    return diff
    

dist_code = """
inline double periodic_dist(const double &x, const double &y, const double &period)
{
    const double half = period*0.5;
    double d = x-y;
    while(d > half)
        d -= period;
    while(d < -half)
        d =+ period;
    return d;
}
"""

def get_displ(start, stop, L=203):
    """Compute displacement with periodic bounday conditions, assuming no particle moved further than L/2"""
    assert start.shape == stop.shape
    displ = np.zeros_like(start)
    code = """
    #pragma omp parallel for
    for(int i=0; i<Nstart[0]; ++i)
    {
        for(int j=0; j<Nstart[1];++j)
            displ(i,j) = periodic_dist(stop(i,j), start(i,j), L);
    }
    """
    weave.inline(
        code,['start', 'stop', 'displ', 'L'],
        type_converters =converters.blitz,
        support_code = dist_code,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return displ

def grv2vtk(pattern, T, radii=None):
    positions = np.asarray(
            [np.loadtxt(pattern%t + '.grv') for t in range(T)]
            )
    vel = np.gradient(positions)[0]
    for t in range(T):
            v= vtk.Polydata()
            v.points = positions[t]
            v.bonds = np.loadtxt(pattern%t + '.bonds', dtype=int)
            cloud = np.loadtxt(pattern%t+'.cloud')
            cloudcg = np.loadtxt(pattern%t+'.grv_space.cloud')
            border = np.zeros(v.bonds.shape[0], bool)
            border[np.where(
                ((v.points[v.bonds[:,1]] - v.points[v.bonds[:,0]])**2).sum(axis=1)>50*50)[0]
                   ]=True
            v.scalars = [
                    ('Q6', cloudcg[:,1]),
                    ('W4', cloudcg[:,4]),
                    ('w6', cloud[:,5]),
                    ]
            if radii != None:
                v.scalars.append(('radius', radii))
            v.vectors = [('vel', vel[t])]
            v.bondsScalars = [('border', border)]
            v.save(pattern%t+'.vtk')

def remove_drift(a, L=203):
    mean_displ = np.cumsum(
        [get_displ(a[i], a[i+1], L).mean(axis=0) for i in range(len(a)-1)],
        axis=0)
    for i, d in enumerate(mean_displ):
        a[i+1] -= d
    a = np.mod(a + 1.5 * L, L) - 0.5 * L

def msd(A, av=10, L=203):
    """
    Mean square displacement in a cubic box of size L
    If av is 0 (Default), the calculation will act greedily,
    averaging over all the avilabe intervals of a given length.
        Example : start=1 stop=4 av=0
            MSD[0] = 1
            MSD[1] = ( msd([1,2]) + msd([2,3]) + msd([3,4]))/3
            MSD[2] = ( msd([1,3]) + msd([2,4]) )/2
            MSD[3] = msd([1,4])
    If av>0, the average will be done over av time intervals starting
    from start, start+1,...,start+av-1
        Example : start=1 stop=4 av=2
            MSD[0] = 1
            MSD[1] = ( msd([1,2]) + msd([2,3]) )/2
            MSD[2] = ( msd([1,3]) + msd([2,4]) )/2
            MSD[3] = ( msd([1,4]) + msd([2,5]) )/2
"""
    msd = np.zeros(len(A)-av+1)
    if av==0:
        for t0, a in enumerate(A):
            for dt, b in enumerate(A[t0+1:]):
                #average is done over all trajectories and the 3 dimensions
                msd[dt+1] += (get_displ(a,b, L)**2).sum()
        msd /= A.shape[1]
        for dt in range(len(A)):
            msd[dt] /= len(A)-dt
        return msd
    else:
        for t0, a in enumerate(A[:av]):
            for dt, b in enumerate(A[t0+1:len(A)-av+t0+1]):
                msd[dt+1] += (get_displ(a,b, L)**2).sum()
        msd /= av * A.shape[1]
        return msd

def self_isf(positions, radius, av, L=203):
    """Self intermediate scattering function in a cubic periodic box"""
    A = np.exp(
        positions * (1j * np.pi * np.floor(L/radius) / L)
        )
    return statistics.time_correlation(A, av).mean(axis=-1).real


def ngp(A, av=10, L=203):
    """
    Non Gaussian parameter in a cubic box of size L
    """
    msd = np.zeros(len(A)-av+1)
    mqd = np.zeros(len(A)-av+1)
    if av==0:
        for t0, a in enumerate(A):
            for dt, b in enumerate(A[t0+1:]):
                #average is done over all trajectories and the 3 dimensions
                diff = (get_displ(a,b, L)**2).sum(axis=-1)
                msd[dt+1] += diff.sum()
                mqd[dt+1] += (diff**2).sum()
        for dt in range(len(A)):
            mqd[dt] *= len(A)-dt
        msd[0]=1
        return A.shape[1] * 3 * mqd / (5*msd**2) -1
    else:
        for t0, a in enumerate(A[:av]):
            for dt, b in enumerate(A[t0+1:len(A)-av+t0+1]):
                diff = (get_displ(a,b, L)**2).sum(axis=-1)
                msd[dt+1] += diff.sum()
                mqd[dt+1] += (diff**2).sum()
        msd[0]=1
        return av * A.shape[1] * 3 * mqd / (5*msd**2) -1
        
def overlap(pos, ov=0.3, av=10, L=203):
    """
    Overlap between configurations function of their time difference.
    If av is 0 (Default), the calculation will act greedily,
    averaging over all the avilabe intervals of a given length.
        Example : start=1 stop=4 av=0
            overlap[0] = 1
            overlap[1] = ( overlap([1,2]) + overlap([2,3]) + overlap([3,4]))/3
            overlap[2] = ( overlap([1,3]) + overlap([2,4]) )/2
            overlap[3] = overlap([1,4])
    If av>0, the average will be done over av time intervals starting
    from start, start+1,...,start+av-1
        Example : start=1 stop=4 av=2
            overlap[0] = 1
            overlap[1] = ( overlap([1,2]) + overlap([2,3]) )/2
            overlap[2] = ( overlap([1,3]) + overlap([2,4]) )/2
            overlap[3] = ( overlap([1,4]) + overlap([2,5]) )/2
    returns Nboverlap, Ntotal
"""
    overlap = np.zeros(len(pos)-av+1, int)
    if av==0:
        code ="""
        const double ovsq = ov*ov;
        #pragma omp parallel for schedule(dynamic)
        for (int t0=0; t0<Npos[0]-1; ++t0)
            for (int t1=t0+1; t1<Npos[0]; ++t1)
            {
                for(int i=0; i<Npos[1]; ++i)
                {
                    double disq = 0.0;
                    for(int dim=0; dim<Npos[2];++dim)
                        disq += pow(periodic_dist(pos(t0,i,dim), pos(t1,i,dim), L), 2);
                    if (disq < ovsq) overlap(t1-t0)++;
                }
            }
        """
        weave.inline(
            code,['pos', 'ov', 'L', 'overlap'],
            type_converters =converters.blitz,
            support_code = dist_code,
            extra_compile_args =['-O3 -fopenmp'],
            extra_link_args=['-lgomp'],
            verbose=2, compiler='gcc')
        overlap[0] = pos.shape[0] * pos.shape[1]
        return overlap, np.arange(1, len(pos))[::-1]*pos.shape[1]
    else:
        code ="""
        const double ovsq = ov*ov;
        const int av = Npos[0] - Noverlap[0]+1;
        #pragma omp parallel for
        for (int t0=0; t0<av; ++t0)
            for (int dt=1; dt<Npos[0]-av; ++dt)
            {
                for(int i=0; i<Npos[1]; ++i)
                {
                    double disq = 0.0;
                    for(int dim=0; dim<Npos[2];++dim)
                        disq += pow(periodic_dist(pos(t0,i,dim), pos(t0+dt,i,dim), L), 2);
                    if (disq < ovsq) overlap(dt)++;
                }
            }
        """
        weave.inline(
            code,['pos', 'ov', 'L', 'overlap'],
            type_converters =converters.blitz,
            support_code = dist_code,
            extra_compile_args =['-O3 -fopenmp'],
            extra_link_args=['-lgomp'],
            verbose=2, compiler='gcc')
        overlap[0] = av * pos.shape[1]
        return overlap, np.ones(len(overlap), int) * av * pos.shape[1]
        
def get_N_ngbs(pos, radii, L=203, N=12, maxdist=8.0):
    assert len(pos) == len(radii)
    #initialize the geometry of each particle
    neighbours = -np.ones([len(pos), N], int)
    code = """
    const double maxdistsq = maxdist*maxdist;
    //look for nearby particles
    #pragma omp parallel for
    for(int p=0; p<Npos[0]; ++p)
    {
        std::multimap<double, int> ngbbydist;
        for(int q=0; q<Npos[0]; ++q)
        {
            if (p==q) continue;
            double disq = 0.0;
            for(int dim=0; dim<3;++dim)
                disq += pow(periodic_dist(pos(p,dim), pos(q,dim), L), 2);
            disq /= pow(radii(p) + radii(q), 2);
            if (disq>maxdistsq) continue;
            ngbbydist.insert(std::make_pair(disq, q));
        }
        std::multimap<double, int>::const_iterator it = ngbbydist.begin();
        for(int i=0; i<Nneighbours[1] && it!=ngbbydist.end(); ++i)
            neighbours(p, i) = (it++)->second;
    }
    """
    weave.inline(
        code,['pos', 'radii', 'maxdist', 'L', 'neighbours'],
        type_converters =converters.blitz,
        support_code = dist_code,
        headers = ['<map>', '<list>'],
        extra_compile_args =['-O3 -fopenmp -march=native'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return neighbours

def loadXbonds(fname, Q6, thr=0.25):
    """Load the bonds linking two MRCO particles"""
    #load all bonds
    bonds = np.loadtxt(fname, dtype=int)
    return bonds[np.where(Q6[bonds].min(axis=1)>thr)]
    
def get_rdf(p, pos, Nbins, L=203.0, maxdist=None):
    if maxdist is None:
        maxdist = L/2.0
    imaxsq = 1.0/maxdist**2
    g = np.zeros(Nbins, int)
    code = """
    #pragma omp parallel for
    for(int j=0; j<Npos[0]; ++j)
    {
        double disq = 0.0;
            for(int dim=0; dim<3;++dim)
                disq += pow(periodic_dist(p, pos(j,dim), L), 2);
        if(disq*imaxsq>=1.0)
            continue;
        const int r = sqrt(disq*imaxsq)*Ng[0];
        #pragma omp atomic
        ++g(r);
    }
    """
    weave.inline(
        code,['p', 'pos', 'imaxsq', 'L', 'g'],
        type_converters =converters.blitz,
        support_code = dist_code,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return g
    
def get_mono_rdf(pos, Nbins, L, maxdist=None):
    if maxdist is None:
        maxdist = L/2.0
    imaxsq = 1.0/(maxdist**2)
    g = np.zeros([len(pos), Nbins], int)
    code = """
    #pragma omp parallel for
    for(int i=0; i<Npos[0]; ++i)
    {
        for(int j=i+1; j<Npos[0]; ++j)
        {
            double disq = 0.0;
            for(int dim=0; dim<3;++dim)
                disq += pow(periodic_dist(pos(i,dim), pos(j,dim), L), 2);
            if(disq*imaxsq>=1)
                continue;
            const int r = sqrt(disq*imaxsq)*Ng[1];
            #pragma omp atomic
            ++g(i, r);
            #pragma omp atomic
            ++g(j, r);
        }
    }
    """
    weave.inline(
        code,['pos', 'imaxsq', 'L', 'g'],
        type_converters =converters.blitz,
        support_code = dist_code,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return g
    
def get_binary_rdf(pos, Nbins, L, sep=800, maxdist=None):
    if maxdist is None:
        maxdist = L/2.0
    imaxsq = 1.0/(maxdist**2)
    g = np.zeros([len(pos), 2, Nbins], int)
    code = """
    #pragma omp parallel for
    for(int i=0; i<Npos[0]; ++i)
    {
        for(int j=0; j<sep; ++j)
        {
            if(i==j)
                continue;
            double disq = 0.0;
            for(int dim=0; dim<3;++dim)
                disq += pow(periodic_dist(pos(i,dim), pos(j,dim), L), 2);
            if(disq*imaxsq>=1)
                continue;
            const int r = sqrt(disq*imaxsq)*Nbins;
            #pragma omp atomic
            ++g(i, 0, r);
        }
        for(int j=sep; j<Npos[0]; ++j)
        {
            if(i==j)
                continue;
            double disq = 0.0;
            for(int dim=0; dim<3;++dim)
                disq += pow(periodic_dist(pos(i,dim), pos(j,dim), L), 2);
            if(disq*imaxsq>=1)
                continue;
            const int r = sqrt(disq*imaxsq)*Nbins;
            #pragma omp atomic
            ++g(i, 1, r);
        }
    }
    """
    weave.inline(
        code,['pos', 'imaxsq', 'L', 'sep', 'g', 'Nbins'],
        type_converters =converters.blitz,
        support_code = dist_code,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return g
    
def get_space_correl(pos, field, Nbins, L):
    maxdist = L/2.0
    imaxsq = 1.0/(maxdist**2)
    g = np.zeros(Nbins, int)
    h = np.zeros([Nbins])
    code = """
    #pragma omp parallel for
    for(int i=0; i<Npos[0]; ++i)
    {
        for(int j=0; j<Npos[0]; ++j)
        {
            if(i==j)
                continue;
            double disq = 0.0;
            for(int dim=0; dim<3;++dim)
                disq += pow(periodic_dist(pos(i,dim), pos(j,dim), L), 2);
            if(disq*imaxsq>=1)
                continue;
            const int r = sqrt(disq*imaxsq)*Ng[0];
            const double v = field(i)*field(j);
            #pragma omp critical
            {
                ++g(r);
                h(r) += v;
            }
        }
    }
    """
    weave.inline(
        code,['pos', 'imaxsq', 'L', 'field', 'g', 'h'],
        type_converters =converters.blitz,
        support_code = dist_code,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return h, g
    
#def get_space_correl(pos, field, Nbins, L):
 #   h, g = get_space_correls(pos, field[:,None], Nbins, L)
  #  return h[:,0], g
    
def get_space_correl_binary(pos, field, Nbins, L, sep=800):
    maxdist = L/2.0
    imaxsq = 1.0/(maxdist**2)
    g = np.zeros([4,Nbins], int)
    h = np.zeros([4,Nbins])
    code = """
    #pragma omp parallel for
    for(int i=0; i<Npos[0]; ++i)
    {
        for(int j=0; j<Npos[0]; ++j)
        {
            if(i==j)
                continue;
            double disq = 0.0;
            for(int dim=0; dim<3;++dim)
                disq += pow(periodic_dist(pos(i,dim), pos(j,dim), L), 2);
            if(disq*imaxsq>=1)
                continue;
            const int r = sqrt(disq*imaxsq)*Ng[1];
            const double v = field(i)*field(j);
            if(i<sep)
            {
                if(j<sep)
                #pragma omp critical
                {
                    ++g(0,r);
                    h(0,r) += v;
                }
                else
                #pragma omp critical
                {
                    ++g(1,r);
                    h(1,r) += v;
                }
            }
            else
            {
                if(j<sep)
                #pragma omp critical
                {
                    ++g(2,r);
                    h(2,r) += v;
                }
                else
                #pragma omp critical
                {
                    ++g(3,r);
                    h(3,r) += v;
                }
            }
        }
    }
    """
    weave.inline(
        code,['pos', 'imaxsq', 'L', 'field', 'sep', 'g', 'h'],
        type_converters =converters.blitz,
        support_code = dist_code,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return h, g

code_kvec = """
inline int generate_kvec(const int nk, blitz::Array<int,2> &kvectors, const int maxNvec=30)
{
    assert(maxNvec%3==0);
    kvectors = 0;
    int nvec = 0;
    for(int i=0; i<nk && nvec<maxNvec; ++i)
        for(int j=0; j*j<nk*nk-i*i && nvec<maxNvec; ++j)
        {
            const int k = sqrt(nk*nk-i*i-j*j);
            if(k != floor(k)) continue;
            kvectors(nvec++, blitz::Range::all()) = i, j, k;
            kvectors(nvec++, blitz::Range::all()) = j, k, i;
            kvectors(nvec++, blitz::Range::all()) = k, i, j;
        }
    return nvec;
}
"""

def get_Sq(pos, Nbins, L=203.0, maxNvec=30, field=None):
    if field is None:
        field = np.ones(len(pos))
    assert len(field) == len(pos)
    Sq = np.zeros(Nbins)
    #cache the complex exponential caculations
    cache = np.ones([Nbins]+list(pos.shape), np.complex128)
    code = """
    //cache the complex exponential caculations
    #pragma omp parallel for
    for(int i=0; i<Npos[0]; ++i)
        for(int j=0; j<Npos[1]; ++j)
        {
            std::complex<double> e = std::polar(1.0, 2*M_PI/L * pos(i,j));
            for(int nk=1; nk<NSq[0]; ++nk)
                cache(nk,i,j) = cache(nk-1,i,j) * e;
        }
    
    #pragma omp parallel for schedule(dynamic)
    for(int nk=1; nk<NSq[0]; ++nk)
    {
        //generate the k-vectors
        blitz::Array<int,2> kvectors(maxNvec, 3);
        const int nvec = generate_kvec(nk, kvectors, maxNvec);
        
        for(int k=0; k<nvec; ++k)
        {
            std::complex<double> sum_rho(0.0,0.0);
            for(int i=0; i<Npos[0]; ++i)
            {
                std::complex<double> prod(field(i));
                for(int j=0; j<Npos[1]; ++j)
                {
                    //random access to cache
                    const int kv = kvectors(k,j);
                    prod *= cache(kv, i,j);
                }
                sum_rho += prod;
            }
            Sq(nk) += std::norm(sum_rho);
        }
        Sq(nk) /= nvec;
    }
    """
    weave.inline(
        code,['pos', 'field', 'Sq', 'cache', 'L', 'maxNvec'],
        support_code = code_kvec,
        type_converters =converters.blitz,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return Sq/len(pos)
    
#same in a rectangular box
code_kvec_rect = """
inline int generate_kvec_rect(const int nk, blitz::Array<double,2> &kvectors, const blitz::Array<double,1> &Lsq, const int maxNvec=30)
{
    kvectors = 0;
    int nvec = 0;
    for(int i=0; i<nk+1 && nvec<maxNvec; ++i)
        for(int j=0; Lsq(0)*j*j<Lsq(1)*(nk*nk-i*i)+1 && nvec<maxNvec; ++j)
        {
            const int ksqmax = Lsq(2) * ((nk*nk-i*i)/Lsq(0) - j*j/Lsq(1));
            const int ksqmin = Lsq(2) * (((nk-1)*(nk-1)-i*i)/Lsq(0) - j*j/Lsq(1));
            for (int k=sqrt(std::max(0, ksqmin)); k<=sqrt(ksqmax) && nvec<maxNvec; ++k)
                if (((nk-1)*(nk-1) - i*i)/Lsq(0) < j*j/Lsq(1) + k*k/Lsq(2))
                    kvectors(nvec++, blitz::Range::all()) = i, j, k;
        }
    return nvec;
}
"""

def rectangular_Sq(positions, Nbins, Ls=[203.0]*3, maxNvec=30, field=None):
    #sort xyz by decreasing box size
    pos = positions[:, np.argsort(Ls)[::-1]]
    dims = np.sort(Ls)[::-1]
    Lsq = dims.astype(float)**2
    if field is None:
        field = np.ones(len(pos))
    assert len(field) == len(pos)
    assert maxNvec%3==0
    Sq = np.zeros(Nbins)
    #cache the complex exponential caculations
    cache = np.ones([Nbins]+list(pos.shape), np.complex128)
    code = """
    //cache the complex exponential caculations
    #pragma omp parallel for
    for(int i=0; i<Npos[0]; ++i)
        for(int j=0; j<Npos[1]; ++j)
        {
            std::complex<double> e = std::polar(1.0, 2*M_PI/dims(j) * pos(i,j));
            for(int nk=1; nk<NSq[0]; ++nk)
                cache(nk,i,j) = cache(nk-1,i,j) * e;
        }
    
    #pragma omp parallel for schedule(dynamic)
    for(int nk=1; nk<NSq[0]; ++nk)
    {
        //generate the k-vectors
        blitz::Array<double,2> kvectors(maxNvec, 3);
        const int nvec = generate_kvec_rect(nk, kvectors, Lsq, maxNvec);
        //const int nvec = generate_kvec(nk, kvectors, maxNvec);
        
        for(int k=0; k<nvec; ++k)
        {
            std::complex<double> sum_rho(0.0,0.0);
            for(int i=0; i<Npos[0]; ++i)
            {
                std::complex<double> prod(field(i));
                for(int j=0; j<Npos[1]; ++j)
                {
                    //random access to cache
                    const int kv = kvectors(k,j);
                    prod *= cache(kv, i,j);
                }
                sum_rho += prod;
            }
            Sq(nk) += std::norm(sum_rho);
        }
        Sq(nk) /= nvec;
    }
    """
    weave.inline(
        code,['pos', 'field', 'Sq', 'dims', 'cache', 'Lsq', 'maxNvec'],
        support_code = code_kvec_rect,
        type_converters =converters.blitz,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return Sq/len(pos)
