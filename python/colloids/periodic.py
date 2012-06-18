import numpy as np
from colloids import vtk, statistics
import os, subprocess, shlex
from scipy import weave
from scipy.weave import converters

def get_displ(start, stop, L=203):
    """Compute displacement with periodic bounday conditions, assuming no particle moved further than L/2"""
    return np.mod(stop - start + 1.5 * L, L) - 0.5 * L

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
        const double disq = blitz::sum(blitz::pow(blitz::fmod(
            blitz::abs(pos(j,blitz::Range::all())-p), 
            0.5*L
            ), 2));
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
                disq += pow(fmod(fabs(pos(i,dim)-pos(j,dim)), 0.5*L), 2);
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
                disq += pow(fmod(fabs(pos(i,dim)-pos(j,dim)), 0.5*L), 2);
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
                disq += pow(fmod(fabs(pos(i,dim)-pos(j,dim)), 0.5*L), 2);
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
                disq += pow(fmod(fabs(pos(i,dim)-pos(j,dim)), 0.5*L), 2);
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
                disq += pow(fmod(fabs(pos(i,dim)-pos(j,dim)), 0.5*L), 2);
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
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return h, g

