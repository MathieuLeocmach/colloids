import numpy as np
from colloids import vtk, statistics
import os, subprocess, shlex

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
        msd /= A.shape[1] * A.shape[2]
        for dt in range(len(A)):
            msd[dt] /= len(A)-dt
        return msd
    else:
        for t0, a in enumerate(A[:av]):
            for dt, b in enumerate(A[t0+1:len(A)-av+t0+1]):
                msd[dt+1] += (get_displ(a,b, L)**2).sum()
        msd /= av * A.shape[1] * A.shape[2]
        return msd

def self_isf(positions, radius, av, L=203):
    """Self intermediate scattering function in a cubic periodic box"""
    A = np.exp(
        positions * (1j * np.pi * np.floor(L/radius) / L)
        )
    return statistics.time_correlation(A, av)


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
                diff = get_displ(a,b, L)**2
                msd[dt+1] += diff.sum()
                mqd[dt+1] += (diff.sum(axis=-1)**2).sum()
        for dt in range(len(A)):
            mqd[dt] *= len(A)-dt
        return A.shape[1] * A.shape[2] * 3 * mqd / (5*msd**2) -1
    else:
        for t0, a in enumerate(A[:av]):
            for dt, b in enumerate(A[t0+1:len(A)-av+t0+1]):
                diff = get_displ(a,b, L)**2
                msd[dt+1] += diff.sum()
                mqd[dt+1] += (diff.sum(axis=-1)**2).sum()
        return av * A.shape[1] * A.shape[2] * 3 * mqd / (5*msd**2) -1

def loadXbonds(fname, Q6, thr=0.25):
    """Load the bonds linking two MRCO particles"""
    #load all bonds
    bonds = np.loadtxt(fname, dtype=int)
    return bonds[np.where(Q6[bonds].min(axis=1)>thr)]

