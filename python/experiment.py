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
from __future__ import with_statement #for python 2.5, useless in 2.6
import scipy as sp
import numpy as np
import scipy.constants as const
from scipy.ndimage.filters import gaussian_filter1d
from scipy.ndimage.morphology import grey_dilation
import os, os.path, subprocess
import re, string
from math import exp

def noNaN(x):
    if x == 'NaN':
        return np.nan
    else:
        return float(x)

class Experiment:
    """A serie of coordinates and it's derived data"""

    def __init__(self, trajPath, T=310, pixel_size=297e-9):
        """Constructor. If the file.traj exist, read it. Else cut and link."""
        self.path,self.trajfile = os.path.split(trajPath)
        self.path = os.path.abspath(self.path)
        self.unitLength = pixel_size
        self.T = T
        if os.path.exists(os.path.join(self.path,self.trajfile)):
            self.read_traj()
        else:
            self.head = os.path.splitext(self.trajfile)[0]
            #find size, token and number of digits out of the .dat filenames
            self.digits = 0
            mask = re.compile(self.head+'(.*?)([0-9]+)\.dat')
            files = filter(mask.match,os.listdir(self.path))
            self.offset = 0
            self.size = len(files)
            for m in map(mask.match,files):
                if len(m.group(2)) > self.digits:
                    self.digits = len(m.group(2))
                    self.token = m.group(1)
            self.cut()
            #guess radius from g(r)
            self.radius = self.rdf_radius()
            #guess time interval from file name
            self.dt = 0.0
            m = re.search('([0-9]+)min',self.head)
            if m:
                self.dt += 60*int(m.group(1))
            m = re.search('([0-9]+)s([0-9]*)',self.head)
            if m:
                self.dt += float(m.group(1)+'.'+m.group(2))
            self.link()

    def read_traj(self):
        with open(os.path.join(self.path,self.trajfile),'r') as f:
            self.radius,self.dt = [
                float(s) for s in re.split("\t",f.readline()[:-1])
                ]
            pattern, = re.split("\t",f.readline()[:-1])
            self.token, = re.split("\t",f.readline()[:-1])
            self.offset,self.size = [
                int(s) for s in re.split("\t",f.readline()[:-1])
                ]
            self.read_pattern(pattern)

    def read_pattern(self,pattern):
        m = re.match('(.*)'+self.token+'([0-9]*)',os.path.splitext(pattern)[0])
        self.head = m.group(1)
        self.digits = len(m.group(2))

    def cut(self):
        """Remove centers closer than 3 pixels in each frame."""
        subprocess.check_call(
                map(str,['cutter',self.get_format_string()%0, 
                   self.token, 5, self.offset, self.size, 0.3])
                )

    def link(self):
        """Link trajectories. Export a file.traj"""
        actual = os.getcwd()
        os.chdir(self.path)
        subprocess.check_call(map(str,
                 ['linker', 
                 self.get_format_string(absPath=False)%0, self.token,
                 self.radius, self.dt,
                 self.offset,self.size])
            )
        os.chdir(actual)

    def get_range(self):
        return range(self.offset,self.offset+self.size)

    def get_format_string(self,postfix='',ext='dat',absPath=True):
        format_string = self.head + postfix + self.token + \
                        '%0' + str(self.digits) + 'd.' + ext
        if absPath:
            format_string = os.path.join(self.path,format_string)
        return format_string

    def mean_Nb(self):
        """Calculate the time averaged number of particles in a frame"""
        if not hasattr(self,'__mean_Nb'):
            nb = 0L
            for t,name in enum(self):
                with open(name,'r') as f:
                    nb += int(re.split("\t",f.readline())[1])
            self.__mean_Nb = float(nb) / self.size
        return self.__mean_Nb

    def mean_V(self):
        """Calculate the time averaged volume"""
        if not hasattr(self,'__mean_V'):
            V=0
            for t,name in enum(self):
                V += np.ptp(
                    np.loadtxt(name, delimiter='\t', skiprows=2),
                    axis=0).prod()
            self.__mean_V = V / self.size
        return self.__mean_V

    def get_VF(self):
        """Calculate an estimte of the volume fraction"""
        return 4 * np.pi * self.mean_Nb() * self.radius**3 / (3 * self.mean_V())

    def rdf_radius(self,force=False):
        """Return the place of the largest (first) peak of the g(r),
        in pixel unit"""
        name = os.path.join(self.path, self.head + '_total.rdf')
        if force or not os.path.exists(name):
            subprocess.check_call(map(str,
                  ['totalRdf',self.get_format_string()%0, self.token,
                   200, 15])
                )
        r,g = np.loadtxt(name, unpack=True)
        return r[np.argmax(g)]
        gm = g.argmax()
        return np.average(r[gm-1:gm+2],weights=g[gm-1:gm+2])

    def get_Nb_density(self, averaged=True):
        nbs = np.empty((size))
        Vs = np.empty((size))
        for t,fname in enum(self):
            coords = np.loadtxt(fname,delimiter='\t', skiprows=2)
            nbs[t-self.offset] = len(coords)
            Vs[t-self.offset] = np.ptp(coords, axis=0).prod()
        if averaged:
            return (nbs/Vs).mean()
        else:
            return (nbs, Vs)
    
    def get_tau(self,thr=exp(-1),force=False):
        """Get the Self-ISF decay time"""
        if not hasattr(self,'__tau'):
            name = os.path.join(self.path, self.head + '.isf')
            if force or not os.path.exists(name):
                subprocess.check_call(
                   ['dynamics',
                   os.path.join(self.path,self.head + '.traj'),
                   '1']
                   )
            t,isf = np.loadtxt(name, delimiter="\t",usecols=[0,4],unpack=True)
            below = np.argwhere(isf<thr)
            if len(below)>0:
                self.__tau = t[below.min()]
            else:
                self.__tau = t[-1]
        return self.__tau

    def boo(self,total=True):
        """Calculate bond orientational order"""
        if total:
            subprocess.check_call([
                   'totalBoo',
                   os.path.join(self.path,self.head + '.traj')
                   ])
        else:
            subprocess.check_call([
                   'boo',
                   os.path.join(self.path,self.head + '.traj')
                   ])
    def ageing(self,intervals):
        """Calculate isf for various time sub-intervals"""
        subprocess.check_call([
                   'ageing',
                   os.path.join(self.path,self.head + '.traj')
                   ]+map(str,np.asarray(intervals).ravel())
                              )
    def get_intervals(self, step=10, average=10):
        """get a few typical intervals"""
        return [(t,self.offset+self.size-(average+1), average) \
                     for t in self.get_range()[:-50:step]]

    def auto_ageing(self, step=10, average=10):
        """calculate isf for a few typical intervals"""
        self.ageing(self.get_intervals(step=10, average=10))

    def get_interval_file(self, inter,ext='isf'):
        return os.path.join(
                self.path,
                self.head + ('_%ifrom_%ito_%iav.' % inter)+ext
                )

    def get_plot_string(self, title, inter):
        return '"%s" using ($1/%f):5 title "%s"' \
               %(self.get_interval_file(inter),
                 br(
                         sigma = 2*self.radius*self.unitLength,
                         T=self.T
                         ),
                 title)
        
        
    def g6(self, Nbins=200, nbDiameters=4.5, force=False):
        """
        Calculate g6 and g for each time step and return the time average
        output is (r,g6,g)
        """
        if not force:
            for t, name in enum(self,ext='g6'):
                if not os.path.exists(name):
                    force=True
                    break
        if force:
            for t, name in enum(self):
                subprocess.check_call(map(str,
                    ['g6', name, self.radius, Nbins, nbDiameters]))
        tot = np.zeros((Nbins,3))
        for t, name in enum(self,ext='g6'):
            tot += np.nan_to_num(
                np.loadtxt(name,converters={0:noNaN,1:noNaN,2:noNaN})
                )
        tot /= self.size
        np.savetxt(
            os.path.join(self.path,self.head + '_total.g6'),
            tot,
            fmt='%f',
            delimiter='\t'
            )
        return tot
    
    def g6_envelope(self, **args):
        """output the envelope of g6"""
        r,g6,g = np.transpose(self.g6(**args))
        #smooth g6
        sg6 = gaussian_filter1d(np.copy(g6),1)
        maxima = np.flatnonzero(
                np.where(np.exp(
                        sg6 - grey_dilation(sg6,size=[3])
                        )>0.9999, 1, 0)
                )
        maxima = [m for m in maxima if g6[m]>0]
        envelope = np.column_stack((r[maxima],g6[maxima]))
        np.savetxt(
            os.path.join(self.path,self.head + '_total_env.g6'),
            envelope,
            fmt='%f',
            delimiter='\t'
            )
        return envelope

class Txp:
    """Implementig time algorithms in python"""
    
    def __init__(self, xp, start=None, size=None):
        self.xp = xp
        if not start or start < self.xp.offset:
            start = self.xp.offset
        if not size or start+size > self.xp.offset+self.xp.size:
            size = self.xp.size + self.xp.offset - start
        self.trajs = self.read_trajs(start, size)
        self.positions = self.read_pos(start, size)
        self.remove_drift()

    def read_trajs(self, start, size):
        """
            Reads the linked trajectories from the .traj file
            Retrieves only trajectories spanning [start,start+size[
            """
        trajs = []
        with open(os.path.join(self.xp.path,self.xp.trajfile),'r') as f:
            for l, line in enumerate(f):
                if l==3:
                    break
            for line in f:
                t0 = int(line[:-1])
                if t0 > start:
                    f.next()
                else:
                    pos = string.split(f.next()[:-1],'\t')[start-t0:]
                    if len(pos)>=size:
                        trajs.append(map(int, pos[:size]))
        return np.asarray(trajs)

    def read_pos(self, start, size):
        """Reads the usefull positions from the .dat files"""
        pos = np.empty((self.trajs.shape[1],self.trajs.shape[0],3))
        for t, fname in enum(self.xp):
            if t<start or t>= start+size:
                continue
            raw_pos = np.loadtxt(fname,delimiter='\t',skiprows=2)
            pos[t-start] = raw_pos[self.trajs[:,t-start]]
        return pos

    def remove_drift(self):
        """Remove the average drift between time steps"""
        drift = np.cumsum(np.diff(self.positions,axis=0).mean(axis=1),axis=0)
        sw = np.swapaxes(self.positions[1:],0,1)
        sw -= drift

    def msd(self,start,stop,av):
        """
        Mean square displacement
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
        A = self.positions[start:stop+av]
        msd = np.zeros((stop-start))
        if av==0:
            for t0, a in enumerate(A):
                for dt, b in enumerate(A[t0+1:]):
                    #average is done over all trajectories and the 3 dimensions
                    msd[dt+1] += ((b-a)**2).sum()
            msd /= A.shape[1] * A.shape[2] * (self.xp.radius*2)**2
            for dt, n in enumerate(range(stop-start,0,-1)):
                msd[dt+1] /= n
            return msd
        else:
            for t0, a in enumerate(A[:av]):
                for dt, b in enumerate(A[t0+1:-av+t0]):
                    msd[dt+1] += ((b-a)**2).sum()
            msd /= av * A.shape[1] * A.shape[2]  * (self.xp.radius*2)**2
            return msd

    def export_msd(self,start,stop,av):
        np.savetxt(
            self.xp.get_interval_file((start,stop,av),ext='msd'),
            np.column_stack((
                np.arange(0,stop-start+1)*self.xp.dt,
                self.msd(start,stop,av),
                )),
            fmt='%f',
            delimiter='\t'
            )

    def self_isf(self,start,stop,av):
        """
        Self intermediate scattering function
        If av is 0 (Default), the calculation will act greedily,
        averaging over all the avilabe intervals of a given length.
            Example : start=1 stop=4 av=0
                ISF[0] = 1
                ISF[1] = ( isf([1,2]) + isf([2,3]) + isf([3,4]))/3
                ISF[2] = ( isf([1,3]) + isf([2,4]) )/2
                ISF[3] = isf([1,4])
	If av>0, the average will be done over av time intervals starting
	from start, start+1,...,start+av-1
            Example : start=1 stop=4 av=2
                ISF[0] = 1
                ISF[1] = ( isf([1,2]) + isf([2,3]) )/2
                ISF[2] = ( isf([1,3]) + isf([2,4]) )/2
                ISF[3] = ( isf([1,4]) + isf([2,5]) )/2
	"""
        A = np.exp(
            self.positions[start:stop+av] * (1j * np.pi / self.xp.radius)
            )
        isf = np.zeros((stop-start+1))
        if av==0:
            for t0, a in enumerate(A):
                for dt, b in enumerate(A[t0+1:]):
                    #average is done over all trajectories and the 3 dimensions
                    isf[dt+1] += np.real(a.conj()*b).sum()
            isf /= A.shape[1] * A.shape[2]
            isf[0]=1
            for dt, n in enumerate(range(stop-start,0,-1)):
                isf[dt+1] /= n
            return isf
        else:
            for t0, a in enumerate(A[:av]):
                for dt, b in enumerate(A[t0+1:-av+t0]):
                    isf[dt+1] += np.real(a.conj()*b).sum()
            isf /= av * A.shape[1] * A.shape[2]
            isf[0]=1
            return isf
        
    def export_self_isf(self,start,stop,av):
        np.savetxt(
            self.xp.get_interval_file((start,stop,av)),
            np.column_stack((
                np.arange(0,stop-start+1)*self.xp.dt,
                self.self_isf(start,stop,av),
                )),
            fmt='%f',
            delimiter='\t'
            )

    def export_dynamics(self,start,stop,av):
        self.export_msd(start,stop,av)
        self.export_self_isf(start,stop,av)
        
    
def enum(xp,postfix='',ext='dat',absPath=True):
    """Generator of couples (time, filename)"""
    format_string = xp.get_format_string(postfix,ext,absPath)
    for t in xp.get_range():
        yield t, (format_string % t)

def br(sigma, T=310, eta=2.22e-3):
        """Brownian time for a particle of diameter sigma (in meters)"""
        return 3 * const.pi * eta * (sigma**3) / (4 * const.k * T)
