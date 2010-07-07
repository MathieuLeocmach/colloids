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
#import scipy.constants as const
from scipy.ndimage.filters import gaussian_filter1d
from scipy.ndimage.morphology import grey_dilation
import os, os.path, subprocess
import re, string, math
from math import exp
from colloids import vtk

def noNaN(x):
    if x == 'NaN':
        return np.nan
    else:
        return float(x)

class Experiment:
    """A serie of coordinates and it's derived data"""

    def __init__(self, trajPath, T=310, pixel_size=297e-9):
        """Constructor. If the file.traj exist, read it. Else launch linkboo."""
        self.path,self.trajfile = os.path.split(trajPath)
        self.path = os.path.abspath(self.path)
        self.unitLength = pixel_size
        self.T = T
        if not os.path.exists(os.path.join(self.path,self.trajfile)):
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
            #guess time interval from file name
            self.dt = 0.0
            m = re.search('([0-9]+)min',self.head)
            if m:
                self.dt += 60*int(m.group(1))
            m = re.search('([0-9]+)s([0-9]*)',self.head)
            if m:
                self.dt += float(m.group(1)+'.'+m.group(2))
            self.linkboo()
        self.read_traj()

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

    def linkboo(self):
        """calculate total g(r), radius, BOO for each time step and link trajectories."""
        actual = os.getcwd()
        os.chdir(self.path)
        subprocess.check_call(map(str,
                 ['linkboo',
                  self.get_format_string(absPath=False)%0,
                  self.token,
                  self.dt,self.size,
                  self.offset])
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

    def enum(self,postfix='',ext='dat',absPath=True):
        """Generator of couples (time, filename)"""
        format_string = self.get_format_string(postfix,ext,absPath)
        for t in self.get_range():
            yield t, (format_string % t)

    def mean_Nb(self):
        """Calculate the time averaged number of particles in a frame"""
        if not hasattr(self,'__mean_Nb'):
            nb = 0L
            for t,name in self.enum():
                with open(name,'r') as f:
                    nb += int(re.split("\t",f.readline())[1])
            self.__mean_Nb = float(nb) / self.size
        return self.__mean_Nb

    def mean_V(self):
        """Calculate the time averaged volume"""
        if not hasattr(self,'__mean_V'):
            V=0
            for t,name in self.enum():
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
        name = os.path.join(self.path, self.head + '.rdf')
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
        nbs = np.empty((self.size))
        Vs = np.empty((self.size))
        for t,fname in self.enum():
            coords = np.loadtxt(fname,delimiter='\t', skiprows=2)
            nbs[t-self.offset] = len(coords)
            Vs[t-self.offset] = np.ptp(coords, axis=0).prod()
        if averaged:
            return (nbs/Vs).mean()
        else:
            return (nbs, Vs)

    def get_zPortion_Nbd(self, lowerMargin=0, upperMargin=0, averaged=True):
        """Get the number density of a z-slab"""
        nbs = np.empty((self.size))
        Vs = np.empty((self.size))
        for t,fname in self.enum():
            coords = np.loadtxt(fname,delimiter='\t', skiprows=2)
            m = np.amin(coords[:,-1])+lowerMargin
            M = np.amax(coords[:,-1])-upperMargin
            coords = coords[np.bitwise_and(m>coords[:,-1], coords[:,-1]<M)]
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
        return [(t,self.offset+self.size-average, average) \
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
            for t, name in self.enum(ext='g6'):
                if not os.path.exists(name):
                    force=True
                    break
        if force:
            for t, name in self.enum():
                subprocess.check_call(map(str,
                    ['g6', name, self.radius, Nbins, nbDiameters]))
        tot = np.zeros((Nbins,3))
        for t, name in self.enum(ext='g6'):
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
        #keep only positive maxima further than the largest maxima
        maxima = [m for m in maxima[g6.argmax():] if g6[m]>0]
        envelope = np.column_stack((r[maxima],g6[maxima]))
        np.savetxt(
            os.path.join(self.path,self.head + '_total_env.g6'),
            envelope,
            fmt='%f',
            delimiter='\t'
            )
        return envelope

    def lost_ngb_profile(self, t, Nbins=50, vf=False):
        """output the lost neighbour profile. Default unit is pixel^-3, or volume fraction"""
        lngb = np.loadtxt(self.get_format_string(ext='lngb')%t)
        pos = np.loadtxt(self.get_format_string()%t, skiprows=2)
        H, xedges, yedges = np.histogram2d(pos[:,-1], lngb, bins=[Nbins,2])
        H /= pos.ptp(axis=0).prod()/Nbins
        if(vf):
            H *= 4 * math.pi * self.radius**3 /3
        np.savetxt(
            self.get_format_string('_lngb','hist')%t,
            np.column_stack((xedges[:-1], H))
            )

class Txp:
    """Implementig time algorithms in python"""
    
    def __init__(self, xp=None, start=None, size=None, copy=None):
        if copy is not None:
            self.xp = copy.xp
            self.trajs = np.copy(copy.trajs)
            self.positions = np.copy(copy.positions)
            return
        if xp is not None:
            self.xp = xp
            if not start or start < self.xp.offset:
                start = self.xp.offset
            if not size or start+size > self.xp.offset+self.xp.size:
                size = self.xp.size + self.xp.offset - start
            self.trajs = self.read_trajs(start, size)
            self.positions = self.read_pos(start, size)
            self.remove_drift()
            
    def __getitem__(self, indices):
        """get a copy of a selection of trajectories"""
        c = Txp(copy=self)
        c.trajs = c.trajs[indices]
        c.positions = c.positions[:, indices]
        return c

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
        for t, fname in self.xp.enum():
            if t<start or t>= start+size:
                continue
            raw_pos = np.loadtxt(fname,delimiter='\t',skiprows=2)
            pos[t-start] = raw_pos[self.trajs[:,t-start]]
        return pos

    def load_bonds(self, t):
        oldbonds = np.loadtxt(self.xp.get_format_string(ext='bonds')%t, dtype=long)
        pos2traj = dict([(p,tr) for tr,p in enumerate(self.trajs[:,t])])
        newbonds = np.asarray([
            np.sort([pos2traj[a], pos2traj[b]])
            for a,b in oldbonds
            if pos2traj.has_key(a) and pos2traj.has_key(b)
            ])
        indices = np.lexsort((newbonds[:,1], newbonds[:,0]))
        return newbonds[indices]

    def remove_drift(self):
        """Remove the average drift between time steps"""
        drift = np.cumsum(np.diff(self.positions,axis=0).mean(axis=1),axis=0)
        sw = np.swapaxes(self.positions[1:],0,1)
        sw -= drift

    def z_slab(self, bottom, top, allTimes=True):
        """remove all trajectories that are not in the slab
            defined by [bottom, top]"""
        if allTimes:
            selection = np.unique1d(np.where(
                np.bitwise_and(
                    self.positions[:,:,-1]>bottom,
                    self.positions[:,:,-1]<top
                    ))[1])
        else:
            selection = np.unique1d(np.where(
                np.bitwise_and(
                    self.positions[:,:,-1].max(axis=0)>bottom,
                    self.positions[:,:,-1].min(axis=0)<top
                    )))
        self.trajs = self.trajs[selection]
        self.positions = self.positions[:,selection]

    def exclude_null(self, postfix='_space',ext='cloud', col=1):
        """Remove trajectories having at least a null value in the field given by postfix, ext and col"""
        field = np.zeros((self.trajs.shape[1], self.trajs.shape[0]))
        for t, fname in self.xp.enum(postfix=postfix, ext=ext):
            field[t] = np.loadtxt(fname, usecols=[col])[self.trajs[:,t]]
        selection = np.where(field.min(axis=0)!=0.0)[0]
        self.trajs = self.trajs[selection]
        self.positions = self.positions[:,selection]
        

    def get_vtk(self, t):
        """Get a view of a time step as a vtk file. remove_drift recommended before use."""
	v=vtk.Polydata()
	v.points=self.positions[t]
	v.bonds=self.load_bonds(t)
	cgcloud = np.loadtxt(
		self.xp.get_format_string(postfix='_space', ext='cloud')%t
		)[self.trajs[:,t]]
	v.scalars.append(('Q6',cgcloud[:,1]))
	v.scalars.append(('W4',cgcloud[:,4]))
	cloud = np.loadtxt(
		self.xp.get_format_string(ext='cloud')%t
		)[self.trajs[:,t]]
	v.scalars.append(('w6',cloud[:,5]))
	v.scalars.append(('w10',cloud[:,7]))
	lngb = np.loadtxt(
		self.xp.get_format_string(postfix='_post', ext='lngb')%t
		)[self.trajs[:,t]]
	v.scalars.append(('lngb',lngb))
	phi = np.loadtxt(
		self.xp.get_format_string(postfix='_space', ext='phi')%t,
		skiprows=2
		)[self.trajs[:,t]]
	v.scalars.append(('phi',phi))
	vel = np.loadtxt(
		self.xp.get_format_string(ext='vel')%t,
		skiprows=1
		)[self.trajs[:,t]]
	v.vectors.append(('vel', vel-vel.mean(axis=0)))
	return v


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
    def nonGaussian(self,start,stop,av):
        """
        Non gaussian parameter
        If av is 0 (Default), the calculation will act greedily,
        averaging over all the avilabe intervals of a given length.
            Example : start=1 stop=4 av=0
                alpha[dt=0] = 1
                alpha[dt=1] = ( alpha([1,2]) + alpha([2,3]) + alpha([3,4]))/3
                alpha[dt=2] = ( alpha([1,3]) + alpha([2,4]) )/2
                alpha[dt=3] = msd([1,4])
        If av>0, the average will be done over av time intervals starting
        from start, start+1,...,start+av-1
            Example : start=1 stop=4 av=2
                alpha[dt=0] = 1
                alpha[dt=1] = ( alpha([1,2]) + alpha([2,3]) )/2
                alpha[dt=2] = ( alpha([1,3]) + alpha([2,4]) )/2
                alpha[dt=3] = ( alpha([1,4]) + alpha([2,5]) )/2
    """
        A = self.positions[start:stop+av]
        msd = np.zeros((stop-start))
        mqd = np.zeros((stop-start))
        if av==0:
            for t0, a in enumerate(A):
                for dt, b in enumerate(A[t0+1:]):
                    #average is done over all trajectories and the 3 dimensions
                    msd[dt+1] += ((b-a)**2).sum()
                    mqd[dt+1] += ((b-a)**4).sum()
            msd /= A.shape[1] * A.shape[2] * (self.xp.radius*2)**2
            mqd /= A.shape[1] * A.shape[2] * (self.xp.radius*2)**4
            for dt, n in enumerate(range(stop-start,0,-1)):
                msd[dt+1] /= n
                mqd[dt+1] /=n
        else:
            for t0, a in enumerate(A[:av]):
                for dt, b in enumerate(A[t0+1:-av+t0]):
                    msd[dt+1] += ((b-a)**2).sum()
                    mqd[dt+1] += ((b-a)**4).sum()
            msd /= av * A.shape[1] * A.shape[2]  * (self.xp.radius*2)**2
            mqd /= av * A.shape[1] * A.shape[2]  * (self.xp.radius*2)**4
        mqd[1:] /= (3 * msd[1:]**2)
        return mqd-1

    def export_dynamics(self,start,stop,av):
        self.export_msd(start,stop,av)
        self.export_self_isf(start,stop,av)

    def time_correlation(self, postfix='',ext='dat', col=0, av=10):
        """read the particle-wise scalar from a time serie of files and compute the time correlation"""
        data = np.zeros((self.trajs.shape[1], self.trajs.shape[0]))
        for t, fname in self.xp.enum(postfix=postfix, ext=ext):
            data[t] = np.loadtxt(fname, usecols=[col])[self.trajs[:,t]]
        data -= data.mean()
        c=np.zeros((data.shape[0]-av+1))
        if av==0:
            for t0, a in enumerate(data):
                for dt, b in enumerate(data[t0:]):
                    #average is done over all trajectories
                    c[dt] += (b*a).mean()
            for dt, n in enumerate(range(c.shape[1],0,-1)):
                c[dt] /= n
            return c/c[0]
        else:
            for t0, a in enumerate(data[:av]):
                for dt, b in enumerate(data[t0:1-av+t0]):
                    c[dt] += (b*a).mean()
            c /= av
            return c/c[0]
        
        
    


def br(sigma, T=310, eta=2.22e-3):
        """Brownian time for a particle of diameter sigma (in meters)"""
        return 3 * const.pi * eta * (sigma**3) / (4 * const.k * T)

def histz(f):
    """export the density histogram of a .dat file into a .z file"""
    hist, bins = np.histogram(np.loadtxt(f, delimiter='\t', skiprows=2, usecols=[2]))
    np.savetxt(f[:-3]+'z', hist/(bins[1]-bins[0]), fmt='%f', delimiter='\t')

def dilate(field, bonds):
	ngb = [[n] for n in range(field.shape[0])]
	for b in bonds:
		ngb[b[0]].append(b[1])
		ngb[b[1]].append(b[0])
	dil = np.zeros_like(field)
	for p,n in enumerate(ngb):
		a = field[n]
		dil[p] = a[np.absolute(a).argmax(axis=0), range(a.shape[1])]
	return dil

def average(field, bonds):
	ngb = [[] for n in range(field.shape[0])]
	for b in bonds:
		ngb[b[0]].append(b[1])
		ngb[b[1]].append(b[0])
	av = np.zeros_like(field)
	for p,n in enumerate(ngb):
		av[p] = (field[p]+field[n].mean(axis=0))/2
	return av
