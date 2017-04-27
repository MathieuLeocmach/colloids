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
from __future__ import print_function #for python 2.x
import scipy as sp
import numpy as np
import scipy.constants as const
from scipy.ndimage.filters import gaussian_filter1d
from scipy.ndimage.morphology import grey_dilation
from scipy import optimize, sparse
from scipy.optimize import leastsq
from scipy.spatial import cKDTree as KDTree
import os, os.path, subprocess
import re, string, math
from math import exp
from colloids import vtk, statistics
from colloids.progressbar import ProgressBar
from colloids.particles import bonds2ngbs_list
import networkx as nx
import numexpr

def noNaN(x):
    if x == 'NaN':
        return np.nan
    else:
        return float(x)

def trajectory_stream(fname):
    """Generator yielding (starting time, list of position indices) for each trajectory"""
    with open(fname,'r') as f:
        for l, line in enumerate(f):
            if l==3:
                break
        for l, line in enumerate(f):
            if l%2 == 0:
                start = int(line[:-1])
            else:
                positions = np.fromstring(line[:-1], int, sep='\t')
                yield start, positions

def load_trajectories(fname):
    """returns (starting times, list of position indices)"""
    starts = []
    positions = []
    for start, pos in trajectory_stream(fname):
        starts.append(start)
        positions.append(pos)
    return starts, positions

class Experiment(object):
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
            files = [name for name in os.listdir(self.path) if mask.match(name) is not None]
            self.offset = 0
            self.size = len(files)
            for m in map(mask.match,files):
                if len(m.group(2)) > self.digits:
                    self.digits = len(m.group(2))
                    self.token = m.group(1)
            #if last frame is empty, reduce size
            lastname = self.get_format_string()%(self.size-1)
            if int(open(lastname,'r').readline().split()[1])==0:
                self.size -= 1
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
        trajfile = os.path.join(self.path,self.trajfile)
        with open(trajfile,'r') as f:
            self.radius,self.dt = [
                float(s) for s in re.split("\t",f.readline()[:-1])
                ]
            pattern, = re.split("\t",f.readline()[:-1])
            self.token, = re.split("\t",f.readline()[:-1])
            self.offset,self.size = [
                int(s) for s in re.split("\t",f.readline()[:-1])
                ]
            self.read_pattern(pattern)
    
    def fill_traj_data(self):
        """read the trajectory data in itself"""
        trajfile = os.path.join(self.path,self.trajfile)
        self._starts, self._trajs = load_trajectories(trajfile)
    
    @property
    def starts(self):
        if not hasattr(self,'_starts'):
            self.fill_traj_data()
        return self._starts
        
    @property
    def trajs(self):
        if not hasattr(self,'_trajs'):
            self.fill_traj_data()
        return self._trajs
        
    @property
    def nb_trajs(self):
        if not hasattr(self,'_nb_trajs'):
            if hasattr(self,'_trajs'):
                self._nb_trajs = len(self._trajs)
            else:
                with open(os.path.join(self.path,self.trajfile),'r') as f:
                    for l, line in enumerate(f):
                        if l==3:
                                break
                        tr=0
                    for line in f:
                        tr += 1
                self._nb_trajs =  tr/2
        return self._nb_trajs
        
    def p2tr(self, t):
        fname = self.get_format_string(ext='p2tr')%t
        if not os.path.exists(fname):
            pos2tr = [np.zeros(n, int) for n in self.get_nb()]
            for tr, (start, pos) in enumerate(trajectory_stream(
                os.path.join(self.path,self.trajfile)
                )):
                for dt, p in enumerate(pos):
                    pos2tr[start+dt][p] = tr
            for (t, name), trs in zip(self.enum(ext='p2tr'), pos2tr):
                np.savetxt(name, trs, fmt='%d')
        return np.loadtxt(fname, dtype=int)
        
    def load_all(self, removedrift=False, showprogress=True):
        """Returns all positions from all trajectories"""
        trajpos = [np.zeros([len(tr), 3]) for tr in self.trajs]
        if showprogress:
            pro = ProgressBar(self.size)
        for t, name in self.enum():
            pos = np.loadtxt(name, skiprows=2)
            p2tr = self.p2tr(t)
            for p, tr in enumerate(p2tr):
                trajpos[tr][t - self.starts[tr]] = pos[p]
            if showprogress:
                pro.animate(t)
        if removedrift:
            #compute drift using all trajectories
            drift = np.zeros((self.size,3))
            nb = np.zeros((self.size), int)
            for tr, start in zip(trajpos, self.starts):
                drift[start+1:start + len(tr)] += tr[1:] - tr[:-1]
                nb[start+1:start + len(tr)] += 1
            drift /= np.maximum(1, nb)[:,None]
            drift = np.cumsum(drift, axis=0)
            #remove drift from each trajectory
            for tr, start in enumerate(self.starts):
                trajpos[tr] -= drift[start:start + len(trajpos[tr])]
        return trajpos

    def read_pattern(self,pattern):
        m = re.match('(.*)'+self.token+'([0-9]*)',os.path.splitext(pattern)[0])
        self.head = m.group(1)
        self.digits = len(m.group(2))

    def cut(self):
        """Remove centers closer than 3 pixels in each frame."""
        subprocess.check_call(
                list(map(str,['cutter',self.get_format_string()%0, 
                   self.token, 5, self.offset, self.size, 0.3]))
                )

    def link(self):
        """Link trajectories. Export a file.traj"""
        actual = os.getcwd()
        os.chdir(self.path)
        subprocess.check_call(list(map(str,
                 ['linker', 
                 self.get_format_string(absPath=False)%0, self.token,
                 self.radius, self.dt,
                 self.offset,self.size]))
            )
        os.chdir(actual)

    def linkboo(self):
        """calculate total g(r), radius, BOO for each time step and link trajectories."""
        actual = os.getcwd()
        os.chdir(self.path)
        subprocess.check_call(list(map(str,
                 ['linkboo',
                  self.get_format_string(absPath=False)%0,
                  self.token,
                  self.dt,self.size,
                  self.offset]))
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
            
    def get_nb(self):
        """Number of particles in each frame"""
        if not hasattr(self,'__Nb'):
            nbname = os.path.join(self.path, self.head + '.nb')
            if not os.path.isfile(nbname):
                self.__Nb = np.array([
                    int(open(name,'r').readline().split()[1])
                    for t,name in self.enum()
                    ], int)
                np.savetxt(nbname, self.__Nb, fmt='%d')
            else:
                self.__Nb = np.loadtxt(nbname).astype(int)
        return self.__Nb
    
    def get_nb_bonds(self):
        """Number of bonds in each frame"""
        if not hasattr(self,'__Nb_bonds'):
            self.__Nb_bonds = np.array([
                len([line for line in open(name,'r')])
                for t,name in self.enum(ext='bonds')
                ], int)
        return self.__Nb_bonds

    def mean_Nb(self):
        """Calculate the time averaged number of particles in a frame"""
        return np.mean(self.get_nb())

    def mean_V(self):
        """Calculate the time averaged volume"""
        if not hasattr(self,'__mean_V'):
            V=0
            for t,name in self.enum():
                V += np.ptp(
                    np.loadtxt(name, skiprows=2),
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
            subprocess.check_call(list(map(str,
                  ['totalRdf',self.get_format_string()%0, self.token,
                   200, 15]))
                )
        r,g = np.loadtxt(name, unpack=True)
        return r[np.argmax(g)]
        gm = g.argmax()
        return np.average(r[gm-1:gm+2],weights=g[gm-1:gm+2])

    def get_Nb_density(self, averaged=True):
        nbs = np.empty((self.size))
        Vs = np.empty((self.size))
        for t,fname in self.enum():
            coords = np.loadtxt(fname, skiprows=2)
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
            coords = np.loadtxt(fname, skiprows=2)
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
            t,isf = np.loadtxt(name,usecols=[0,4],unpack=True)
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
        
    def msd(self, t0=0, dts=None, showProgress=False):
        """Compute mean square displacement from one initial time, given an array of interval times. By default dts are logarythmically spaced."""
        assert t0 >= 0
        assert t0 < self.size-1
        #generate default dts
        if dts is None:
            dts = np.unique(np.logspace(0, np.log10(self.size-t0-1),100).astype(int))
        assert t0 + max(dts) < self.size
        msd = np.zeros(len(dts))
        if showProgress:
            pro = ProgressBar(len(msd))
        #load initial time
        p0 = np.loadtxt(self.get_format_string()%(t0), skiprows=2)
        p2tr0 = self.p2tr(t0)
        for idt, dt in enumerate(dts):
            p1 = np.loadtxt(self.get_format_string()%(t0+dt), skiprows=2)
            p2tr1 = self.p2tr(t0+dt)
            #select trajectories that span both time steps
            span = np.zeros(max([p2tr0.max(), p2tr1.max()])+1, bool)
            span[np.intersect1d(p2tr0, p2tr1)] = True
            #displacement between t0 and t1
            dx = np.zeros((len(span),3))
            dx[p2tr1[span[p2tr1]]] = p1[span[p2tr1]]
            dx[p2tr0[span[p2tr0]]] -= p0[span[p2tr0]]
            #remove drift
            dx[span] -= dx[span].mean()
            #msd for this dt
            msd[idt] = np.mean(dx[span]**2)
            if showProgress:
                pro.animate(idt)
        return dts, msd

        
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
                subprocess.check_call(list(map(str,
                    ['g6', name, self.radius, Nbins, nbDiameters])))
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
        lngb = np.loadtxt(self.get_format_string('_post', ext='lngb')%t)
        pos = np.loadtxt(self.get_format_string()%t, skiprows=2)
        H, xedges, yedges = np.histogram2d(pos[:,-1], lngb, bins=[Nbins,2])
        H /= pos.ptp(axis=0).prod()/Nbins
        if(vf):
            H *= 4 * math.pi * self.radius**3 /3
        np.savetxt(
            self.get_format_string('_lngb','hist')%t,
            np.column_stack((xedges[:-1], H))
            )
        
    def traj_crea_dest(self):
        """The creation and destruction time for each trajectory"""
        crea = []
        dest = []
        with open(os.path.join(self.path,self.trajfile),'r') as f:
            for l, line in enumerate(f):
                if l==3:
                    break
            for line in f:
                t0 = int(line[:-1])
                crea.append(t0)
                dest.append(t0+len(string.split(f.next()[:-1],'\t')))
        return np.column_stack((crea, dest))
                            
    def vtk_tracking(self):
        #read all trajectories
        trajs = []
        with open(os.path.join(self.path,self.trajfile),'r') as f:
            for l, line in enumerate(f):
                if l==3:
                    break
            for line in f:
                t0 = int(line[:-1])
                pos = string.split(f.next()[:-1],'\t')
                trajs.append((t0, list(map(int, pos))))
        #mapping (t, pos)->traj
        framesizes = np.zeros(self.size, dtype=int)
        for t0, pos in trajs:
            framesizes[t0:t0+len(pos)]+=1
        pos2traj = [np.empty(s, dtype=int) for s in framesizes]
        for tr, (t0, pos) in enumerate(trajs):
            for t, p in enumerate(pos):
                pos2traj[t0+t][p]=tr
        #creating and saving vtk
        for t, fname in self.enum():
            v = vtk.Polydata()
            v.bonds = []
            v.points = np.loadtxt(fname, skiprows=2)
            v.scalars.append((
                'traj_length',
                np.asarray([len(trajs[tr][1]) for tr in pos2traj[t]])
                ))
            v.save(self.get_format_string('_check', ext='vtk')%t)
            
    def broken_bonds(self, t0, t1):
        """bonds (between trajectories) existing at t0 but no more at t1 with both members still existing at t1 = broken bonds - lost trajectories"""
        #load trajectory data
        p2tr0 = self.p2tr(t0)
        p2tr1 = self.p2tr(t1)
        #which trajectories exist at both times step
        existin0 = np.zeros(max(p2tr0.max(), p2tr1.max())+1, bool)
        existin1 = np.zeros_like(existin0)
        existin0[p2tr0] = True
        existin1[p2tr1] = True
        existinboth = existin0 & existin1
        #load bonds, translate to trajectory indices
        gb0, gb1 = (
            p2tr[np.loadtxt(self.get_format_string(ext='bonds')%t, dtype=int)]
            for t, p2tr in zip([t0,t1], [p2tr0, p2tr1])
            )
        #restrict to these trajectories, sort for unicity
        gb0, gb1 = (
            np.sort(gb[existinboth[gb].min(-1)], 1)
            for gb in (gb0, gb1)
            )
        #bonds existing at t but not t+dt
        return set([(a,b) for a,b in gb0]) - set([(a,b) for a,b in gb1])
    
    def bond_life(self):
        """When do bonds first form and are last seen"""
        bstart = sparse.dok_matrix(tuple([self.nb_trajs]*2), int)
        blast = sparse.dok_matrix(tuple([self.nb_trajs]*2), int)
        for t,name in self.enum(ext='bonds'):
            bonds = np.sort(self.p2tr(t)[np.loadtxt(name, int)], 1)
            for a,b in bonds:
                blast[a,b] = t
                if not bstart.has_key((a,b)):
                    bstart[a,b] = t
        return bstart, blast
        
    def bond_last(self, persistence=10):
        blast = dict()
        for t,name in self.enum(ext='bonds'):
            bonds = np.sort(self.p2tr(t)[np.loadtxt(name, int)], 1)
            for a,b in bonds:
                blast[(a,b)] = t
            out = dict()
            for k, tf in blast.iteritems():
                if tf+persistence<t:
                    out[k]=tf
            for k in out.iterkeys():
                del blast[k]
            yield out
        yield blast
        
    def neverreform(self, L=3, mintrajlength=None, showprogress=True):
        """All the broken bonds (in term of trajectories) with a on graph post breaking distance larger than L, associated with the last time they break"""
        try:
            ret = []
            for t, name in self.enum('_L%d'%L, ext='nev'):
                if os.stat(name).st_size:
                    ret.append(np.loadtxt(name, dtype=int))
                else:
                    ret.append(np.array((0,2), dtype=int))
        except IOError:
            pass
        if showprogress:
            pro = ProgressBar(self.size)
        result = dict()
        for t, name in self.enum(ext='bonds'):
            try:
                bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
            except UserWarning:
                bonds1 = np.zeros([0,2], int)
            p2tr1 = self.p2tr(t)
            if t>0:
                for path in broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1):
                    if len(path)<=L:
                        continue
                    a, b = sorted([path[0], path[-1]])
                    if mintrajlength is not None and min(len(self.trajs[u]) for u in [a,b]) <= mintrajlength:
                        continue
                    result[(a,b)] = t #create or replace if was already broken before
                #remove all bonds that exist at t1
                for a,b in np.sort(p2tr1[bonds1], 1):
                    result.pop((a,b),-1) 
            bonds0 = bonds1
            p2tr0 = p2tr1
            if showprogress:
                pro.animate(t)
        #translate as a list of arrays of bonds, an array per time
        tdl = [
            np.zeros([nbb, 2], int) 
            for nbb in np.histogram(result.values(), bins=np.arange(0, self.size+1))[0]
            ]
        ns = np.zeros(self.size, int)
        for bond, t in result.iteritems():
            tdl[t][ns[t]] = bond
            ns[t] += 1
        #save a cache
        for t, name in self.enum('_L%d'%L, ext='nev'):
            np.savetxt(name, tdl[t], fmt='%d')
        return tdl
        
    def never_closer(self, L=3, mintrajlength=None, showprogress=True):
        """All the broken bonds (in term of trajectories) with a on graph distance that stays larger than L, associated with the last time they break"""
        try:
            ret = []
            for t, name in self.enum('_L%d'%L, ext='nevc'):
                if os.stat(name).st_size:
                    ret.append(np.loadtxt(name, dtype=int))
                else:
                    ret.append(np.array((0,2), dtype=int))
        except IOError:
            pass
        if showprogress:
            pro = ProgressBar(self.size)
        result = dict()
        for t, name in self.enum(ext='bonds'):
            try:
                bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
            except UserWarning:
                bonds1 = np.zeros([0,2], int)
            p2tr1 = self.p2tr(t)
            if t>0:
                for path in broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1):
                    if len(path)<=L:
                        continue
                    a, b = sorted([path[0], path[-1]])
                    if mintrajlength is not None and min(len(self.trajs[u]) for u in [a,b]) <= mintrajlength:
                        continue
                    result[(a,b)] = t #create or replace if was already broken before
                #graph of the bonds between trajectories at t1
                g = nx.Graph()
                g.add_nodes_from(p2tr1)
                g.add_edges_from(p2tr1[bonds1])
                #remove all bonds that are closer than L at t1
                existing = [[a,b] for a,b in result.iterkeys()]
                for path in shortest_path(g, existing):
                    if len(path) <= L:
                        a,b = sorted([path[0], path[-1]])
                        result.pop((a,b),-1)
            bonds0 = bonds1
            p2tr0 = p2tr1
            if showprogress:
                pro.animate(t)
        #translate as a list of arrays of bonds, an array per time
        tdl = [
            np.zeros([nbb, 2], int) 
            for nbb in np.histogram(result.values(), bins=np.arange(0, self.size+1))[0]
            ]
        ns = np.zeros(self.size, int)
        for bond, t in result.iteritems():
            tdl[t][ns[t]] = bond
            ns[t] += 1
        #save a cache
        for t, name in self.enum('_L%d'%L, ext='nevc'):
            np.savetxt(name, tdl[t], fmt='%d')
        return tdl
        
    def neighbouring_bonds(self, t, bonds_of_interest):
        """Returns the set of bonds that are neighbouring the bonds of interest at time t. Input and output are in terms of trajectories."""
        p2tr = self.p2tr(t)
        #load all the bonds at time t
        bonds = np.loadtxt(self.get_format_string(ext='bonds')%t, int)
        #construct neighbours for each particle
        ngb = bonds2ngbs_list(bonds, p2tr.shape[0])
        result = []
        for tra, trb in bonds_of_interest:
            for tr in [tra, trb]:
                a = self.trajs[tr][t - self.starts[tr]]
                for ngb_tr in p2tr[ngb[a]]: #iterate on neighbours
                    #exclude the bonds of interest
                    if ngb_tr in [tra, trb]: continue
                    result.append(sorted([tr, ngb_tr]))
        #convert into a nice array of bonds
        return np.array(result).reshape([len(result), 2]).astype(int)
        

class Txp:
    """Contains spanning trajectories as a (time, particle, dimention) array. Implements time-dependent algorithms."""
    
    def __init__(self, xp=None, start=None, size=None, copy=None):
        if copy is not None:
            self.xp = copy.xp
            self.trajs = np.copy(copy.trajs)
            self.positions = np.copy(copy.positions)
            return
        if xp is not None:
            self.xp = xp
            if start is None or start < self.xp.offset:
                self.start = self.xp.offset
            else:
                self.start=start
            if size is None or self.start+size > self.xp.offset+self.xp.size:
                self.size = self.xp.size + self.xp.offset - self.start
            else:
                self.size= size
            self.trajs = self.read_trajs()
            assert len(self.trajs)>0
            self.positions = self.read_pos()
            self.remove_drift()
            
    def __getitem__(self, indices):
        """get a copy of a selection of trajectories"""
        c = Txp(copy=self)
        c.trajs = c.trajs[indices]
        c.positions = c.positions[:, indices]
        return c

    def read_trajs(self):
        """
            Reads the linked trajectories from the .traj file
            Retrieves only trajectories spanning [start,start+size[
            """
        starts, positions = load_trajectories(
            os.path.join(self.xp.path,self.xp.trajfile)
            )
        lengths = np.array([len(p) for p in positions])
        begin = self.start - np.array(starts)
        size = lengths - begin
        return np.asarray([
            pos[b:b+self.size]
            for pos, b, s in zip(positions, begin, size)
            if b>=0 and s>=self.size
            ])

    def read_pos(self):
        """Reads the usefull positions from the .dat files"""
        pos = np.empty((self.trajs.shape[1],self.trajs.shape[0],3))
        for t, fname in self.xp.enum():
            if t<self.start or t>= self.start+self.size:
                continue
            raw_pos = np.loadtxt(fname, skiprows=2)
            pos[t-self.start] = raw_pos[self.trajs[:,t-self.start]]
        return pos

    def load_bonds(self, t):
        oldbonds = np.loadtxt(self.xp.get_format_string(ext='bonds')%t, dtype=int)
        pos2traj = -np.ones((oldbonds.max()+1), dtype=int)
        pos2traj[self.trajs[:,t]] = range(self.trajs.shape[0])
        newbonds = pos2traj[oldbonds]
        newbonds = newbonds[np.where(newbonds.min(axis=1)>-1)]
        newbonds.sort(axis=1)
        indices = np.lexsort((newbonds[:,1], newbonds[:,0]))
        return newbonds[indices]

    def load_q6m(self):
        q6m = np.zeros((self.positions.shape[0], self.positions.shape[1], 7), dtype=complex)
        for t, fname in self.xp.enum(ext='qlm'):
            a = np.loadtxt(fname, usecols=range(18,32))[self.trajs[:, t]]
            q6m[t] = a[:,::2] + 1j * a[:,1::2]
        return q6m

    def load_Q6m(self):
        q6m = np.zeros((self.positions.shape[0], self.positions.shape[1], 7), dtype=complex)
        Q6m = np.zeros_like(q6m)
        for t, fname in self.xp.enum(ext='qlm'):
            A = np.loadtxt(fname, usecols=range(18,32))
            B = np.copy(A)
            nb = np.ones((len(B),1), dtype=int)
            bonds = np.loadtxt(self.xp.get_format_string(ext='bonds')%t, dtype=int)
            B[bonds[:,0]] += A[bonds[:,1]]
            B[bonds[:,1]] += A[bonds[:,0]]
            nb[bonds.ravel()] +=1
            B /= nb
            q6m[t] = A[self.trajs[:,t]][:, ::2] + 1j * A[self.trajs[:,t]][:, 1::2]
            Q6m[t] = B[self.trajs[:,t]][:, ::2] + 1j * B[self.trajs[:,t]][:, 1::2]
        return q6m, Q6m

    def remove_drift(self, on=True):
        """Remove the average drift between time steps"""
        drift = np.cumsum(np.diff(self.positions,axis=0).mean(axis=1),axis=0)
        sw = np.swapaxes(self.positions[1:],0,1)
        sw -= drift
#        if on:
 #           if not hasattr(self, 'drift'):
  #              self.drift = np.cumsum(np.diff(self.positions,axis=0).mean(axis=1),axis=0)
   #         self.positions[1:] -= self.drift[:,None,:]
    #    else:
     #       if hasattr(self, 'drift'):
      #          self.positions[1:] += self.drift[:,None,:]
       #     else:
        #        raise UserWarning("Can't add back a drift that was not previously removed")

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
            if t<self.trajs.shape[1]:
                field[t] = np.loadtxt(fname, usecols=[col])[self.trajs[:,t]]
        selection = np.where(field.min(axis=0)!=0.0)[0]
        self.trajs = self.trajs[selection]
        self.positions = self.positions[:,selection]

    def get_Nb_density(self, t):
        """Get the number density (px**-3) in the volume actually used in the time experiment"""
        pos = np.loadtxt(self.xp.get_format_string()%t, skiprows=2)
        inself = pos[self.trajs[:,t]]
        inside = pos[np.bitwise_and(
            (pos>=inself.min(axis=0)).min(axis=-1),
            (pos<=inself.max(axis=0)).min(axis=-1)
            )]
        return len(inside)/inside.ptp(axis=0).prod()
        

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
            msd /= A.shape[1] * (self.xp.radius*2)**2
            for dt in range(len(A)):
                msd[dt] /= len(A)-dt
            return msd
        else:
            for t0, a in enumerate(A[:av]):
                for dt, b in enumerate(A[t0+1:-av+t0]):
                    msd[dt+1] += ((b-a)**2).sum()
            msd /= av * A.shape[1] * (self.xp.radius*2)**2
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
        A = numexpr.evaluate('exp(r * q)',{
            'r':self.positions[start:stop+av-1],
            'q': (1j * np.pi / self.xp.radius)}
            )
        return statistics.time_correlation(A, av).mean(axis=-1).real
        
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
                    diff = ((b-a)**2).sum(axis=-1)
                    msd[dt+1] += diff.sum()
                    mqd[dt+1] += (diff**2).sum()
            for dt in range(len(mqd)):
                mqd[dt] *= (len(mqd)-dt) * A.shape[1]
        else:
            for t0, a in enumerate(A[:av]):
                for dt, b in enumerate(A[t0+1:-av+t0]):
                    diff = ((b-a)**2).sum(axis=-1)
                    msd[dt+1] += diff.sum()
                    mqd[dt+1] += (diff**2).sum()
            mqd *= av * A.shape[1]
        mqd[1:] /= (5 * msd[1:]**2)/3.
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

    def cor_q6m(self, q6m, start, stop):
        A = np.exp(
            self.positions * (1j * np.pi / self.xp.radius)
            )
        B = np.atleast_3d(A.mean(axis=-1)) * q6m
        C = np.real(B[start].conj() * B[start:stop])
        return (C.sum(axis=2)-C[:,:,0]).mean(axis=1)
        
    def get_Brownian_corr(self, t0, dts, rmin, rmax, nbins, inside=None, verbose=False):
        """Calculates 'Brownian correlation' components.
        See JC Crocker, MT Valentine, ER Weeks, T Gisler, PD Kaplan, AG Yodh, and DA Weitz, Phys. Rev. Lett. 85, 888 (2000).

        Input
            - t0: initial time
            - dts: 1D-array of time lags
            - rmin: minimim distance between particles
            - rmax: maximum distance between particles
            - nbins: number of bins
        
        Output:
            - r dependent sums of r-r, theta-theta and phi-phi components and sums of squares
            - number of pairs used for the sum
            - r bins
        
        Code Inspired by John C. Crocker IDL code http://www.physics.emory.edu/faculty/weeks//idl/rheo.html
        """
        # generate the 'r' partition in log scale
        lrmin, lrmax = np.log([rmin, rmax])
        lrbinsize = (lrmax - lrmin) / nbins
        lrbins = np.arange(nbins+1) * lrbinsize
        #accumulate both sums and squared sums
        bres = np.zeros((2*self.positions.shape[-1], len(dts), nbins))
        nb = np.zeros(nbins, int)
        #vector displacements
        pos0 = self.positions[t0]
        displs = self.positions[t0+np.array(dts)] - pos0
        #3D only
        assert self.positions.shape[-1] == 3, "Implemented only for 3D data"
        
        if inside is None:
            inside = np.arange(len(pos0))
        
        #Spatial indexing
        tree = KDTree(pos0, 12)
        if verbose:
            pro = ProgressBar(len(inside)-1)
        for it,i in enumerate(inside):
            if verbose:
                pro.animate(it)
            js = tree.query_ball_point(pos0[i], rmax)
            rs = pos0[js] - pos0[i]
            distsq = np.sum(rs**2, -1)
            ok = (distsq >= rmin**2) #& (distsq < rmax**2)
            js = np.array(js)[ok]#np.where(ok)[0] #+ i + 1
            #distances
            dists = np.sqrt(distsq[ok])
            #vector separations
            rs = rs[ok]
            #unit vectors along r
            ns = rs / dists[:,None]
            #unit vectors in xy plane: projection
            ps = rs[:, :-1] / np.sqrt(np.sum(rs[:, :-1]**2, -1))[:,None]
            ps = np.column_stack((ps[:,1], -ps[:,0])) #orthogonal
            #third unit vectors
            qs = np.cross(ns, np.column_stack([ps[:,0], ps[:,1], np.zeros(len(ps))]))
            #randomize the directions of the unit vectors to help with numerical errors
            #actually useless, commented
            #ns *= (2*np.random.randint(2, size=len(ns)) -1)[:,None]
            #ps *= (2*np.random.randint(2, size=len(ns)) -1)[:,None]
            #qs *= (2*np.random.randint(2, size=len(ns)) -1)[:,None]
            #calculate the longitudinal part
            ddn = (displs[...,None,i,:] * ns).sum(-1) * (displs[...,js,:] * ns).sum(-1)
            # calculate the theta transverse part
            ddq = (displs[...,None,i,:] * qs).sum(-1) * (displs[...,js,:] * qs).sum(-1)
            # calculate the phi transverse part
            ddp = (displs[...,None,i,:-1] * ps).sum(-1) * (displs[...,js, :-1] * ps).sum(-1)
            # accumulate
            lrs = np.floor((np.log(dists)-lrmin)/lrbinsize).astype(int)
            bres[..., lrs] += [ddn, ddq, ddp, ddn**2, ddq**2, ddp**2]
            nb[lrs] += 1
        return bres, nb, np.exp(lrbins+lrmin)
    
    def MSD_twopoint(self, bres, nb):
        """Convert 'Brownian correlation' components (e.g. the output of get_Brownian_corr) into two-point mean squared displacement.
        
        Output :
            - r dependent two-point MSD components along r-r, theta-theta and phi-phi
            - the respective errors"""
        msds = bres[:3] / nb
        errors = np.sqrt(bres[3:] /nb - msds**2)
        return msds, errors
        
        
    def get_clusters(self, t, isNode):
        """Regroup the given particles into connected clusters.

Keyword arguments:
t -- time step where to look for clusters
isNode -- a 1d boolean array of dimension (number of trajectories) containing False where the particle should not be part of a cluster.
Single particle clusters are not returned, even if the particle was declared as node.

Return a dictionary (particle id -> cluster id)

"""
        nodes = np.where(isNode)
        posbonds = np.loadtxt(self.xp.get_format_string(ext='bonds')%t, dtype=int)
        pos2traj = -np.ones((posbonds.max()+1), dtype=int)
        pos2traj[self.trajs[nodes, t]] = nodes
        trajbonds = pos2traj[posbonds]
        trajbonds = trajbonds[np.where(trajbonds.min(axis=1)>-1)]
        gr = nx.Graph()
        gr.add_nodes_from(np.unique(trajbonds))
        gr.add_edges_from(trajbonds)
        return nx.connected_components(gr)

    def get_time_clusters(self, isNode):
        clusters = [self.get_clusters(t, isNode[t]) for t in range(self.xp.size)]
        timegr = nx.Graph()
        for t, cluster in enumerate(clusters):
            timegr.add_nodes_from(
        [(t,k) for k in np.unique(cluster.values())]
        )
        for t, frame in enumerate(clusters[:-1]):
                for tr, k in frame.iteritems():
                    c = clusters[t+1].get(tr,-1)
                    if c>-1 and not timegr.has_edge(((t, k), (t+1, c))):
                        timegr.add_edge(((t, k), (t+1, c)))
        return clusters, nx.connected_components(timegr)
    


def br(radius, T=28, eta28C=2.00139e-3, detadT=-0.03e-3):
        """Brownian time is the time for a particle to diffuse over it\'s own radius (in meters)"""
        return const.pi * (eta28C+(T-28)*detadT) * (radius**3) / (const.k * const.C2K(T))

def histz(f):
    """export the density histogram of a .dat file into a .z file"""
    hist, bins = np.histogram(np.loadtxt(f, skiprows=2, usecols=[2]))
    np.savetxt(f[:-3]+'z', hist/(bins[1]-bins[0]), fmt='%f', delimiter='\t')

def dilate(field, bonds):
    ngb = bonds2ngbs_list(bonds, field.shape[0])
    dil = np.zeros_like(field)
    for p,n in enumerate(ngb):
        a = field[n]
        dil[p] = a[np.absolute(a).argmax(axis=0), range(a.shape[1])]
    return dil

def average(field, bonds):
    ngb = bonds2ngbs_list(bonds, field.shape[0])
    av = np.zeros_like(field)
    for p,n in enumerate(ngb):
        av[p] = (field[p]+field[n].mean(axis=0))/2
    return av

OrnsteinZernike3D = lambda p, r: p[0]/r * np.exp(-r/p[1])
singleExponentialDecay = lambda p, t: p[0] * np.exp(-(t/p[1])**p[2])
multipleStrechedExp = lambda p, ts, ys: np.asarray([
    u for t, y, tau, beta in zip(ts, ys, p[1::2], p[2::2])
    for u in p[0] * np.exp(-(t/tau)**beta) - y
    ])
VogelFulcherTammann = lambda p, phi: p[0]*np.exp(p[1]*phi/(p[2]-phi))
errfun_msd = lambda D, t, y: np.log(6*D) + np.log(t) - np.log(y)

def fit_OrnsteinZernike(g6, envelope, p0=[1,1]):
    fitfunc = OrnsteinZernike3D
    errfunc = lambda p, x, y: fitfunc(p, x) - y
    logerrfunc = lambda p, x, y: np.log(fitfunc(p, x)) - np.log(y)
    p1, success = optimize.leastsq(errfunc, p0[:], args=(g6[:,0][envelope], g6[:,1][envelope]))
    p2, success = optimize.leastsq(logerrfunc, p0[:], args=(g6[:,0][envelope], g6[:,1][envelope]))
    return p1, p2

def fit_decay(isf, p0 = [0.95, 100, 1]):
    fitfunc = singleExponentialDecay
    errfunc = lambda p, x, y: fitfunc(p, x) - y
    if isf[0]==1:
        data = isf[1:]
    else:
        data = isf
    p1, success = optimize.leastsq(errfunc, p0, args=(np.arange(len(data))+1, data))
    return p1

def fit_multipleDecays(isfs, p0=None):
    """fit simultaneously many steched exponential decays having the same prefactor (plateau)
isfs - a list of 1D data to fit (not necessary the same length). The first value is not taken into account.
p0 - parameters as [prefactor, tau1, beta1, tau2, beta2, ...]"""
    if not p0:
        p0 = [0.95]+len(isfs)*[10, 1]
    p, success = optimize.leastsq(
        multipleStrechedExp, p0,
        ([np.arange(1, len(k)) for k in isfs], [k[1:] for k in isfs])
        )
    if not success==1: print("Fitting failure")
    return p

def fit_vft(tau, p0 = [30, 0.5, 0.62]):
    fitfunc = VogelFulcherTammann
    errfunc = lambda p, x, y: np.log(fitfunc(p, x)) - np.log(y)
    p1, success = optimize.leastsq(errfunc, p0[:], args=(tau[:,0], tau[:,1]))
    return p1

def fit_normaized_vft(tau, p0 = [0.5, 0.62]):
    errfunc = lambda p, phi, y: p[0]*phi/(p[1]-phi) - np.log(y)
    p1, success = optimize.leastsq(errfunc, p0[:], args=(tau[:,0], tau[:,1]))
    return p1

def normalize_rdf(a):
    """Removes the linear component of the unnormalized rdf by fitting its second half"""
    u, v = leastsq(
            lambda p, x, y: p[0]+p[1]*x-y,
            [40, 1],
            args=(a[len(a)/2:,0], a[len(a)/2:,1]))[0]
    b = np.copy(a)
    b[:,1] /= (u + v * a[:,0])
    return b

def envelope(g6, smooth=1.0):
    #smooth g6
    sg6 = gaussian_filter1d(np.copy(g6[:,1]),smooth)
    #find local maxima
    dg6 = np.gradient(sg6)
    env = np.where(np.bitwise_and(dg6[:-1]>0, dg6[1:]<0))[0]
    #remove the points before the first peak and with negative values
    env = env[g6[:,0][env]>1.5]
    env = env[g6[:,1][env]>0]
    #from peak to peak, the function should be decreasing
    #denv = np.gradient(sg6[env])
    #sel = list(np.where(denv[1:]<0)[0])+[len(env)-1]
    #env = env[sel]
    return env

def fit_cg6(infile, crop=None, remove_peaks=[], rmin=2, rmax=10):
    #a file.cg6 has 3 columns : r, N, N*G_6
    data = np.loadtxt(infile)
    #remove the end of the table if long time correlations are wrong
    if crop is not None:
        data = data[:crop]
    #remove the begining of the table with N=0
    blank = np.where(data[:,1]==0)[0]
    if len(blank)>0 and (blank[-1] < len(data)-1):
            data = data[blank[-1]+1:]
    #remove all values for r<1.5
    data = data[np.searchsorted(data[:,0],1):]
    #smooth the data
    cg6 = gaussian_filter1d(data[:,-1]/data[:,1], rmin)
    #high-pass filter to transform the function into an oscillatory decaying function
    #extract peaks from it
    peaks, mins = find_peak_mins(cg6*data[:,0]-gaussian_filter1d(cg6*data[:,0], rmax))
    #remove the peaks with negative values
    peaks = [p for p in peaks if cg6[p]>0]
    peaks = [p for i,p in enumerate(peaks) if i not in remove_peaks]
    params = leastsq(
        lambda p, x, y: np.log(p[0]/x*np.exp(-x/p[1]))-np.log(y),
        [2,3],
        args=(data[peaks,0], cg6[peaks])
        )[0]
    return params, np.column_stack((data[:,0], cg6)), peaks

def rdf2Sq(rdf, rho):
    """Calculate the radial Fourier transform of rdf(r) and normalize it to get the structure factor S(q)"""
    s = np.zeros_like(rdf)
    s[:,0] = np.fft.fftfreq(2*len(rdf)+1, rdf[1,0])[1:len(rdf)+1]
    s[:,1] = (rdf[:,0] * np.sin(np.outer(s[:,0], rdf[:,0])) * (rdf[:,1]-1)).sum(1) * (4*np.pi*rho)/s[:,0] * rdf[1,0]
    s[:,1] += 1
    return s

def find_peak_mins(a):
    """Find global maxima of decreasing intensity separated by a global minimum. Usefull for rdf"""
    peaks = [a.argmax()]
    mins = []
    while peaks[-1]<len(a)-1:
            mins.append(peaks[-1] + a[peaks[-1]:].argmin())
            if mins[-1]==len(a)-1:
                    break
            peaks.append(mins[-1] + a[mins[-1]:].argmax())
    return peaks, mins

def get_clusters(bonds):
    """Returns a list of clusters"""
    gr = nx.Graph()
    gr.add_nodes_from(np.unique1d(bonds.ravel()))
    gr.add_edges_from(bonds)
    return nx.connected_components(gr)


def shortest_path(g, trs):
    """From a graph g, whoes nodes are trajectory indices and edges corresponds to connectivity at time t1,
    and couples of nodes trs, first at time t0 second at time t1,
    yield the on-graph shortest path between each couple in term of trajectory indices. 
    
    trs is a iterable of couples of integers."""
    for i, (a,b) in enumerate(trs):
        if not g.has_node(a) or not g.has_node(b): continue
        try:
            yield nx.shortest_path(g, a,b)
        except nx.NetworkXNoPath:
            pass
            
def broken_bonds_path(bonds0, bonds1, p2tr0, p2tr1):
    """Paths on graph at t1 between particles involved in a bond at t0 and no more at t1.
    
    bonds0, bonds1 are respectively the bonds at t0 and 1 in terms of position
    p2tr0, p2tr1 are respectively the position to trajectory relationship at t0 and t1
    
    generate paths in terms of trajectory indices"""
    #bonds (between trajectories) existing at t but no more at t+dt
    # = broken bonds + lost trajectories
    trbonds = set([(a,b) for a,b in np.sort(p2tr0[bonds0], 1)]) - set([(a,b) for a,b in np.sort(p2tr1[bonds1], axis=1)])

    #graph of the bonds between trajectories at t+dt
    g = nx.Graph()
    g.add_nodes_from(p2tr1)
    g.add_edges_from(p2tr1[bonds1])
    return shortest_path(g, trbonds)
