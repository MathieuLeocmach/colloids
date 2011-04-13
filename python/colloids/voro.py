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
from __future__ import with_statement
import numpy as np
import shlex, subprocess, os
from os.path import splitext

def saveVoroFormat(fname, a):
    """Save a numpy array to a format readable to voro++"""
    np.savetxt(
        fname,
        np.column_stack((np.arange(len(a)), a)),
        fmt='%d %g %g %g'
        )

def volume(fname):
    """Use voro++ to output to disk the volume of the voronoi cell of each particle in file.vol"""
    outName = splitext(fname)[0]
    dat = np.loadtxt(fname, skiprows=2)
    bb = np.column_stack([dat.min(axis=0)-0.1, dat.max(axis=0)+0.1])
    saveVoroFormat(outName, dat)
    subprocess.check_call(
        shlex.split(
            'voro++ -c "%i %v" '+('10 %g %g %g %g %g %g' % tuple(bb.ravel()))
            )+[outName]
        )
    vol = np.sort(
            np.loadtxt(
                outName+".vol",
                dtype=[('id', int),('V', float)]
                ),
            order='id'
            )['V']
    np.savetxt(outName+'.vol', vol, fmt='%g')
    os.remove(outName)
    return vol

def faces(fname):
    """Use voro++ to output in file.face the area of each face of each voronoi cell

    Output format is:
    ID_particle Nb_neighbours List_neigbours_ID List_face_areas

    Neighbours with negative ID are walls
    """
    outName = splitext(fname)[0]+'.face'
    dat = np.loadtxt(fname, skiprows=2)
    bb = np.column_stack([dat.min(axis=0)-0.1, dat.max(axis=0)+0.1])
    saveVoroFormat(outName, dat)

    subprocess.check_call(
        shlex.split(
            'voro++ -c "%i %s %n %g" '+('10 5 %g 5 %g 5 %g' % tuple(bb))
        )+[outName]
        )
    os.remove(outName)
    os.rename(outName+'.vol', outName)

def cgVolume(fname, Return=False):
    """Use voro++ to output to disk the volume of the voronoi cell of each particle in file.vol"""
    outName = '_space_t'.join(splitext(fname)[0].split('_t'))
    dat = np.loadtxt(fname, skiprows=2)
    size = len(dat)
    bb = np.column_stack([dat.min(axis=0)-0.1, dat.max(axis=0)+0.1])
    saveVoroFormat(outName, dat)
    subprocess.check_call(
        shlex.split(
            'voro++ -c "%i %v %n" '+('10 %g %g %g %g %g %g' % tuple(bb.ravel()))
            )+[outName]
        )
    vol = np.zeros(size)
    nbNgb = np.zeros(size, dtype=int)
    with open(outName+'.vol') as f:
        for line in f:
            l = line.split()
            i = int(l[0])
            v = float(l[1])
            ngb = [n for n in map(int, l[2:]) if n>=0]
            #the center cell volume is significant only if not near a wall
            if len(ngb)+2 == len(l):
                vol[i] += v
                nbNgb[i] += 1
                for n in ngb:
                    vol[n] += v
                    nbNgb[n] += 1
    vol[np.nonzero(vol*nbNgb)] /= nbNgb[np.nonzero(vol*nbNgb)]
    np.savetxt(outName+'.vol', vol, fmt='%g')
    os.remove(outName)
    if Return :
        return np.column_stack((dat,vol))

def load_vorobonds(fname):
    """load the bond network from a custom output of Voro++ '%i %v %s %n'"""
    #load all bonds
    bonds = np.vstack([np.column_stack((
	int(line.split()[0])*np.ones(int(line.split()[2]), int),
	map(int, line.split()[3:])
	)) for line in open(fname)])
    walls = np.signbit(bonds.min(axis=-1))
    outside = np.unique1d(bonds[walls].max(axis=-1))
    #remove the walls and the duplicates
    bonds = bonds[np.bitwise_and(
	np.diff(bonds, axis=-1)[:,0]>0,
	np.bitwise_not(walls)
	)]
    #sort by second then first column
    return bonds[np.lexsort(bonds.T[::-1].tolist())], outside
