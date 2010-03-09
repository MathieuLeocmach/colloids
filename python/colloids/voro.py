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
    with open(fname, 'w') as f:
        for i, p in enumerate(a):
            f.write('%i %s\n' % (i, np.array_str(p)[1:-1]))

def volume(fname):
    """Use voro++ to output to disk the volume of the voronoi cell of each particle in file.vol"""
    outName = splitext(fname)[0]
    with open(fname) as f:
        f.readline()
        bb = [float(x)-5 for x in f.readline()[:-1].split()]
        with open(outName, 'w') as out:
            for i, p in enumerate(f):
                out.write('%i %s' % (i, p))

    subprocess.check_call(
        shlex.split(
            'voro++ -c "%i %v" '+('10 5 %g 5 %g 5 %g' % tuple(bb))
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
    with open(fname) as f:
        f.readline()
        bb = [float(x)-5 for x in f.readline()[:-1].split()]
        with open(outName, 'w') as out:
            for i, p in enumerate(f):
                out.write('%i %s' % (i, p))

    subprocess.check_call(
        shlex.split(
            'voro++ -c "%i %s %n %g" '+('10 5 %g 5 %g 5 %g' % tuple(bb))
        )+[outName]
        )
    os.remove(outName)
    os.rename(outName+'.vol', outName)

def cgVolume(fname):
    """Use voro++ to output to disk the volume of the voronoi cell of each particle in file.vol"""
    outName = '_space_t'.join(splitext(fname)[0].split('_t'))
    with open(fname) as f:
        f.readline()
        bb = [float(x)-5 for x in f.readline()[:-1].split()]
        with open(outName, 'w') as out:
            for i, p in enumerate(f):
                out.write('%i %s' % (i, p))
            size = i+1

    subprocess.check_call(
        shlex.split(
            'voro++ -c "%i %v %n" '+('10 5 %g 5 %g 5 %g' % tuple(bb))
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
            nbNgb[i] = len(ngb)
            if len(ngb)+2 == len(l):
                #the center cell volume itself is counted only if not near a wall
                vol[i] += v
                nbNgb[i] += 1
            for n in ngb:
                vol[n] += v
    vol /= nbNgb
    np.savetxt(outName+'.vol', vol, fmt='%g')
    os.remove(outName)
