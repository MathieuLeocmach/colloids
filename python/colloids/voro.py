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

def voroVolume(fname):
    """Use voro++ to output to disk the volume of the voronoi cell of each particle in file.vol"""
    outName = splitext(fname)[0]
    with open(fname) as f:
        f.readline()
        bb = [float(x)-6 for x in f.readline()[:-1].split()]
        with open(outName, 'w') as out:
            for i, p in enumerate(f):
                out.write('%i %s' % (i, p))

    subprocess.check_call(
        shlex.split(
            'voro++ -c "%i %v" '+('10 6 %f 6 %f 6 %f' % tuple(bb))
        )+[outName]
        )
    vol = np.sort(
            np.loadtxt(
                outName+".vol",
                dtype=[('id', int),('V', float)]
                ),
            order='id'
            )['V']
    np.savetxt(outName+'.vol', vol, fmt='%f')
    os.remove(outName)
    return vol

def voroFaces(fname):
    """Use voro++ to output in file.face the area of each face of each voronoi cell

    Output format is:
    ID_particle Nb_neighbours List_neigbours_ID List_face_areas

    Neighbours with negative ID are walls
    """
    outName = splitext(fname)[0]+'.face'
    with open(fname) as f:
        f.readline()
        bb = [float(x)-6 for x in f.readline()[:-1].split()]
        with open(outName, 'w') as out:
            for i, p in enumerate(f):
                out.write('%i %s' % (i, p))

    subprocess.check_call(
        shlex.split(
            'voro++ -c "%i %s %n %f" '+('10 6 %f 6 %f 6 %f' % tuple(bb))
        )+[outName]
        )
    os.remove(outName)
    os.rename(outName+'.vol', outName)
