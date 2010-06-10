#! /usr/bin/env python
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
import lif
import sys
from os.path import splitext,isfile

if len(sys.argv) >1:
    filename = sys.argv[1]
else:
    while(True):
        filename = raw_input("filename --> ")
        if isfile(filename):
            break
        else:
            print '"%s" is not an existing file' % filename

reader = lif.Reader(filename)
hout = open(splitext(filename)[0]+'.xml', 'w')
reader.xmlHeader.writexml(hout)
hout.close()

if len(sys.argv)>2:
    serie = reader.chooseSerie(int(sys.argv[1]))
else:
    serie = reader.chooseSerie()
#geometrical information
shape = serie.getFrameShape()
extent = [serie.getVoxelSize(d+1)*l for d,l in enumerate(shape)]
flatdims = 3-len(shape)
nbpixels = serie.getNbPixelsPerFrame()
nbchannels = max(1, len(serie.getChannels()))

nbdigits = len('%s'%serie.getNbFrames())
output = u'%s_%s_t%%0%dd.vtk' %(splitext(filename)[0], serie.getName(), nbdigits)

serie.f.seek(serie.getOffset())
for t in range(serie.getNbFrames()):
    out = open(output%t, 'wb')
    out.write('# vtk DataFile Version 2.0\n')
    out.write('%s %s\n' % (reader.getName(), serie.getName()))
    out.write('BINARY\nDATASET STRUCTURED_POINTS\n')
    out.write('DIMENSIONS %d %d %d\nORIGIN 0 0 0\n' % tuple(shape+flatdims*[1]))
    out.write(
        'SPACING %g %g %g\n' % tuple(
            [serie.getVoxelSize(d+1) for d in range(len(shape))]+flatdims*[1]
            )
        )
    out.write('POINT_DATA %d\n' % nbpixels)
    #out.write('SCALARS intensity unsigned_char %d\n' % nbchannels)
    #out.write('LOOKUP_TABLE default\n')
    out.write('COLOR_SCALARS intensity %d\n' % nbchannels)
    out.write(serie.f.read(nbpixels*nbchannels))
    out.close()
