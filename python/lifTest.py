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
import sys
from os.path import splitext,isfile
import lif
import Image

#reader = lif.LifReader("D:\Users\ishizaka\Documents\dataMathieu\Tsuru11dm_phi=52.53_J36.lif")
if len(sys.argv) >1:
    filename = sys.argv[1]
    if not isfile(filename):
        print "%s is not an existing file" % filename
        exit
else:
    while(True):
        filename = raw_input("filename --> ")
        if isfile(filename):
            break
        else:
            print "%s is not an existing file" % filename

reader = lif.LifReader(filename)

#export of the XML header
fHeader = open(splitext(filename)[0]+"_header.xml","w")
reader.xmlHeader.writexml(fHeader)
fHeader.close()

#Content short description
serie = reader.chooseSerie()
print "%s : %i frames X %i pixels per frames" % (serie.getName(),serie.getNbFrames(),serie.getNbPixelsPerFrame())

#just getting a 2D slice
serie.get2DImage().show()

#loading a full time step into a numpy array
a = serie.getFrame(0)
#extracting the slice Z==0 into an image
if len(serie.getDimensions()) >2:
    im = Image.fromarray(a[0,:,:])
else:
    im = Image.fromarray(a)
im.show()
