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
from numpy import savetxt

if len(sys.argv) >1:
    filename = sys.argv[1]
else:
    while(True):
        filename = raw_input("filename --> ")
        if isfile(filename):
            break
        else:
            print '"%s" is not an existing file' % filename

reader = lif.LifReader(filename)
for serie in reader.getSeries():
    if serie.getNbFrames() > 1:
        displ = serie.getDisplacements2D()
        savetxt(
            splitext(filename)[0]+"_"+serie.getName()+".displ",
            displ,
            fmt='%1d',
            delimiter='\t'
            )
