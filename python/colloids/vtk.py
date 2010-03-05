from __future__ import with_statement #for python 2.5, useless in 2.6
import numpy as np

class Polydata:
    """A VTK polydata dataset"""

    def __init__(self, fileName):
        """Constructor from vtk legacy format"""
        self.scalars = {}
        self.vectors = {}
        with open(fileName) as f:
            f.seek(0, 2)
            filesize = f.tell()
            f.seek(0)
            while f.tell() < filesize:
                line = f.readline().split()
                if line[0] == 'POINTS':
                    size = int(line[1])
                    self.points = np.fromfile(
                        f, count=size*3, sep=" "
                        ).reshape((size,3))
                elif line[0] == 'LINES':
                    size = int(line[1])
                    self.bonds = np.fromfile(
                        f, dtype=int, count=size*3, sep=" "
                        ).reshape((size,3))[:,1:]
                elif line[0] == 'POINT_DATA':
                    size = int(line[1])
                elif line[0] == 'SCALARS':
                    f.readline() #discard LOOKUP_TABLE
                    self.scalars[line[1]] = np.fromfile(
                        f, count=size, sep=" "
                        )
                elif line[0] == 'VECTORS':
                    self.vectors[line[1]] = np.fromfile(
                        f, count=size*3, sep=" "
                        ).reshape((size,3))
            
                            
