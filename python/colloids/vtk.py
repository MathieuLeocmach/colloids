from __future__ import with_statement #for python 2.5, useless in 2.6
import numpy as np

class Polydata:
    """A VTK polydata dataset"""

    def __init__(self, fileName):
        """Constructor from vtk legacy format"""
        self.scalars = []
        self.vectors = []
        self.name = ''
        with open(fileName) as f:
            f.seek(0, 2)
            filesize = f.tell()
            f.seek(0)
            while f.tell() < filesize:
                line = f.readline().split()
                if line[0] == '#':
                    self.name = f.readline()[:-1]
                elif line[0] == 'POINTS':
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
                    self.scalars.append((
                        line[1],np.fromfile(
                            f, count=size, sep=" "
                            )
                        ))
                elif line[0] == 'VECTORS':
                    self.vectors.append((
                        line[1],
                        np.fromfile(
                            f, count=size*3, sep=" "
                            ).reshape((size,3))
                        ))
            
    def save(self, fileName):
        """save to vtk legacy format"""
        with open(fileName, 'w') as f:
            f.write(
                ('# vtk DataFile Version 3.0\n%s\n' % self.name)+
                'ASCII\nDATASET POLYDATA\n'
                )
            
            f.write('POINTS %i double\n' % len(self.points))
            #self.points.tofile(f, sep=' ', format = '%f')
            for p in self.points:
                p.tofile(f, sep=' ', format = '%f')
                f.write('\n')

            size = len(self.bonds)
            f.write('LINES %i %i\n' % (size, 3*size))
            for b in self.bonds:
                f.write('2 ')
                b.tofile(f, sep=' ', format = '%i')
                f.write('\n')
            #np.hstack((
             #   np.ones((size,1), dtype=int)*2,
              #  self.bonds
               # )).tofile(f, sep=' ', format = '%i')

            f.write('POINT_DATA %i\n' % len(self.points))

            for name, field in self.scalars:
                f.write('SCALARS %s double\n'%name)
                f.write('LOOKUP_TABLE default\n')
                field.tofile(f, sep='\n', format = '%f')
            for name, field in self.vectors:
                f.write('VECTORS %s double\n'%name)
                #field.tofile(f, sep=' ', format = '%f')
                for v in field:
                    v.tofile(f, sep=' ', format = '%f')
                    f.write('\n')


