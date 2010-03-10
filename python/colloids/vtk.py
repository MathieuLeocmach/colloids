from __future__ import with_statement #for python 2.5, useless in 2.6
import numpy as np
from scipy.spatial import KDTree
import math

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

    def coarseGrainField(self, field):
        """coarse grain the field by averaging the value at p over the neighbourhood of p"""
        cg = np.copy(field)
        cg[self.bonds[:,0]] += field[self.bonds[:,1]]
        cg[self.bonds[:,1]] += field[self.bonds[:,0]]
        div = np.ones_like(field)
        div[:self.bonds[:,0].max()+1] += np.bincount(self.bonds[:,0])
        div[:self.bonds[:,1].max()+1] += np.bincount(self.bonds[:,1])
        return cg/div

    def spatialCorelation(self, Nbins=200, maxDist=50.0):
        """Compute the spatial corellation of each field"""
        return spatialCorelation(
            self.points,
            np.hstack(tuple(
                [np.column_stack(tuple(
                    [s for n, s in self.scalars]
                    ))]+[v for n, v in self.vectors]
                      )),
            np.array(
                range(len(self.scalars))
                + [i/3 + len(self.scalars) for i in range(3*len(self.vectors))]
                ),
            Nbins,
            maxDist
            )
        

def spatialCorelation(points, fields, vectorColumns=None, Nbins=200, maxDist=50.0):
    """Compute the spatial corellation of each field

    points -- 2D array of points coordinates. Shape is (N,d) with d the number of spatial dimensions.
    fields -- 2D array of scalar field or of coordinates of vector fields. Shape is (N, F)
    with F the sum of the dimensions of each field.
    vectorColumns -- 1D array indexing the columns of fields into vector fields.
    for example [0, 1, 1, 1] means that the first column of fields is the scalar field 0 and
    the next 3 columns are the coordinates of a 3D vector field.
    Nbins -- The number of bins of the histogram
    maxLength -- The maximum distance between a pair of points taken into account in the histogram
    """
    if len(points) != len(fields):
        raise ValueError(
            'You must have exactly one field value per point\n'
            + 'Here points id %i and fieds is %i'%(len(points), len(fields))
            )
    if vectorColumns==None:
        vectorColumns = np.arange(fields.shape[1])
    if len(vectorColumns) != fields.shape[1]:
        vectorColumns = np.concatenate((
            vectorColumns,
            np.arange(vectorColumns.max()+1,fields.shape[1])
            ))
    slices = [np.where(vectorColumns==v)[0] for v in range(vectorColumns.max()+1)]
    bins = np.zeros((Nbins, 2 + vectorColumns.max()))
    lowerBound = points.min(axis=0) + maxDist/2
    upperBound = points.max(axis=0) - maxDist/2
    maxDistSq = maxDist**2
    for p, f in zip(points, fields):
        if (p >= lowerBound).min() and (p <= upperBound).min():
            for q, g in zip(points, fields):
                dsq = np.dot(p, q)
                if dsq < maxDistSq:
                    r = math.sqrt(dsq) * Nbins / maxDist
                    toAdd = [1.0] + [
                        np.dot(f[cols], g[cols])
                        for cols in slices
                                     ]
                    bins[r] += np.asarray(toAdd)
    bins[:,1:] /= bins[:,0:1]
    bins[:,0] = np.arange(Nbins)/maxDist
    return bins


