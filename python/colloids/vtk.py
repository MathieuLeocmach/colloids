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
import numpy as np
from scipy.spatial import KDTree
import math

class Polydata:
    """A VTK polydata dataset"""

    def __init__(self, fileName=None):
        """Constructor from vtk legacy format"""
        self.name = ''
        self.points = np.empty((0,3))
        self.bonds = np.empty((0,2), dtype=int)
        self.scalars = []
        self.vectors = []
        self.bondsScalars = []
        if not fileName ==None:
            self.load(fileName)

    def load(self, fileName=None):
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
        self.__fields = dict(self.scalars + self.vectors)
            
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
                p.tofile(f, sep=' ', format = '%g')
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
                if field.dtype.kind=='i':
                    f.write('SCALARS %s int\n'%name)
                    f.write('LOOKUP_TABLE default\n')
                    field.tofile(f, sep='\n', format = '%i')
                else:
                    f.write('SCALARS %s double\n'%name)
                    f.write('LOOKUP_TABLE default\n')
                    field.tofile(f, sep='\n', format = '%g')
                f.write('\n')
            for name, field in self.vectors:
                f.write('VECTORS %s double\n'%name)
                #field.tofile(f, sep=' ', format = '%f')
                for v in field:
                    v.tofile(f, sep=' ', format = '%g')
                    f.write('\n')

            if len(self.bondsScalars)>0:
                f.write('CELL_DATA %i\n' % len(self.bonds))
                for name, field in self.bondsScalars:
                    if field.dtype.kind=='i':
                        f.write('SCALARS %s int\n'%name)
                        f.write('LOOKUP_TABLE default\n')
                        field.tofile(f, sep='\n', format = '%i')
                    else:
                        f.write('SCALARS %s double\n'%name)
                        f.write('LOOKUP_TABLE default\n')
                        field.tofile(f, sep='\n', format = '%g')
                    f.write('\n')

    def getField(self, name):
        """Return the field corresponding to name"""
        self.__fields = dict(self.scalars + self.vectors)
        return self.__fields[name]

    def coarseGrainField(self, field):
        """coarse grain the field by averaging the value at p over the neighbourhood of p"""
        cg = np.copy(field)
        cg[self.bonds[:,0]] += field[self.bonds[:,1]]
        cg[self.bonds[:,1]] += field[self.bonds[:,0]]
        div = np.ones_like(field)
        div[:self.bonds[:,0].max()+1] += np.bincount(self.bonds[:,0])
        div[:self.bonds[:,1].max()+1] += np.bincount(self.bonds[:,1])
        return cg/div

    def dilateField(self, field):
        """coarse grain the field by averaging the value at p over the neighbourhood of p"""
        dil = np.zeros_like(field)
        for b in self.bonds:
            dil[b[0]] = dil[b[1]] = max([field[b[0]], field[b[1]]])
        return dil

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

def export_structured_points(fname, data, name="", spacing=np.ones(3), origin=np.zeros(3)):
    """Export grid data to a vtk file"""
    shape = list(data.shape)
    shape.reverse()
    with open(fname, 'wb') as f:
        f.write(
            ('# vtk DataFile Version 3.0\n%s\n' % name)+
            'BINARY\nDATASET STRUCTURED_POINTS\n'+
            ('DIMENSIONS %d %d %d\n'%tuple(shape))+
            ('ORIGIN %g %g %g\n'%tuple(list(origin)))+
            ('SPACING %g %g %g\n'%tuple(list(spacing)))+
            ('POINT_DATA %d\n'%data.size)+
            'SCALARS Intensity unsigned_char\nLOOKUP_TABLE default\n'
            )
        if data.dtype.itemsize > 8:
             #Paraview reads only bigendian by default
            np.array(data, np.dtype('>'+data.dtype.str[1:])).tofile(f)
        else:
            data.tofile(f)
        

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
    #parameters parsing
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
    #spatial query
    lowerBound = points.min(axis=0) + maxDist/2
    upperBound = points.max(axis=0) - maxDist/2
    inside_id = [
        i for i, p in enumerate(points)
        if (p>= lowerBound).all() and (p <= upperBound).all()
        ]
    tree = KDTree(points)
    inside_tree = KDTree(points[inside_id])
    pairs = inside_tree.query_ball_tree(tree, maxDist)
    #binning
    coord_bins = np.zeros((Nbins, fields.shape[1]))
    nb_bins = np.zeros((Nbins), dtype=int)
    for p, qs in zip(inside_id, pairs):
        qs.remove(p)
        rs = np.asarray(
            np.sqrt(
                ((points[qs] - points[p])**2).sum(axis=1)
                ) * Nbins / maxDist,
            dtype=int)
        nb_bins[rs] += 1 
        coord_bins[rs] += fields[qs]*fields[p]
    bins = np.column_stack([coord_bins[:,cols].sum(axis=1) for cols in slices])
    bins[np.nonzero(nb_bins)] /= nb_bins[np.nonzero(nb_bins)][:,np.newaxis]
    return np.column_stack((np.arange(Nbins, dtype=float)/maxDist,bins))


