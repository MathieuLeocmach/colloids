#
#    Copyright 2011 Mathieu Leocmach
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
import numpy as np
import rtree.index

class Particles:
    """Positions of the particles at a givent time"""
    def __init__(self, positions, radius=0.0):
        self.pos = positions
        if np.isscalar(radius):
            self.radii = np.ones(len(self.pos))*radius
        else:
            if radius.ndim>1 or len(radius)!=len(self.pos):
                raise ValueError("""radius must be a one dimensional array with the same length as positions""")
            self.radii = radius
        self.indexing()
        self.__maxRad = self.radii.max()

    def generated_boundingbox(self):
        for i, bb in enumerate(
            np.column_stack((
                (self.pos.T-self.radii).T,
                (self.pos.T+self.radii).T))
            ):
            yield i, bb, None
        
    def indexing(self):
        p = rtree.index.Property()
        p.dimension = self.pos.shape[1]
        self.index = rtree.index.Index(
            self.generated_boundingbox(),
            properties=p,
            interleaved=True)

    def maxRad(self):
        return self.__maxRad

    def get_sqdist(self, i, j=None, q=None, r=None):
        if q is None or r is None:
            q = self.pos[j]
            r = self.radii[j]
        return np.sum((self.pos[i]-q)**2, axis=-1)/(self.radii[i]+r)**2

    def get_ngbs(self, i, maxlength):
        rnge = self.maxRad()*(maxlength-1)+self.radii[i]*maxlength
        ngb1 = np.asarray([j for j in self.index.intersection(np.concatenate((
            self.pos[i]-rnge,
            self.pos[i]+rnge
            ))) if j!=i])
        return ngb1[self.get_sqdist(i, ngb1)<maxlength**2]

    def get_bonds(self, maxlength):
        return np.vstack([
            [i, n] for i in range(len(self.pos))
            for n in self.get_ngbs(i, maxlength)
            if n>i])
        
