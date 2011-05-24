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

    def link(self, other, maxdist):
        """Try to find the best correpondence between the particles of the present object and those of another one"""
        #list the possible links
        links = [
            (i, j, self.get_sqdist(i, q=other.position[j], r=other.position[j]))
            for j in other.index.intersection(np.concatenate((p-r,p+r)))
            for i, (p,r) in enumerate(zip(
                self.pos, maxdist*(self.radii+other.maxRad())
                ))
            ]
        #keep only the links shorter than maxdist
        maxsqd = maxdist**2
        links = np.asarray(
            [(i,j,l) for i,j,l in links if l<maxsqd],
            dtype=np.dtype("i4, i4, f8")
            )
        links = links[np.argsort(links['f2'])]
        usedi = np.zeros(len(self.pos), bool)
        usedj = np.zeros(len(other.pos), bool)
        goodlinks = []
        for i,j,l in links:
            if usedi[i] or usedj[j]:
                continue
            goodlinks.append((i,j))
            usedi[i] = True
            usedj[j] = True
        return np.asarray(goodlinks)
        
def non_overlapping(positions, radii):
    """Give the indices of non-overlapping particles. Early bird."""
    assert len(positions)==len(radii)
    pr = rtree.index.Property()
    pr.dimension = positions.shape[1]
    tree = rtree.index.Index(properties=pr, interleaved=True)
    for i,(p,r) in enumerate(zip(positions, radii)):
        bb = np.concatenate((p-r, p+r))
        for j in tree.intersection(bb):
            if np.sum((p-positions[j])**2) < (r + radii[j])**2:
                break
        else:
            tree.insert(i, bb)
    return [i for i in tree.intersection(tree.bounds)]
