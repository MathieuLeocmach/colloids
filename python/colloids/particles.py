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
from scipy.special import sph_harm
import numexpr



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

    def inside_box(self, margin):
        bb = np.asarray(self.index.bounds)
        bb[:len(bb)/2] += margin
        bb[len(bb)/2:] -= margin
        return bb

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

    def get_N_ngbs(self, i, maxlength, N=12):
        """Get the first Nth neighbours of particle i"""
        rnge = self.maxRad()*(maxlength-1)+self.radii[i]*maxlength
        ngb1 = np.asarray([j for j in self.index.intersection(np.concatenate((
            self.pos[i]-rnge,
            self.pos[i]+rnge
            ))) if j!=i])
        assert len(ngb1)>=N, "Only %d<12 neighbours within maxlength for particle %d, increase maxlength"%(len(ngb1), i)
        distsq = self.get_sqdist(i, ngb1)
        sorted_by_distsq = np.argsort(distsq)
        assert distsq[sorted_by_distsq[N-1]]<maxlength**2, "Only %d<12 neighbours within maxlength for particle %d, increase maxlength"%(np.sum(distsq<maxlength**2), i)
        return ngb1[np.argsort(distsq)[:N]]    

    def get_bonds(self, maxlength):
        return np.vstack([
            [i, n] for i in range(len(self.pos))
            for n in self.get_ngbs(i, maxlength)
            if n>i])

    def scaled_rdf(self, maxlength, nbins=200):
        """Construct the rdf but scales each distance by the sum of the radii of the pair"""
        bins = np.linspace(0, maxlength, nbins)
        margin = 2 * self.maxRad() * maxlength
        bb = self.inside_box(margin)
        inside = [i for i in self.index.intersection(bb)]
        bonds = np.vstack([
            [i, n] for i in inside
            for n in self.get_ngbs(i, maxlength)
            if n!=i])
        bondlengths = np.sqrt(np.sum(np.diff(self.pos[bonds], axis=1)[:,0]**2,-1))
        bondlengths /= self.radii[bonds].sum(-1)
        nbdens = (self.radii[inside]**3).sum() * 8.0 / self.pos[inside].ptp(0).prod()
        g = np.histogram(bondlengths, bins=bins)[0] / (
            4*np.pi * bins[1:]**2 * bins[1] * nbdens * (len(np.unique1d(bonds[:,0]))-1)
            )
        return g

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
    
    def get_qlm(self, i, ngb_indices, ls=[6]):
        #vectors to neighbours
        cartesian = self.pos[ngb_indices]-self.pos[i]
        spherical = cart2sph(cartesian)
        return [
            np.asarray([
                sph_harm(m, l, spherical[:,2], spherical[:,1]).mean()
                for m in range(l+1)])
            for l in ls
            ]
        
        
        
        
        
        
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

def cart2sph(cartesian):
    """Convert Cartesian coordinates [[x,y,z],] to spherical coordinates [[r,phi,theta],]
phi is cologitudinal and theta azimutal"""
    spherical = np.zeros_like(cartesian)
    #distances
    spherical[:,0] = np.sum(cartesian**2, axis=-1)
    sel = np.asarray([
        r!=z or r+1.0>1.0
        for z, r in zip(cartesian[:,-1]**2, spherical[:,0])
        ])
    spherical[:,0] = np.sqrt(spherical[:,0])
    #colatitudinal phi [0, pi[
    spherical[sel,1] = np.arccos(cartesian[sel,-1]/spherical[sel,0])
    #azimutal (longitudinal) theta [0, 2pi[
    spherical[sel,2] = np.arctan2(cartesian[sel,1], cartesian[sel,0])
    spherical[spherical[:,2]<0, 2] += 2*np.pi
    return spherical
    

def ql(qlm):
    l = qlm.shape[1]-1
    n = numexpr.evaluate("""abs(qlm).real**2""")
    return numexpr.evaluate(
        """sqrt(4*pi/(2*l+1)*(2*a+b))""",
        {'a':np.sum(n[:,1:], -1), 'b':n[:,0], 'l':qlm.shape[1]-1, 'pi':np.pi}
        )
##factor 10 speed using numpexr
##    s = 2*np.sum(np.abs(qlm[:,1:])**2, -1) + np.abs(qlm[:,0])**2
##    return np.sqrt(4*np.pi/(2*l+1)*s)
def wl(qlm):
    l = qlm.shape[1]-1
    w = np.zeros(qlm.shape[0])
    for m1 in range(-l, l+1):
        qlm1 = get_qlm(qlm, m1)
        for m2 in range(-l, l+1):
            m3 = -m1-m2
            if -l<=m3 and m3<=l:
                w += numexpr.evaluate(
                    """real(w3j * qlm1 * qlm2 * qlm3)""",
                    {
                        'w3j': get_w3j(l, [m1, m2, m3]),
                        'qlm1': qlm1,
                        'qlm2': get_qlm(qlm, m2),
                        'qlm3': get_qlm(qlm, m3)
                        })
                #factor 2 speed using numexpr
##                w += get_w3j(l, [m1, m2, m3]) * get_qlm(qlm, m1) * get_qlm(qlm, m2) * get_qlm(qlm, m3)
    return w
    
def get_qlm(qlms, m):
    if m>=0:
        return qlms[:,m]
    if (-m)%2 == 0:
        return np.conj(qlms[:,-m])
    return -np.conj(qlms[:,-m])
		
_w3j = [
    [1],
    np.sqrt([2/35., 1/70., 2/35., 3/35.])*[-1,1,1,-1],
    np.sqrt([
        2/1001., 1/2002., 11/182., 5/1001.,
        7/286., 5/143., 14/143., 35/143., 5/143.,
        ])*[3, -3, -1/3.0, 2, 1, -1/3.0, 1/3.0, -1/3.0, 1],
    np.sqrt([
        1/46189., 1/46189.,
        11/4199., 105/46189.,
        1/46189., 21/92378.,
        1/46189., 35/46189., 14/46189.,
        11/4199., 21/4199., 7/4199.,
        11/4199., 77/8398., 70/4199., 21/4199.
        ])*[-20, 10, 1, -2, -43/2.0, 3, 4, 2.5, -6, 2.5, -1.5, 1, 1, -1, 1, -2],
    np.sqrt([
        10/96577., 5/193154.,
        1/965770., 14/96577.,
        1/965770., 66/482885.,
        5/193154., 3/96577., 77/482885.,
        65/14858., 5/7429., 42/37145.,
        65/14858., 0.0, 3/7429., 66/37145.,
        13/74290., 78/7429., 26/37145., 33/37145.,
        26/37145., 13/37145., 273/37145., 429/37145., 11/7429.,
        ])*[
            7, -7, -37, 6, 73, -3,
            -5, -8, 6, -1, 3, -1,
            1, 0, -3, 2, 7, -1, 3, -1,
            1, -3, 1, -1, 3],
    np.sqrt([
        7/33393355., 7/33393355.,
        7/33393355., 462/6678671.,
        7/33393355., 1001/6678671.,
        1/233753485., 77/6678671., 6006/6678671.,
        1/233753485., 55/46750697., 1155/13357342.,
        1/233753485., 2926/1757545., 33/46750697., 3003/6678671.,
        119/1964315., 22/2750041., 1914/474145., 429/5500082.,
        17/13750205., 561/2750041., 77/392863., 143/27500410., 2002/392863.,
        323/723695., 1309/20677., 374/144739., 143/144739., 1001/206770.,
        323/723695., 7106/723695., 561/723695., 2431/723695., 2002/103385., 1001/103385.
        ])*[
            -126, 63, 196/3.0, -7, -259/2.0, 7/3.0,
            1097/3.0, 59/6.0, -2,
            4021/6.0, -113/2.0, 3,
            -914, 1/3.0, 48, -3,
            -7/3.0, 65/3.0, -1, 3,
            214/3.0, -3, -2/3.0, 71/3.0, -1,
            3, -1/3.0, 5/3.0, -2, 1/3.0,
            2/3.0, -1/3.0, 2, -4/3.0, 2/3.0, -1]
    ]
_w3j_m1_offset = np.array([0,1,2,4,6,9,12,16,20,25,30], int)

def get_w3j(l, ms):
    sm = np.sort(np.abs(ms))
    return _w3j[l/2][_w3j_m1_offset[sm[-1]]+sm[0]]
