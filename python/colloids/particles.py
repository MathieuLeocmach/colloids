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
import numexpr
import subprocess, shlex, StringIO
from scipy import weave
from scipy.weave import converters


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
        self.__index = rtree.index.Index(
            self.generated_boundingbox(),
            properties=p,
            interleaved=True)

    def get_index(self):
        if not hasattr(self, "__index"):
            self.indexing()
        return self.__index

    def maxRad(self):
        return self.__maxRad

    def inside_box(self, margin):
        bb = np.asarray(self.get_index().bounds)
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
        ngb1 = np.asarray([j for j in self.get_index().intersection(np.concatenate((
            self.pos[i]-rnge,
            self.pos[i]+rnge
            ))) if j!=i])
        return ngb1[self.get_sqdist(i, ngb1)<maxlength**2]

    def get_N_ngbs(self, i, maxlength, N=12):
        """Get the first Nth neighbours of particle i"""
        rnge = self.maxRad()*(maxlength-1)+self.radii[i]*maxlength
        ngb1 = np.asarray([j for j in self.get_index().intersection(np.concatenate((
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
        inside = [i for i in self.get_index().intersection(bb)]
        bonds = np.vstack([
            [i, n] for i in inside
            for n in self.get_ngbs(i, maxlength)
            if n!=i])
        bondlengths = np.sqrt(numexpr.evaluate(
                                               '(a-b)**2', 
                                               {'a': self.pos[bonds[:,0]], 'b': self.pos[bonds[:,1]]}
                                               ).sum(-1))
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
            for j in other.get_index().intersection(np.concatenate((p-r,p+r)))
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

    def get_voro(self):
        """Interface to the stream version of voro++"""
        bb = np.asarray(self.get_index().bounds)
        pr = subprocess.Popen(
            [
                'voro++', '-r', '-c', '"%i %v %n"', '%g'%(self.maxRad()*2)
                ]+['%g'%b for b in bb.reshape([2,3]).T.ravel()]+['-'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE
        )
        out = pr.communicate(''.join([
            '%d %f %f %f %f '%(i, p[0], p[1], p[2], np.sqrt(2)*r)
            for i, (p, r) in enumerate(zip(self.pos, self.radii))
            ]))[0]
        vol = np.zeros(len(self.pos))
        ngbs = [[] for i in range(len(self.pos))]
        for line in out.split("\n"):
            if len(line)==0:
                continue
            l = line[1:-1].split()
            i = int(l[0])
            vol[i] = float(l[1])
            ngbs[i] = map(int, l[2:])
        return vol, ngbs
        
        
        
        
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

def weave_non_overlapping(positions, radii):
    """Give the mask of non-overlapping particles. Early bird."""
    assert len(positions)==len(radii)
    support = """
    typedef RStarTree<int, 3, 4, 32, double> RTree;
	struct Gatherer {
		std::list<int> *gathered;
		bool ContinueVisiting;

		Gatherer(std::list<int> &result) : gathered(&result), ContinueVisiting(true) {};

		void operator()(const typename RTree::Leaf * const leaf)
		{
			gathered->push_back(leaf->leaf);
		}
	};
    """
    good = np.zeros(len(positions), dtype=bool)
    code = """
    RTree tree;
    for(int p=0; p<Npositions[0]; ++p)
    {
        //norm 1 overlapping
		typename RTree::BoundingBox bb;
		for(int d=0; d<3; ++d)
		{
			bb.edges[d].first = positions(p,d) - radii(p);
			bb.edges[d].second = positions(p,d) + radii(p);
		}
		std::list<int> overlapping;
		tree.Query(typename RTree::AcceptOverlapping(bb), Gatherer(overlapping));
		bool is_overlapping = false;
		//norm 2 overlapping
		for(std::list<int>::const_iterator q= overlapping.begin(); q!=overlapping.end(); ++q)
		{
		    double disq = 0;
		    for(int d=0; d<3; ++d)
		        disq += pow(positions(p,d)-positions(*q,d), 2);
			if(disq < pow(radii(p) + radii(*q) ,2))
			{
				is_overlapping = true;
				break;
			}
		}
		if(!is_overlapping)
		{
			tree.Insert(p, bb);
			good(p) = true;
		}
    }
    """
    weave.inline(
        code,['positions', 'radii', 'good'],
        type_converters =converters.blitz,
        support_code = support,
        include_dirs = ['/home/mathieu/src/colloids/multiscale/RStarTree'],
        headers = ['"RStarTree.h"'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return good

def get_bonds(positions, radii, maxdist=3.0):
    support = """
    typedef RStarTree<int, 3, 4, 32, double> RTree;
	struct Gatherer {
		std::list<int> *gathered;
		bool ContinueVisiting;

		Gatherer(std::list<int> &result) : gathered(&result), ContinueVisiting(true) {};

		void operator()(const typename RTree::Leaf * const leaf)
		{
			gathered->push_back(leaf->leaf);
		}
	};
    """
    pairs = []
    dists = []
    code = """
    //spatial indexing
    RTree tree;
    for(int p=0; p<Npositions[0]; ++p)
    {
        typename RTree::BoundingBox bb;
		for(int d=0; d<3; ++d)
		{
			bb.edges[d].first = positions(p,d) - maxdist*radii(p);
			bb.edges[d].second = positions(p,d) + maxdist*radii(p);
		}
        tree.Insert(p, bb);
    }
    //look for nearby particles
    #pragma omp parallel for
    for(int p=0; p<Npositions[0]; ++p)
    {
        double rsq = 9.0*radii(p)*radii(p);
        std::list<int> overlapping;
        typename RTree::BoundingBox bb;
        for(int d=0; d<3; ++d)
		{
			bb.edges[d].first = positions(p,d) - maxdist*radii(p);
			bb.edges[d].second = positions(p,d) + maxdist*radii(p);
		}
		tree.Query(typename RTree::AcceptOverlapping(bb), Gatherer(overlapping));
        overlapping.sort();
        overlapping.unique();
        for(std::list<int>::const_iterator it=std::lower_bound(
            overlapping.begin(), overlapping.end(), p+1
            ); it!=overlapping.end(); ++it)
        {
            const int q = *it;
            double dsq = 0;
            for(int d=0; d<3; ++d)
                dsq += pow(positions(p,d)-positions(q,d), 2);
            if(dsq < pow(maxdist*(radii(p)+radii(q)) ,2))
            #pragma omp critical
            {
                 pairs.append(p);
                 pairs.append(q);
                 dists.append(sqrt(dsq));
            }
        }
    }
    """
    weave.inline(
        code,['positions', 'radii', 'maxdist', 'pairs', 'dists'],
        type_converters =converters.blitz,
        support_code = support,
        include_dirs = ['/home/mathieu/src/colloids/multiscale/RStarTree'],
        headers = ['"RStarTree.h"','<deque>', '<list>'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return np.resize(pairs, [len(pairs)/2,2]), np.asarray(dists)
