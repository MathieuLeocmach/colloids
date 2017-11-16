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
#import rtree.index
import numexpr
import subprocess, shlex, os, os.path
from scipy.spatial import cKDTree as KDTree
#from scipy import weave
#from scipy.weave import converters
import itertools
import os.path
from colloids import periodic

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
        raise DeprecationWarning('This class is an ugly bit of code not in use since 2010')

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
        
        
class Grid:
    """Spatial indexing structure that stores the indices of objects"""
    def __init__(self, dimension=3, sizes=[256,256,256], divisions=[16,16,16]):
        try:
            assert len(sizes) == dimension
        except TypeError:
            sizes = [float(sizes)]*dimension
        try:
            assert len(divisions) == dimension
        except TypeError:
            divisions = [float(divisions)]*dimension
        self.__sizes = np.array(sizes)
        self.__divisions = np.array(divisions, int)
        self.__steps = np.cumprod(self.__divisions)
        self.__steps[1:] = self.__steps[:-1]
        self.__steps[0] = 1
        self.__data = [[] for i in range(np.prod(divisions))]
        self.__index2cell = []
        
    def __len__(self):
        return len(self.__index2cell)
        
    def __get_cell_id(self, cell_coords):
        return np.sum(self.__steps * cell_coords)
        
    def __get_cell(self, cell_coords):
        return self.__data[self.__get_cell_id(cell_coords)]
        
    def add_point(self, point):
        coords = np.array(point*self.__divisions/self.__sizes, int)
        cell = self.__get_cell_id(coords)
        self.__data[cell].append(len(self))
        self.__index2cell.append(cell)
    
    def get_near(self, point, d=1):
        coords = np.array(point*self.__divisions/self.__sizes, int)
        return np.unique1d([
            i for co in itertools.product(*[
                range(max(0, c-d), min(div, c+d+1))
                for c, div in zip(coords, self.__divisions)
                ]) for i in self.__get_cell(co)
            ])
            
    def get_inside(self, d=1):
        return [
            i for co in itertools.product(*[
                range(d, div-d)
                for div in self.__divisions
                ]) for i in self.__get_cell(co)
            ]
    def iter_pairs(self, d=1):
        for cp in itertools.product(*[range(d, div-d) for div in self.__divisions]):
            for cq in itertools.product(*[
                range(c-d, c+d+1) 
                for c, div in zip(cp, self.__divisions)
                ]):
                for p,q in itertools.product(self.__get_cell(cp), self.__get_cell(cq)):
                    if(p!=q):
                        yield p,q
        

def non_overlapping(positions, radii):
    """Give the mask of non-overlapping particles. Early bird."""
    assert len(positions)==len(radii)
    tree = KDTree(positions)
    rmax = radii.max()
    good = np.ones(len(positions), dtype=bool)
    for i, (p,r) in enumerate(zip(positions, radii)):
        if not good[i]:
            continue
        for j in tree.query_ball_point(p, rmax + r):
            if j==i or not good[j]:
                continue
            #the python loop is actually faster than numpy on small(3) arrays
            s = 0.0
            for pd, qd in zip(positions[j], p):
                s = (pd - pd)**2
            if s < (radii[j] + r)**2:
                good[j] = False
    return good
    
def non_halfoverlapping(positions, radii):
    """Give the mask of non-half overlapping particles. Early bird.
    Half overlap is when $r_{ij} < \max(R_i, Rj)$"""
    assert len(positions)==len(radii)
    tree = KDTree(positions)
    rmax = radii.max()
    good = np.ones(len(positions), dtype=bool)
    for i, (p,r) in enumerate(zip(positions, radii)):
        if not good[i]:
            continue
        for j in tree.query_ball_point(p, rmax + r):
            if j==i or not good[j]:
                continue
            #the python loop is actually faster than numpy on small(3) arrays
            s = 0.0
            for pd, qd in zip(positions[j], p):
                s = (pd - pd)**2
            if s < max((radii[j], r))**2:
                good[j] = False
    return good
    
def get_bonds(positions, radii, maxdist=3.0):
    """Bonds by relative distances, such that $r_{ij} < maxdist (R_i + R_j)$.
    
    Returns pairs, distances. Pairs are sorted and unique."""
    assert len(positions)==len(radii)
    tree = KDTree(positions)
    rmax = radii.max()
    #fetch all potential pairs, already sorted
    pairs = np.array([
        [i,j] 
        for i, js in enumerate(tree.query_ball_tree(tree, 2*rmax*maxdist))
        for j in sorted(js)
        if i<j
        ])
    if len(pairs) == 0:
        return np.zeros((0,2), int), np.zeros(0)
    #compute all pair's square distances via numpy
    dists = np.sum((positions[pairs[:,0]] - positions[pairs[:,1]])**2, -1)
    #filter out the pairs that are too far
    good = dists < maxdist**2 * radii[pairs].sum(-1)**2
    return pairs[good], np.sqrt(dists[good])
    
def get_N_ngbs(positions, radii, N=12, maxdist=3.0, edge = None):
    """N first neighbours, with a maximum relative distances, such that $r_{ij} < maxdist (R_i + R_j)$.
    If a potential neighbour is further away than the distance to the edge of the field of view, 
    the current particle of interest is considered as "on the edge" and the neighbour not taken into account.
    
    Returns neighbours, inside"""
    assert len(positions)==len(radii)
    if edge is None:
        edge = (positions.min(0), positions.max(0))
    #initialize the geometry of each particle
    to_edge = np.minimum((positions - edge[0]).min(-1), (edge[0] - positions).min(-1))**2
    inside = np.full(len(positions), True, dtype=bool)
    neighbours = np.full([len(positions), N], -1, dtype=int)
    tree = KDTree(positions)
    rmax = radii.max()
    for i, js in enumerate(tree.query_ball_tree(tree, 2*rmax*maxdist)):
        disq = np.sum((positions[js] - positions[i])**2, -1)
        ags = np.argsort(disq)[:N]
        if disq[ags[-1]] < to_edge[i]:
            neighbours[i, :len(js)] = np.array(js)[ags]
        else:
            inside[i] = False
            N2 = np.where(disq[ags] < to_edge[i])[0][0]+1
            neighbours[i, :N2] = np.array(js)[ags[:N2]]
    return neighbours, inside
    
def ngbN2bonds(ngbs):
    """Convert the result of get_N_ngbs into a sorted array of bonds"""
    bonds = []
    for i,l in enumerate(ngbs):
        for n in l:
            if n<0: continue
            bonds.append(tuple(sorted([i,n])))
    return np.array([[a,b] for a,b in set(bonds)])
    
    
def bonds2ngbs_list(bonds, N):
    """Returns a list of arrays of neighbours from bond data. N is the number of particles"""
    ngbs = [[] for i in range(N)]
    for a,b in bonds:
        ngbs[a].append(b)
        ngbs[b].append(a)
    return map(np.array, ngbs)
    
def bonds2ngbs(bonds, N):
    """Returns an array of neighbours from bond data. N is the number of particles"""
    ngbs = -np.ones([N, np.histogram(bonds, bins=np.arange(N+1))[0].max()], int)
    if bonds.shape[-1]>0:
        for a,b in bonds:
            ngbs[a, np.where(ngbs[a]==-1)[0][0]] = b
            ngbs[b, np.where(ngbs[b]==-1)[0][0]] = a
    return ngbs

rstartree_path = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    '../../multiscale/RStarTree'
    ))
     
support_Rtree = """
template <typename Leaf>
struct Gatherer {
    std::list<int> *gathered;
    bool ContinueVisiting;

    Gatherer(std::list<int> &result) : gathered(&result), ContinueVisiting(true) {};

    void operator()(const Leaf * const leaf)
    {
	    gathered->push_back(leaf->leaf);
    }
};
"""

def get_rdf(pos, inside, Nbins=250, maxdist=30.0):
    """Radial distribution function, not normalised.
    For each particle tagged as inside, count the particles around and bin then with respect to distance. Need to be normalised by inside.sum() and density x volume of the spherical shell between r and r+maxdist/Nbins.
    
     - pos is a Nxd array of coordinates, with d the dimension of space
     - inside is a N array of booleans. For example all particles further away than maxdist from any edge of the box.
     - Nbins is the number of bins along r
     - maxdist is the maximum distance considered"""
    g = np.zeros(Nbins, int)
    #conversion factor between indices and bins
    l2r = (Nbins-1)/maxdist
    #spatial indexing
    tree = KDTree(pos, 12)
    centertree = KDTree(pos[inside], 12)
    centerindex = np.where(inside)[0]
    #all pairs of points closer than maxdist with their distances in a record array
    query = centertree.sparse_distance_matrix(tree, maxdist, output_type='ndarray')
    #keep only pairs where the points are distinct
    query['i'] = centerindex[query['i']]
    good = query['i'] != query['j']
    query = query[good]
    #binning
    rs = (query['v'] * l2r).astype(int)
    np.add.at(g, rs, 1)
    return g

def structure_factor(positions, Nbins, Ls=[203.0]*3, maxNvec=30, field=None):
    """Compute the structure factor of the positions in a non-periodic rectinilear box of dimensions Ls. 
    To avoid the effects of box window, a Hanning window is applied. The three first coefficients are probably affected by this windowing. Be careful of the boundaries: do not leave empty margins around the data."""
    if field is None:
        field = np.ones(len(positions))
    assert len(field) == len(positions)
    #use the periodic version of get_Sq, but first the field is
    # windowed by a Hanning
    window = np.prod(0.5 * (1 + np.sin(2*np.pi*positions/Ls)), -1)
    return periodic.rectangular_Sq(positions, Nbins, Ls, maxNvec, field*window) / (window**2).mean()
    
def S2d_3d(S, q, width, dq=None):
    """Approximate restoration of a 3D structure factor from a structure factor taken in 2D on a slice of finite width. See Hiroo Totsuji and Chieko Totsuji, Physical Review E 85, 031139 (2012)."""
    if dq is None:
        if q[0]==0:
            dq = q[1]
        else: 
            dq = q[0]
    return 1 + np.asarray([
        np.sum((
            -2*np.gradient(S)[1+i:]/width + 
            dq*q[1+i:]/6*(width - width**3*(q[1+i:]**2-k**2)/60)*(S[1+i:]-1)
            )/np.sqrt(q[i+1:]**2-k**2))
        for i,k in enumerate(q)
        ])
    
def get_srdf(pos, radii, inside, Nbins=250, maxdist=3.0):
    g = np.zeros(Nbins, int)
    code = """
    //spatial indexing
    typedef RStarTree<int, %(dim)d, 4, 32, double> RTree;
    RTree tree;
    for(int p=0; p<Npos[0]; ++p)
    {
        typename RTree::BoundingBox bb;
		for(int d=0; d<Npos[1]; ++d)
		{
			bb.edges[d].first = pos(p,d) - maxdist*radii(p);
			bb.edges[d].second = pos(p,d) + maxdist*radii(p);
		}
        tree.Insert(p, bb);
    }
    const double imaxsq = 1.0 / pow(maxdist, 2);
    #pragma omp parallel for
    for(int i=0; i<Npos[0]; ++i)
    {
        if(!inside(i))
            continue;
        std::list<int> overlapping;
        typename RTree::BoundingBox bb;
        for(int d=0; d<Npos[1]; ++d)
		{
			bb.edges[d].first = pos(i,d) - maxdist*radii(i);
			bb.edges[d].second = pos(i,d) + maxdist*radii(i);
		}
		tree.Query(typename RTree::AcceptOverlapping(bb), Gatherer<RTree::Leaf>(overlapping));
        overlapping.sort();
        overlapping.unique();
        for(std::list<int>::const_iterator it=overlapping.begin(); it!=overlapping.end(); ++it)
        {
            const int j = *it;
            if(i==j)
                continue;
            const double disq = blitz::sum(blitz::pow(pos(j,blitz::Range::all())-pos(i,blitz::Range::all()), 2)) / pow(radii(i)+radii(j), 2);
            if(disq*imaxsq>=1.0)
                continue;
            const int r = sqrt(disq*imaxsq)*Ng[0];
            #pragma omp atomic
            ++g(r);
        }
    }
    """%{'dim':pos.shape[1]}
    weave.inline(
        code,['pos', 'radii', 'inside', 'maxdist', 'g'],
        type_converters =converters.blitz,
        support_code = support_Rtree,
        include_dirs = [rstartree_path],
        headers = ['"RStarTree.h"','<deque>', '<list>'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return g


def get_planar_rdf(pos, inside, Nbins=250, maxdist=30.0, maxangle=np.pi/3):
    """Construct radial distribution function considering only the bonds in the XY plane"""
    g = np.zeros(Nbins, int)
    maxsin = float(np.sin(maxangle)**2)
    code = """
    //spatial indexing
    typedef RStarTree<int, %(dim)d, 4, 32, double> RTree;
    RTree tree;
    for(int p=0; p<Npos[0]; ++p)
    {
        typename RTree::BoundingBox bb;
		for(int d=0; d<Npos[1]; ++d)
		{
			bb.edges[d].first = pos(p,d) - maxdist;
			bb.edges[d].second = pos(p,d) + maxdist;
		}
        tree.Insert(p, bb);
    }
    const double imaxsq = 1.0 / pow(maxdist, 2);
    #pragma omp parallel for
    for(int i=0; i<Npos[0]; ++i)
    {
        if(!inside(i))
            continue;
        std::list<int> overlapping;
        typename RTree::BoundingBox bb;
        for(int d=0; d<Npos[1]; ++d)
		{
			bb.edges[d].first = pos(i,d) - maxdist;
			bb.edges[d].second = pos(i,d) + maxdist;
		}
		tree.Query(typename RTree::AcceptOverlapping(bb), Gatherer<RTree::Leaf>(overlapping));
        overlapping.sort();
        overlapping.unique();
        for(std::list<int>::const_iterator it=overlapping.begin(); it!=overlapping.end(); ++it)
        {
            const int j = *it;
            if(i==j)
                continue;
            const double disq = blitz::sum(blitz::pow(pos(j,blitz::Range::all())-pos(i,blitz::Range::all()), 2));
            if(disq*imaxsq>=1.0)
                continue;
			if(pow((double)(pos(j,2)-pos(i,2)), 2) > maxsin*disq)
				continue;
            const int r = sqrt(disq*imaxsq)*Ng[0];
            #pragma omp atomic
            ++g(r);
        }
    }
    """%{'dim':pos.shape[1]}
    weave.inline(
        code,['pos', 'inside', 'maxdist', 'maxsin','g'],
        type_converters =converters.blitz,
        support_code = support_Rtree,
        include_dirs = [rstartree_path],
        headers = ['"RStarTree.h"','<deque>', '<list>'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return g
    
def get_slice_rdf(pos, inside, Nbins=250, maxdist=30.0, width=10.0):
    """Construct radial distribution function in a thin slice around each central particle"""
    assert pos.shape[1]==3
    g = np.zeros(Nbins, int)
    code = """
    //spatial indexing
    typedef RStarTree<int, %(dim)d, 4, 32, double> RTree;
    RTree tree;
    for(int p=0; p<Npos[0]; ++p)
    {
        typename RTree::BoundingBox bb;
		for(int d=0; d<Npos[1]; ++d)
		{
			bb.edges[d].first = pos(p,d) - maxdist;
			bb.edges[d].second = pos(p,d) + maxdist;
		}
        tree.Insert(p, bb);
    }
    const double imaxsq = 1.0 / pow(maxdist, 2);
    #pragma omp parallel for
    for(int i=0; i<Npos[0]; ++i)
    {
        if(!inside(i))
            continue;
        std::list<int> overlapping;
        typename RTree::BoundingBox bb;
        for(int d=0; d<2; ++d)
		{
			bb.edges[d].first = pos(i,d) - maxdist;
			bb.edges[d].second = pos(i,d) + maxdist;
		}
		bb.edges[2].first = pos(i,2) - 0.5*width;
		bb.edges[2].second = pos(i,2) + 0.5*width;
		tree.Query(typename RTree::AcceptOverlapping(bb), Gatherer<RTree::Leaf>(overlapping));
        overlapping.sort();
        overlapping.unique();
        for(std::list<int>::const_iterator it=overlapping.begin(); it!=overlapping.end(); ++it)
        {
            const int j = *it;
            if(i==j)
                continue;
            if(abs((double)(pos(j,2)-pos(i,2))) > 0.5*width)
				continue;
            const double disq = blitz::sum(blitz::pow(pos(j,blitz::Range(0,1))-pos(i,blitz::Range(0,1)), 2));
            if(disq*imaxsq>=1.0)
                continue;
            const int r = sqrt(disq*imaxsq)*Ng[0];
            #pragma omp atomic
            ++g(r);
        }
    }
    """%{'dim':pos.shape[1]}
    weave.inline(
        code,['pos', 'inside', 'maxdist', 'width','g'],
        type_converters =converters.blitz,
        support_code = support_Rtree,
        include_dirs = [rstartree_path],
        headers = ['"RStarTree.h"','<deque>', '<list>'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return g

def get_rdf_window(pos, box, Nbins=250, maxdist=30.0):
    assert len(box)==2
    assert len(box[0]) == 3
    L = (box[1]-box[0])
    g = np.zeros(Nbins)
    #compute the weight of each particle according to the Hamming window function
    #weights = np.prod(0.54 - 0.46*np.sin(2*np.pi*(pos-box[0])/(box[1]-box[0])), -1)
    #compute the weight of each particle according to the Hann window function
    weights = np.prod(0.5*(1-np.sin(2*np.pi*(pos-box[0])/(box[1]-box[0]))), -1)
    code = """
    //spatial indexing
    const double imaxsq = 1.0 / pow(maxdist, 2);
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<Npos[0]-1; ++i)
    {
        for(int j=i+1; j<Npos[0]; ++j)
        {
            const double disq = blitz::sum(blitz::pow(blitz::fmod(pos(j,blitz::Range::all())-pos(i,blitz::Range::all())+1.5*L, L)-0.5*L, 2));
            if(disq*imaxsq>=1.0)
                continue;
            const int r = sqrt(disq*imaxsq)*Ng[0];
            #pragma omp atomic
            g(r) += weights(i);
        }
    }
    """
    weave.inline(
        code,['pos', 'weights', 'L','maxdist', 'g'],
        type_converters =converters.blitz,
        #support_code = support_Rtree,
        include_dirs = [rstartree_path],
        headers = ['"RStarTree.h"','<deque>', '<list>'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return g
    
def get_s2(rdfbin, nd, av=5000, binlength=0.02):
    assert(rdfbin.ndim==2)
    s2 = np.zeros(len(rdfbin))
    factors = 1.0/(av*nd*np.diff(
        4*np.pi/3.*(np.arange(rdfbin.shape[1]+1)*binlength)**3
        ))
    code = """
    #pragma omp parallel for
    for(int i=0; i<Nrdfbin[0]; ++i)
    {
        for(int j=0; j<Nrdfbin[1]; ++j)
        {
            if(rdfbin(i,j)<1)
                s2(i) += 1;
            else
            {
                const double g = rdfbin(i,j) * factors(j);
                s2(i) += g * log(g) - g + 1;
            }
        }
        s2(i) *= -nd/2; 
    }
    """
    weave.inline(
        code,['rdfbin', 'factors', 'nd', 's2'],
        type_converters =converters.blitz,
        extra_compile_args =['-O3 -fopenmp -mtune=native'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return s2
    
def spatial_correlation(pos, is_center, values, Nbins, maxdist):
    """
    Spatial correlation of the qlms and the Qlms
    """
    assert len(pos) == len(values)
    assert len(is_center) == len(pos)
    maxsq = float(maxdist**2)
    h = np.zeros(Nbins)
    g = np.zeros(Nbins, int)
    code = """
    #pragma omp parallel for
    for(int i=0; i<Npos[0]; ++i)
    {
        if(!is_center(i)) 
            continue;
        for(int j=0; j<Npos[0]; ++j)
        {
            if(i==j) continue;
            double disq = 0.0;
            for(int dim=0; dim<Npos[1];++dim)
                disq += pow(pos(i,dim)-pos(j,dim), 2);
            if(disq>=(double)maxsq)
                continue;
            const int r = sqrt(disq/(double)maxsq)*Nbins;
            const double prod = values(i)*values(j);
            #pragma omp critical
            {
                ++g(r);
                h(r) += prod;
            }
        }
    }
    """
    weave.inline(
        code,['values', 'pos', 'maxsq', 'Nbins', 'h', 'g','is_center'],
        type_converters =converters.blitz,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return h, g

def get_links(pos0, radii0, pos1, radii1, maxdist=1.0):
    """Get the pairs of particles closer than maxdist in two configurations and their distances"""
    pairs = []
    dists = []
    maxdist_=maxdist
    code = """
    const double maxdist=maxdist_;
    //spatial indexing
    typedef RStarTree<int, %(dim)d, 4, 32, double> RTree;
    RTree tree;
    for(int p=0; p<Npos0[0]; ++p)
    {
        typename RTree::BoundingBox bb;
        for(int d=0; d<Npos1[1]; ++d)
        {
            const double r = radii0(p);
            bb.edges[d].first = pos0(p,d) - maxdist*r;
            bb.edges[d].second = pos0(p,d) + maxdist*r;
        }
        tree.Insert(p, bb);
    }
    //look for nearby particles
    #pragma omp parallel for
    for(int p=0; p<Npos1[0]; ++p)
    {
        std::list<int> overlapping;
        typename RTree::BoundingBox bb;
        for(int d=0; d<Npos1[1]; ++d)
        {
            const double r = radii1(p);
            bb.edges[d].first = pos1(p,d) - maxdist*r;
            bb.edges[d].second = pos1(p,d) + maxdist*r;
        }
        tree.Query(typename RTree::AcceptOverlapping(bb), Gatherer<RTree::Leaf>(overlapping));
        overlapping.sort();
        overlapping.unique();
        for(std::list<int>::const_iterator it=overlapping.begin(); it!=overlapping.end(); ++it)
        {
            const int q = *it;
            double dsq = 0;
            for(int d=0; d<Npos1[1]; ++d)
                dsq += pow(pos1(p,d)-pos0(q,d), 2);
            if(dsq < pow(maxdist*(radii1(p)+radii0(q)), 2))
            #pragma omp critical
            {
                 pairs.append(q);
                 pairs.append(p);
                 dists.append(sqrt(dsq));
            }
        }
    }
    """%{'dim':pos1.shape[1]}
    weave.inline(
        code,['pos0', 'radii0', 'pos1', 'radii1', 'maxdist_', 'pairs', 'dists'],
        type_converters =converters.blitz,
        support_code = support_Rtree,
        include_dirs = [rstartree_path],
        headers = ['"RStarTree.h"','<deque>', '<list>'],
        extra_compile_args =['-O3 -fopenmp -mtune=native'],
        extra_link_args=['-lgomp'],
        verbose=2)
    return np.resize(pairs, [len(pairs)/2,2]), np.asarray(dists)
    
def get_links_size(pos0, radii0, pos1, radii1, maxdist=1.0):
    pairs, distances = get_links(pos0, radii0, pos1, radii1, maxdist)
    code = """
    #pragma omp parallel for
    for(int i=0; i<Npairs[0]; ++i)
    {
        const int p = pairs(i,0), q = pairs(i,1);
        distances(i) += pow(radii0(p) - radii1(q), 2);
    }
    """
    weave.inline(
        code, ['pairs', 'radii0', 'radii1', 'distances'],
        type_converters =converters.blitz,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return pairs, distances  

class Linker:
    def __init__(self, nb_initial_pos):
        self.tr2pos = [[p] for p in range(nb_initial_pos)]
        self.pos2tr = [np.arange(nb_initial_pos)]
        self.trajstart = [0 for p in range(nb_initial_pos)]
        self.nbtrajs = [nb_initial_pos]
        
    def loadFrame(self, frame):
        """Load a precalculated pos2tr frame."""
        t = len(self.pos2tr)
        self.pos2tr.append(frame)
        for p,tr in enumerate(frame):
            if tr<self.nbtrajs[-1]:
                self.tr2pos[tr].append(p)
            else:
                self.trajstart.append(t)
                self.tr2pos.append([p])
        self.nbtrajs.append(len(self.trajstart))
        
    def addFrame(self, frame_size, pairs, distances):
        assert len(pairs) == len(distances)
        if len(pairs)==0:
            self.tr2pos += [[p] for p in range(frame_size)]
            self.trajstart += [len(self.pos2tr) for p in range(frame_size)]
            self.pos2tr.append(np.arange(frame_size)+self.nbtrajs[-1])
            self.nbtrajs.append(len(self.trajstart))
            return
        assert pairs[:,1].max() < frame_size, "The largest particle index in the new frame is larger than the new frame size"
        #sort the possible links by increasing distances
        pairs = pairs[np.argsort(distances)]
        #any position can be linked only once. At init none are linked
        from_used = np.zeros(len(self.pos2tr[-1]), bool)
        to_used = np.zeros(frame_size, bool)
        #create the new frame
        newframe = np.zeros(frame_size, int)
        code = """
        for(int i=0; i<Npairs[0]; ++i)
        {
            const int p = pairs(i,0), q = pairs(i,1);
            if(from_used(p) | to_used(q)) continue;
            from_used(p) = true;
            to_used(q) = true;
            const int tr = pos2tr(p);
            newframe(q) = tr;
            py::object o(tr2pos[tr]);
            ((py::list *)&o)->append(q);
        }
        """
        pos2tr = self.pos2tr[-1]
        tr2pos = self.tr2pos
        weave.inline(
            code, ['pairs', 'from_used', 'to_used', 'pos2tr', 'tr2pos', 'newframe'],
            type_converters =converters.blitz,
            extra_compile_args =['-O3 -fopenmp -mtune=native'],
            extra_link_args=['-lgomp'],
            verbose=2)
        #link the bounded positions into trajectories
        #for p,q in pairs:
         #   if(from_used[p] or to_used[q]):continue
          #  from_used[p] = True
           # to_used[q] = True
            #tr = self.pos2tr[-1][p]
            #newframe[q] = tr
            #self.tr2pos[tr].append(q);
        #the trajectories of the previous frame that are not linked in the new frame are terminated by construction
        #but the trajectories starting in the new frame have to be created
        notused = np.where(np.bitwise_not(to_used))[0]
        newframe[notused] = np.arange(len(notused)) + len(self.tr2pos)
        self.tr2pos += [[p] for p in notused]
        self.trajstart += [len(self.pos2tr) for p in notused]
        self.nbtrajs.append(len(self.trajstart))
        #add the new frame
        self.pos2tr.append(newframe)
        
    def update(self, pairs, distances):
        """Continue trajectories terminated two steps ago (not existing one step ago). The first column of pairs contains trajectory indices, the second column contains present (last frame) positions indices."""
        assert len(self.pos2tr)>2
        assert len(pairs) == len(distances)
        #sort the possible links by increasing distances
        pairs = pairs[np.argsort(distances)]
        #any position can be linked only once. Initialize the ones allready linked.
        from_used = np.zeros(self.nbtrajs[-2], bool)
        from_used[self.pos2tr[-2]] = True
        to_used = self.pos2tr[-1] < self.nbtrajs[-2]
        #filter the links
        links = []
        for tr, p in pairs:
            if from_used[tr] or to_used[p]: continue
            from_used[tr] = True
            to_used[p] = True
            links.append((tr, p))
        links = np.array(links, dtype=int)
        links = links[np.argsort(links[:,0])]
        #create intermediate position indices
        self.pos2tr[-2] = np.concatenate((self.pos2tr[-2], links[:,0]))
        #set grown trajectories in the present frame
        self.pos2tr[-1][links[:,1]] = links[:,0]
        #remove the trajectories previously starting in the present frame
        del self.trajstart[self.nbtrajs[-2]:]
        del self.tr2pos[self.nbtrajs[-2]:]
        #grow trajectories
        for (tr, p), inter in zip(links, range(len(self.pos2tr[-2]) - len(links), len(self.pos2tr[-2]))):
            self.tr2pos[tr] += [inter, p]
        #recreate the trajectories now starting in the present frame
        notused = np.where(np.bitwise_not(to_used))[0]
        self.pos2tr[-1][notused] = np.arange(len(notused)) + self.nbtrajs[-2]
        self.tr2pos += [[p] for p in notused]
        self.trajstart += [len(self.pos2tr)-1 for p in notused]
        self.nbtrajs[-1] = len(self.trajstart)
        #note the trajectories that have grown in the process
        has_grown = np.zeros(len(self.trajstart), bool)
        has_grown[links[:,0]] = True
        return has_grown
        
    def save(self, f):
        """write the trajectory data to an opend file"""
        for start, tr in zip(self.trajstart, self.tr2pos):
            f.write('%d\n'%start)
            f.write('\t'.join(['%d'%p for p in tr])+'\n')
            
def link_save(path, dt, size, radius=4.32692):
    """Link and save the trajectories"""
    pattern = path+'_t%0'+('%dd.dat'%len('%d'%size))
    #last frame is often empty when acquisition was interrupted
    if len([line for line in open(pattern%(size-1))])<3 or len(np.loadtxt(pattern%(size-1), skiprows=2))==0:
        size -=1
    #linking
    pos1 = np.loadtxt(pattern%0, skiprows=2)
    linker = Linker(len(pos1))
    np.savetxt((pattern[:-3]+'p2tr')%0, linker.pos2tr[0], fmt='%d')
    for t in range(size-1):
        pos0 = pos1
        pos1 = np.loadtxt(pattern%(t+1), skiprows=2)
        pairs, distances = get_links(pos0, np.ones(len(pos0)), pos1, np.ones(len(pos1)), maxdist=5)
        linker.addFrame(len(pos1), pairs, distances)
        #saving p2tr
        np.savetxt((pattern[:-3]+'p2tr')%(t+1), linker.pos2tr[-1], fmt='%d')
    #saving in the same format as the c++
    with open(path+'.traj', 'w') as f:
        f.write('%g\t%g\n'%(radius, dt))
        f.write('%s\n_t\n0\t%d\n'%(os.path.basename(pattern%0), size))
        linker.save(f)
    
