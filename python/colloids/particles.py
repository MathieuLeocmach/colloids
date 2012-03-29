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
import subprocess, shlex, StringIO, os, os.path
from scipy import weave
from scipy.weave import converters
import itertools

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

rstartree_path = os.path.join(
    os.path.dirname(__file__),
    '../../multiscale/RStarTree'
    )
     
support_Rtree = """
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

def weave_non_overlapping(positions, radii):
    """Give the mask of non-overlapping particles. Early bird."""
    assert len(positions)==len(radii)
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
        support_code = support_Rtree,
        include_dirs = [rstartree_path],
        headers = ['"RStarTree.h"'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return good

def get_bonds(positions, radii, maxdist=3.0):
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
        support_code = support_Rtree,
        include_dirs = [rstartree_path],
        headers = ['"RStarTree.h"','<deque>', '<list>'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return np.resize(pairs, [len(pairs)/2,2]), np.asarray(dists)
    
def get_rdf(pos, inside, Nbins=250, maxdist=30.0):
    g = np.zeros(Nbins, int)
    code = """
    //spatial indexing
    RTree tree;
    for(int p=0; p<Npos[0]; ++p)
    {
        typename RTree::BoundingBox bb;
		for(int d=0; d<3; ++d)
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
        for(int d=0; d<3; ++d)
		{
			bb.edges[d].first = pos(i,d) - maxdist;
			bb.edges[d].second = pos(i,d) + maxdist;
		}
		tree.Query(typename RTree::AcceptOverlapping(bb), Gatherer(overlapping));
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
            const int r = sqrt(disq*imaxsq)*Ng[0];
            #pragma omp atomic
            ++g(r);
        }
    }
    """
    weave.inline(
        code,['pos', 'inside', 'maxdist', 'g'],
        type_converters =converters.blitz,
        support_code = support_Rtree,
        include_dirs = [rstartree_path],
        headers = ['"RStarTree.h"','<deque>', '<list>'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return g

def get_Sq(pos, inside, qmax=10.0, maxdist=30.0, rate=1):
    assert rate>0
    assert inside.ndim == 1
    assert len(inside) == len(pos)
    qns = np.arange(0, qmax, np.pi/maxdist/rate)[rate:]
    qphis = np.linspace(0, np.pi, 4, False)
    qths = np.linspace(0, 2*np.pi, 8, False)
    qas = np.column_stack((
        np.outer(np.sin(qphis), np.cos(qths)).ravel(),
        np.outer(np.sin(qphis), np.sin(qths)).ravel(),
        np.repeat(np.cos(qphis), len(qths))
        ))
    #qas = np.eye(3)
    S = np.zeros(len(qns), float)
    code = """
    //spatial indexing
    RTree tree;
    for(int p=0; p<Npos[0]; ++p)
    {
        typename RTree::BoundingBox bb;
		for(int d=0; d<3; ++d)
		{
			bb.edges[d].first = pos(p,d) - maxdist;
			bb.edges[d].second = pos(p,d) + maxdist;
		}
        tree.Insert(p, bb);
    }
    const double imaxsq = 1.0 / pow(maxdist, 2);
    blitz::Array<std::complex<double>, 1> A(NS[0]);
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<Npos[0]; ++i)
    {
        if(!inside(i))
            continue;
        blitz::Array<std::complex<double>, 2> a(NS[0], Nqas[0]);
        a=0;
        std::list<int> overlapping;
        typename RTree::BoundingBox bb;
        for(int d=0; d<3; ++d)
		{
			bb.edges[d].first = pos(i,d) - maxdist;
			bb.edges[d].second = pos(i,d) + maxdist;
		}
		tree.Query(typename RTree::AcceptOverlapping(bb), Gatherer(overlapping));
        overlapping.sort();
        overlapping.unique();
        for(std::list<int>::const_iterator it=overlapping.begin(); it!=overlapping.end(); ++it)
        {
            const int j = *it;
            if(i==j)
                continue;
            blitz::Array<double, 1> diff(pos(j,blitz::Range::all())-pos(i,blitz::Range::all()));
            for(int l=0; l<Nqas[0]; ++l)
            {
                const double th = blitz::sum(qas(l, blitz::Range::all())*diff);
                a(blitz::Range::all(), l) += blitz::polar(1.0, qns*th);
            }
        }
        #pragma omp critical
        {
            blitz::secondIndex l;
            S += blitz::sum(blitz::norm(a), l)/(overlapping.size()-1);
            A += blitz::sum(a, l)/(double)(overlapping.size()-1);
        }
    }
    S = (S-blitz::norm(A))/Nqas[0];
    //S /= Nqas[0];
    """
    weave.inline(
        code,['pos', 'inside', 'maxdist', 'qns', 'qas', 'S'],
        type_converters =converters.blitz,
        support_code = support_Rtree,
        include_dirs = [rstartree_path],
        headers = ['"RStarTree.h"','<deque>', '<list>'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return 1+S/np.sum(inside), qns
    
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
    RTree tree;
    for(int p=0; p<Npos[0]; ++p)
    {
        typename RTree::BoundingBox bb;
		for(int d=0; d<3; ++d)
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
        for(int d=0; d<3; ++d)
		{
			bb.edges[d].first = pos(i,d) - maxdist*radii(i);
			bb.edges[d].second = pos(i,d) + maxdist*radii(i);
		}
		tree.Query(typename RTree::AcceptOverlapping(bb), Gatherer(overlapping));
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
    """
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
    RTree tree;
    for(int p=0; p<Npos[0]; ++p)
    {
        typename RTree::BoundingBox bb;
		for(int d=0; d<3; ++d)
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
        for(int d=0; d<3; ++d)
		{
			bb.edges[d].first = pos(i,d) - maxdist;
			bb.edges[d].second = pos(i,d) + maxdist;
		}
		tree.Query(typename RTree::AcceptOverlapping(bb), Gatherer(overlapping));
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
    """
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
    g = np.zeros(Nbins, int)
    code = """
    //spatial indexing
    RTree tree;
    for(int p=0; p<Npos[0]; ++p)
    {
        typename RTree::BoundingBox bb;
		for(int d=0; d<3; ++d)
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
		tree.Query(typename RTree::AcceptOverlapping(bb), Gatherer(overlapping));
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
    """
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
            for(int dim=0; dim<3;++dim)
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
    code = """
    //spatial indexing
    RTree tree;
    for(int p=0; p<Npos0[0]; ++p)
    {
        typename RTree::BoundingBox bb;
        for(int d=0; d<3; ++d)
        {
            bb.edges[d].first = pos0(p,d) - maxdist*radii0(p);
            bb.edges[d].second = pos0(p,d) + maxdist*radii0(p);
        }
        tree.Insert(p, bb);
    }
    //look for nearby particles
    #pragma omp parallel for
    for(int p=0; p<Npos1[0]; ++p)
    {
        std::list<int> overlapping;
        typename RTree::BoundingBox bb;
        for(int d=0; d<3; ++d)
        {
            bb.edges[d].first = pos1(p,d) - maxdist*radii1(p);
            bb.edges[d].second = pos1(p,d) + maxdist*radii1(p);
        }
        tree.Query(typename RTree::AcceptOverlapping(bb), Gatherer(overlapping));
        overlapping.sort();
        overlapping.unique();
        for(std::list<int>::const_iterator it=overlapping.begin(); it!=overlapping.end(); ++it)
        {
            const int q = *it;
            double dsq = 0;
            for(int d=0; d<3; ++d)
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
    """
    weave.inline(
        code,['pos0', 'radii0', 'pos1', 'radii1', 'maxdist', 'pairs', 'dists'],
        type_converters =converters.blitz,
        support_code = support_Rtree,
        include_dirs = [rstartree_path],
        headers = ['"RStarTree.h"','<deque>', '<list>'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return np.resize(pairs, [len(pairs)/2,2]), np.asarray(dists)
    
def get_links_size(pos0, radii0, pos1, radii1, maxdist=1.0):
    pairs, distances = particles.get_links(pos0, radii0, pos1, radii1, maxdist)
    code = """
    #pragma omp parallel for
    for(int i=0; i<Npairs[0]; ++i)
    {
        const int p = pairs(i,0), q = pairs(i,1);
        distances(i) += pow(radii0(i) - radii1(j), 2);
    }
    """
    weave.inline(
        code, ['pairs', 'radii0', 'radii1', 'distances'],
        type_converters =converters.blitz,
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return pairs, dists  

class Linker:
    def __init__(self, nb_initial_pos):
        self.tr2pos = [[p] for p in range(nb_initial_pos)]
        self.pos2tr = [np.arange(nb_initial_pos)]
        self.trajstart = [0 for p in range(nb_initial_pos)]
        
    def addFrame(self, frame_size, pairs, distances):
        assert len(pairs) == len(distances)
        assert pairs[:,1].max() < frame_size, "The largest particle index in the new frame is larger than the new frame size"
        #sort the possible links by increasing distances
        pairs = pairs[np.argsort(distances)]
        #any position can be linked only once. At init none are linked
        from_used = np.zeros(len(self.pos2tr[-1]), bool)
        to_used = np.zeros(frame_size, bool)
        #create the new frame
        newframe = np.zeros(frame_size, int)
        #link the bounded positions into trajectories
        for p,q in pairs:
            if(from_used[p] or to_used[q]):continue
            from_used[p] = True
            to_used[q] = True
            tr = self.pos2tr[-1][p]
            newframe[q] = tr
            self.tr2pos[tr].append(q);
        #the trajectories of the previous frame that are not linked in the new frame are terminated by construction
        #but the trajectories starting in the new frame have to be created
        notused = np.where(np.bitwise_not(to_used))[0]
        newframe[notused] = np.arange(len(notused)) + len(self.tr2pos)
        self.tr2pos += [[p] for p in notused]
        self.trajstart += [len(self.pos2tr) for p in notused]
        #add the new frame
        self.pos2tr.append(newframe)

