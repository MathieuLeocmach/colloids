import numpy as np
import colloids.particles
from  colloids.particles import Particles
from colloids import boo, track
from colloids import experiment as xp
import colloids.statistics as sta
import colloids.povray as pov
import os.path, subprocess
from tables import *
import numexpr
from scipy import weave
from scipy.weave import converters
from pygraph.classes.graph import graph
from pygraph.algorithms.accessibility import connected_components
import anfft

def export_ngb_q6(x,mrads,pdata):
    for t, (pd, rm) in enumerate(zip(pdata, mrads)):
        p = Particles(np.ascontiguousarray(pd[:,:3]), rm)
        sel = [i for i in p.index.intersection(p.inside_box(8.0))]
        ngbs = -1*np.ones([len(pd), 12], int)
        for i in sel:
            ngbs[i] = p.get_N_ngbs(i, 5.0)
        np.save(x.get_format_string('ngbs','npy')%t, ngbs)
        q6ms = np.zeros([len(pd), 7], np.complex128)
        for i in sel:
            q6ms[i] = p.get_qlm(i, ngbs[i], ls=[6])[0]
        np.save(x.get_format_string('_q6m','npy')%t, q6ms)
        Q6ms = np.zeros_like(q6ms)
        for i in sel:
            Q6ms[i] = (q6ms[i] + q6ms[ngbs[i]].sum(0)) / 13.0
        np.save(x.get_format_string('_Q6m','npy')%t, Q6ms)

class Positions(IsDescription):
    idnumber  = Int64Col(pos=1)
    trnumber  = Int64Col(pos=2, dflt=-1)
    x         = Float64Col(pos=3)
    y         = Float64Col(pos=4)
    z         = Float64Col(pos=5)
    r         = Float64Col(pos=6)

class bonds12(IsDescription):
    idnumber  = Int64Col(pos=1)
    ngbs      = Int64Col(pos=2, shape=12)
    inside    = BoolCol(pos=3)

class Bond_order(IsDescription):
    idnumber  = Int64Col(pos=1)
    q6m       = ComplexCol(16, pos=2, shape=7)
    Q6m       = ComplexCol(16, pos=3, shape=7)
    q6        = Float64Col(pos=4)
    w6        = Float64Col(pos=5)
    Q6        = Float64Col(pos=6)

#create the group that will contain the whole sample data
#sample = h5file.createGroup("/", 'LS%d'%phi, 'Large Sample %d'%phi)
def fill_positions_trajs(h5file, x, group):
    xyzr = [np.load(fname, 'r') for t, fname in x.enum('_mrad', 'npy')]
    trids = [np.empty(len(p), int) for p in xyzr]
    trajs_start, trajs = xp.load_trajectories(
        xp.os.path.join(x.path, x.head+'.traj')
        )
    for trid, (t0, tr) in enumerate(zip(trajs_start, trajs)):
        for t, i in enumerate(tr):
            trids[t0+t][i] = trid
    group._v_attrs.unit_size = x.unitLength
    group._v_attrs.unit_time = x.dt
    group._v_attrs.n_time_steps = x.size
    group._v_attrs.Temperature = x.T
    meanrad = np.mean([np.mean(p[:,3]**3) for p in xyzr])**(1/3.)
    group._v_attrs.mean_diam = 2*np.sqrt(2) * meanrad
    group._v_attrsBrownian_time = xp.br(meanrad*np.sqrt(2)*x.unitLength, T)
    #trajectories are a ragged array
    trajectory_array = h5file.createVLArray(
        group, 'trajectories',
        Int64Atom(shape=()),
        "starting time and then positions id at each subsequent times")
    for t0, tr in zip(trajs_start, trajs):
        trajectory_array.append(np.array([t0]+tr))
    trajectory_array.flush()
    #data at each time
    for t, part in enumerate(xyzr):
        if hasattr(group, 't%03d'%t):
            time_step = getattr(group, 't%03d'%t)
            if hasattr(time_step, 'positions'):
                continue
        else:
            #create the group that will contain the data a the given time step
            time_step = h5file.createGroup(group, 't%03d'%t, 'All data for time %03d'%t)
            time_step._v_attrs.t = t
        #position table
        pos = h5file.createTable(
            time_step, 'positions',
            Positions, 'position and size of each particle'
        )
        pos.flush()
        prow = pos.row
        for i in range(len(part)):
            prow['idnumber'] = i
            prow['trnumber'] = trids[t][i]
            prow['x'], prow['y'], prow['z'], prow['r'] = part[i]
            prow.append()
        pos.flush()
    h5file.flush()

def fill_geometry12(h5file, sample_group):
    for time_step in h5file.walkGroups(sample_group):
        if not hasattr(time_step._v_attrs, 't') or hasattr(time_step, 'neighbours12'):
            continue
        #load positions
        pos = np.column_stack([time_step.positions.col(c)[:] for c in 'xyzr'])
        neighbours, inside = colloids.particles.get_N_ngbs(pos[:,:-1], pos[:,-1])
        #create geometry table
        table = h5file.createTable(
            time_step, 'neighbours12',
            bonds12, 'the first 12 neighbours of each particle'
        )
        table.flush()
        row = table.row
        for i, (ins, ngbs) in enumerate(zip(inside, neighbours)):
            row['idnumber'] = i
            row['inside'] = ins
            row['ngbs'] = ngbs
            row.append()
        table.flush()
    h5file.flush()

def fill_boo12_l(h5file, sample_group, l=4):
    class Bond_order_l(IsDescription):
        qlm       = ComplexCol(16, pos=1, shape=l+1)
        Qlm       = ComplexCol(16, pos=2, shape=l+1)
        ql        = Float64Col(pos=3)
        wl        = Float64Col(pos=4)
        Ql        = Float64Col(pos=5)
    for time_step in h5file.walkGroups(sample_group):
        if not hasattr(time_step._v_attrs, 't') or hasattr(time_step, 'boo12_l%d'%l):
            continue
        coords = np.column_stack([time_step.positions.col(c)[:] for c in 'xyz'])
        inside = time_step.neighbours12.cols.inside[:]
        ngbs = time_step.neighbours12.cols.ngbs[:]
        qlms = boo.weave_qlm(coords, ngbs, inside, l)
        Qlms = np.copy(qlms)
        code = """
        #pragma omp parallel for
        for(int i=0; i<Nqlms[0]; ++i)
        {
            if(!inside(i)) continue;
            int nb = 0;
            for(int j=0; j<Nngbs[1]; ++j)
            {
                int q = ngbs(i, j);
                if(q<0 || q>=Nqlms[0] || !inside(q)) continue;
                nb++;
                for(int m=0; m<Nqlms[1]; ++m)
                    Qlms(i, m) += qlms(q,m);
            }
            for(int m=0; m<Nqlms[1]; ++m)
                Qlms(i, m) /= nb+1;
        }
        """
        weave.inline(
            code,['Qlms', 'ngbs', 'inside', 'qlms'],
            type_converters =converters.blitz,
            extra_compile_args =['-O3 -fopenmp'],
            extra_link_args=['-lgomp'],
            verbose=2, compiler='gcc')
        qls = boo.ql(qlms)
        wls = boo.wl(qlms)
        Qls = boo.ql(Qlms)
        #bond orientational order table
        table = h5file.createTable(
            time_step, 'boo12_l%d'%l,
            Bond_order_l, '%d-fold symmetry'%l
        )
        table.flush()
        row = table.row
        for i, (qlm, Qlm, ql, wl, Ql) in enumerate(zip(qlms, Qlms, qls, wls, Qls)):
            row['qlm'] = qlm
            row['Qlm'] = Qlm
            row['ql'] = ql
            row['wl'] = wl
            row['Ql'] = Ql
            row.append()
        table.flush()
    h5file.flush()
    
def periodic_ngb12(positions, radii, L):
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
    ngbs = np.zeros([len(positions), 12], int)
    code = """
    //spatial indexing
    RTree tree;
    for(int p=0; p<Npositions[0]; ++p)
    {
        typename RTree::BoundingBox bb;
		for(int d=0; d<3; ++d)
		{
			bb.edges[d].first = positions(p,d) - radii(p);
			bb.edges[d].second = positions(p,d) + radii(p);
		}
        tree.Insert(p, bb);
    }
    //look for neighbours
    for(int p=0; p<Npositions[0]; ++p)
    {
        double rsq = 9.0*radii(p)*radii(p);
        std::list<int> overlapping;
        Gatherer ga(overlapping);
        typename RTree::BoundingBox bb;
        for(int i=-1; i<2; ++i)
        {
            bb.edges[0].first = positions(p,0) - r + i*L;
            bb.edges[0].second = positions(p,0) + r + i*L;
            for(int j=-1; j<2; ++j)
            {
                bb.edges[1].first = positions(p,1) - r + j*L;
	            bb.edges[1].second = positions(p,1) + r + j*L;
                for(int k=-1; k<2; ++k)
	            {
		            bb.edges[2].first = positions(p,2) - r + k*L;
		            bb.edges[2].second = positions(p,2) + r + k*L;
		            tree.Query(typename RTree::AcceptOverlapping(bb), ga);
	            }
	        }
        }
        overlapping.sort();
        overlapping.unique();
        std::multimap<double, int> bysqdist;
        for(std::list<int>::const_iterator it=overlapping.begin(); it!=overlapping.end(); ++it)
        {
            const int j = *it;
            if(j==i)
                continue;
            double dsq = 0;
            for(int d=0; d<3; ++d)
                dsq += pow(fmod(positions(i,d)-positions(j,d)+1.5*L, L), 2);
            if(dsq<rsq)
                bysqdist.insert(std::make_pair(dsq, j));
        }
	    assert(bysqdist.size()>=12);
        for(int n=0; n<12; ++n)
        {
            ngbs(p,n) = bysqdist.begin()->second;
            bysqdist.erase(bysqdist.begin());
        }
    }
    """
    weave.inline(
        code,['positions', 'radii', 'ngbs', 'L'],
        type_converters =converters.blitz,
        support_code = support,
        include_dirs = ['/home/mathieu/src/colloids/multiscale/RStarTree'],
        headers = ['"RStarTree.h"','<map>', '<list>'],
        extra_compile_args =['-O3 -fopenmp'],
        extra_link_args=['-lgomp'],
        verbose=2, compiler='gcc')
    return ngbs
    

def fill_nXbonds(h5file, sample_group, threshold=0.7):
    for time_step in h5file.walkGroups(sample_group):
        if not hasattr(time_step._v_attrs, 't'):
            continue
        nXbonds = np.asarray([
            np.sum(colloids.particles.boo_normed_product(
                time_step.boo12.cols.q6m[i],
                time_step.boo12.readCoordinates(ngbs[ngbs>-1], "q6m")
                )>threshold)
            for i, ngbs in enumerate(time_step.neighbours12.cols.ngbs)
            ])
        h5file.createArray(
            time_step, 'nXbonds',
            nXbonds, 'number of crystalline bonds per particles'
            )
        time_step.nXbonds._v_attrs.threshold=threshold
    h5file.flush()
        

def fill_ROI(h5file, sample_group, zmin, zmax):
    roi = h5file.createGroup(
        sample_group, "ROI",
        "Region of interest, discarding walls and sedimented regions"
        )
    roi._v_attrs.zmin = zmin
    roi._v_attrs.zmax = zmax
    inside = h5file.createVLArray(
        roi, 'inside',
        Int64Atom(shape=()),
        "ID of included position in each time step")
    for t in range(sample_group._v_attrs.n_time_steps):
        time_step = getattr(sample_group, 't%03d'%t)
        inslab = time_step.positions.getWhereList(
            '(z<=%f) & (z>=%f)'%(roi._v_attrs.zmax, roi._v_attrs.zmin)
            )
        inside.append(inslab[
            time_step.neighbours12.readCoordinates(inslab, 'inside')
            ])
    inside.flush()

def fill_ROI_traj(h5file, sample_group, tmin, tmax):
    roi = sample_group.ROI
    roi._v_attrs.tmin = tmin
    roi._v_attrs.tmax = tmax
    #select trajectories spanning [tmin, tmax]
    trajs = np.asarray([
        tr[tmin+1-tr[0]:tmax-tr[0]+2] for tr in sample_group.trajectories.iterrows()
        if tr[0]<=tmin and len(tr)-1+tr[0]>=tmax+1
        ])
    #remove trajectories in the margin
    trajs = trajs[np.min([
        getattr(sample_group, 't%03d'%(t+tmin)).neighbours12.readCoordinates(tr, 'inside')
        for t, tr in enumerate(trajs.T)
        ], axis=0)]
    #remove trajectories that get out of the z slab
    z = np.asarray([
        getattr(sample_group, 't%03d'%(t+tmin)).positions.readCoordinates(tr, 'z')
        for t, tr in enumerate(trajs.T)]
        )
    trajs = trajs[numexpr.evaluate(
        '(m<=zmax) & (M>=zmin)', {
            'm': z.min(0),
            'M': z.max(0),
            'zmin': roi._v_attrs.zmin,
            'zmax': roi._v_attrs.zmax
            })]
    h5file.createArray(
        sample_group.ROI, 'trajectories',
        trajs, 'trajectories spanning [%d, %d] in the ROI'%(tmin, tmax)
        )

def fill_dynamics(h5file, sample_group):
    roi = sample_group.ROI
    trajs = roi.trajectories[:]
    positions = np.asarray([
        np.column_stack([getattr(
            sample_group,
            't%03d'%(t+roi._v_attrs.tmin)
            ).positions.readCoordinates(tr, c) for c in 'xyz'])
        for t, tr in enumerate(trajs.T)])
    radii = sample_group.t000.positions.readCoordinates(tr, 'r')*np.sqrt(2)
    #remove drift
    drift = np.mean(positions-positions[0], 1)
    positions -= drift[:,None,:]
    #Self isf
    A = numexpr.evaluate(
        'exp(complex(0, a) * pi / r)',
        {'a':positions, 'pi':np.pi, 'r':radii[None,:,None]})
    isf = sta.time_correlation(A,0).mean(axis=-1).real
    h5file.createArray(
        sample_group.ROI, 'isf',
        isf, 'Ensemble averaged self intermediate scattering function.'
        )
    #MSD and non-Gaussian parameter
    msd = np.zeros(len(positions))
    mqd = np.zeros(len(positions))
    for t0, p in enumerate(positions):
        for dt, q in enumerate(positions[t0+1:]):
            diff = numexpr.evaluate(
                '((q-p)/r)**2',
                {'p':p, 'q':q, 'r':radii[None,:,None]}
                ).sum(axis=-1)
            msd[dt+1] += diff.sum()
            mqd[dt+1] += (diff**2).sum()
    mqd *= (len(mqd) - np.arange(len(mqd))) * positions.shape[1]
    mqd[1:] = numexpr.evaluate(
        '3 * a / (5 * b**2) - 1',
        {'a': mqd[1:], 'b': msd[1:]})
    mqd[0] = 0
    h5file.createArray(
        sample_group.ROI, 'ngp',
        mqd, 'Evolution of the non Gaussian parameter.'
        )
    msd /= (len(mqd) - np.arange(len(mqd))) * positions.shape[1]
    h5file.createArray(
        sample_group.ROI, 'msd',
        msd, 'Evolution of the non Mean Square Displacement.'
        )

def fill_S(h5file, sample_group, shape=[256]*3):
    im = np.zeros(shape)
    #mask of all the "distances" in (half)fourier space
    dists = np.fft.fftshift(np.sqrt(
        (np.arange(-shape[0]/2, shape[0]/2, dtype=np.float32)**2)[:,None,None] + 
        (np.arange(-shape[1]/2, shape[1]/2, dtype=np.float32)**2)[None,:,None] + 
        (np.arange(shape[2]/2+1, dtype=np.float32)**2)[None,None,:]
        ), [0,1])
    #do not take pure directions into account to avoid window harmonics
    dists[0]=-1
    dists[:,0]=-1
    dists[:,:,0]=-1
    roi = sample_group.ROI
    nbtot = 0
    S = np.zeros(min(shape))
    for t, inside in enumerate(roi.inside.iterrows()):
        if t<roi._v_attrs.tmin or t>roi._v_attrs.tmax:
            continue
        pos = np.column_stack([getattr(getattr(sample_group, 't%03d'%t).positions.cols, c)[:] for c in 'xyz'])
        nbtot += len(pos)
        #draw a white pixel at the position of each particle
        im.fill(0)
        for x, y, z in pos0:
            im[x,y,z] = 1
        #do the (half)Fourier transform
        spectrum = numexpr.evaluate('abs(a).real**2', {'a':anfft.rfftn(im, 3, measure=True)})
        #radial average (sum)
        S += np.histogram(dists.ravel(), np.arange(len(S)+1), weights=spectrum.ravel())[0]
    #normalize by the total number of particles
    S /= nbtot
    nb = np.histogram(dists.ravel(), np.arange(len(S4)+1))[0]
    #radial average (division)
    S[nb>0] /= nb[nb>0]
    S_array = h5file.createArray(
        sample_group.ROI, 'Sq',
        S, 'Structure factor'
        )
    S_array._v_attrs.bins = np.arange(len(S))
    
def fill_G6(h5file, sample_group, bins=np.linspace(0, 256., 1024.)):
    maxsq = 1.0 * bins[-1]**2
    G6_table = np.zeros([sample_group._v_attrs.n_time_steps, len(bins)-1, 2])
    for t, inside in enumerate(sample_group.ROI.inside.iterrows()):
        time_step = getattr(sample_group, 't%03d'%t)
        ngbs = time_step.neighbours12.cols.ngbs[:]
        isin = time_step.neighbours12.cols.inside[:][ngbs[inside]].min(-1)
        pos = np.column_stack([
            time_step.positions.readCoordinates(inside[isin], c)
            for c in 'xyz'])
        if hasattr(time_step, 'boo12'):
            Q6ms = time_step.boo12.readCoordinates(inside[isin], 'Q6m')
        else:
            Q6ms = time_step.boo12_l6.readCoordinates(inside[isin], 'Qlm')
        htot = np.zeros(len(bins)-1)
        gtot = np.zeros(len(bins)-1, int)
        code="""
        #pragma omp parallel for
        for(int i=0; i<Npos[0]; ++i)
        {
            blitz::Array<int,1> g(Ngtot[0]);
            g=0;
            blitz::Array<double,1> h(Nhtot[0]);
            h=0;
            for(int j=0; j<Npos[0]; ++j)
            {
                double disq = 0.0;
                for(int dim=0; dim<3;++dim)
                    disq += pow(pos(i,dim)-pos(j,dim), 2);
                if(disq>= (double)maxsq)
                    continue;
                const int d = sqrt(disq/(double)maxsq)*Nhtot[0];
                ++g(d);
                double pq = real(Q6ms(i,0)*conj(Q6ms(j,0)));
                for(int m=1; m<NQ6ms[1]; ++m)
                    pq += 2.0*real(Q6ms(i,m)*conj(Q6ms(j,m)));
                h(d) += pq;
            }
            #pragma omp critical
            {
                for(int b=0; b<Ngtot[0]; ++b)
                    if(g(b)>0)
                        htot(b) += h(b)/g(b);
                gtot += g;
            }
        }
        """
        weave.inline(
            code,['Q6ms', 'pos', 'maxsq', 'htot', 'gtot'],
            type_converters =converters.blitz,
            extra_compile_args =['-O3 -fopenmp'],
            extra_link_args=['-lgomp'],
            verbose=2, compiler='gcc')
        G6_table[t] = np.column_stack((htot, gtot))/ len(pos)
    G6_array = h5file.createArray(
        sample_group.ROI, 'G6',
        G6_table, 'Spatial correlation of the Q6m'
        )
    G6_array._v_attrs.bins = bins

def fill_Glgl(h5file, sample_group, maxdist=256.0, Nbins =1024, l=6):
    maxsq = float(maxdist**2)
    G6_table = np.zeros([sample_group._v_attrs.n_time_steps, Nbins, 3])
    for t, inside in enumerate(sample_group.ROI.inside.iterrows()):
        time_step = getattr(sample_group, 't%03d'%t)
        pos = np.column_stack((
            time_step.positions.readCoordinates(inside, 'x'),
            time_step.positions.readCoordinates(inside, 'y'),
            time_step.positions.readCoordinates(inside, 'z')
            ))
        qlms = getattr(time_step, 'boo12_l%d'%l).readCoordinates(inside, 'qlm')[:]
        Qlms = getattr(time_step, 'boo12_l%d'%l).readCoordinates(inside, 'Qlm')[:]
        hQtot = np.zeros(Nbins)
        hqtot = np.zeros(Nbins)
        gtot = np.zeros(Nbins, int)
        code = """
        #pragma omp parallel for
        for(int i=0; i<Npos[0]; ++i)
        {
            blitz::Array<int,1> g(Nbins);
            g = 0;
            blitz::Array<double,1> hq(Nbins);
            hq = 0.0;
            blitz::Array<double,1> hQ(Nbins);
            hQ = 0.0;
            for(int j=0; j<Npos[0]; ++j)
            {
                if(i==j) continue;
                double disq = 0.0;
                for(int dim=0; dim<3;++dim)
                    disq += pow(pos(i,dim)-pos(j,dim), 2);
                if(disq>=(double)maxsq)
                    continue;
                const int d = sqrt(disq/(double)maxsq)*Nbins;
                ++g(d);
                double pq = real(qlms(i,0)*conj(qlms(j,0)));
                for(int m=1; m<Nqlms[1]; ++m)
                    pq += 2.0*real(qlms(i,m)*conj(qlms(j,m)));
                pq *= 4.0*M_PI/(2.0*(Nqlms[1]-1)+1);
                hq(d) += pq;
                double pQ = real(Qlms(i,0)*conj(Qlms(j,0)));
                for(int m=1; m<NQlms[1]; ++m)
                    pQ += 2.0*real(Qlms(i,m)*conj(Qlms(j,m)));
                pQ *= 4.0*M_PI/(2.0*(NQlms[1]-1)+1);
                hQ(d) += pQ;
            }
            #pragma omp critical
            {
                gtot += g;
                hqtot += hq;
                hQtot += hQ;
            }
            /*for(int b=0; b<Nbins; ++b)
                if(g(b)>0)
                {
                    #pragma omp atomic
                    hqtot(b) += hq(b)/g(b);
                    #pragma omp atomic
                    hQtot(b) += hQ(b)/g(b);
                }*/
        }
        """
        weave.inline(
            code,['qlms', 'Qlms', 'pos', 'maxsq', 'Nbins', 'hQtot', 'hqtot', 'gtot'],
            type_converters =converters.blitz,
            extra_compile_args =['-O3 -fopenmp'],
            extra_link_args=['-lgomp'],
            verbose=2, compiler='gcc')
        #return np.column_stack((hqtot, hQtot))/ gtot
        G6_table[t] = np.column_stack((hqtot, hQtot, gtot))/ len(pos)
    G6_array = h5file.createArray(
        sample_group.ROI, 'gG_l%d'%l,
        G6_table, 'Spatial correlation of the q%dm and the Q%dm'%(l,l)
        )
    G6_array._v_attrs.bins = np.linspace(0, maxdist, Nbins)
    
def fill_Glgl_conservative(h5file, sample_group, maxdist=30.0, Nbins =100, l=6):
    maxsq = float(maxdist**2)
    G6_table = np.zeros([sample_group._v_attrs.n_time_steps, Nbins, 3])
    nbdens =  np.zeros(sample_group._v_attrs.n_time_steps)
    for t, inside in enumerate(sample_group.ROI.inside.iterrows()):
        time_step = getattr(sample_group, 't%03d'%t)
        pos = np.column_stack((
            time_step.positions.readCoordinates(inside, 'x'),
            time_step.positions.readCoordinates(inside, 'y'),
            time_step.positions.readCoordinates(inside, 'z')
            ))
        bounds = np.vstack((pos.min(0)+maxdist, pos.max(0)-maxdist))
        is_center = (pos>bounds[0]).min(1) & (pos<bounds[1]).min(1)
        nb = np.sum(is_center)
        assert nb>0, "maxdist is too long ! No particle selected."
        nbdens[t] = nb/np.prod(bounds[1]-bounds[0])
        if l==6 and hasattr(time_step, 'boo12'):
            qlms = time_step.boo12.readCoordinates(inside, 'q6m')[:]
            Qlms = time_step.boo12.readCoordinates(inside, 'Q6m')[:]
        else:
            qlms = getattr(time_step, 'boo12_l%d'%l).readCoordinates(inside, 'qlm')[:]
            Qlms = getattr(time_step, 'boo12_l%d'%l).readCoordinates(inside, 'Qlm')[:]
        G6_table[t] = np.column_stack(boo.gG_l(pos, qlms, Qlms, is_center, Nbins, maxdist))
        if nb>0:
            G6_table[t]/=nb
    G6_array = h5file.createArray(
        sample_group.ROI, 'gG_l%d_conservative'%l,
        G6_table, 'Spatial correlation of the q%dm and the Q%dm using only particles inside the ROI'%(l,l)
        )
    G6_array._v_attrs.bins = np.linspace(0, maxdist, Nbins)
    G6_array._v_attrs.nbdens = nbdens

def fill_Glgl_extended(h5file, sample_group, maxdist=75.0, Nbins =200, l=6):
    maxsq = float(maxdist**2)
    G6_table = np.zeros([sample_group._v_attrs.n_time_steps, Nbins, 3])
    nbdens =  np.zeros(sample_group._v_attrs.n_time_steps)
    roi = sample_group.ROI
    for t, inside in enumerate(roi.inside.iterrows()):
        time_step = getattr(sample_group, 't%03d'%t)
        #select all particles within maxdist of the ROI
        outslab = time_step.positions.getWhereList(
            '(z<=%f) & (z>=%f)'%(roi._v_attrs.zmax+maxdist, roi._v_attrs.zmin-maxdist)
            )
        #restrict to the particles having BOO data
        outslab = outslab[time_step.neighbours12.readCoordinates(outslab, 'inside')]
        #load the coordinates
        pos = np.column_stack((
            time_step.positions.readCoordinates(outslab, 'x'),
            time_step.positions.readCoordinates(outslab, 'y'),
            time_step.positions.readCoordinates(outslab, 'z')
            ))
        #select particles in the ROI
        is_in = np.zeros(len(time_step.positions), bool)
        is_in[inside] = True
        #and within maxdist from the boundaries of outslab
        bounds = np.vstack((pos.min(0)+maxdist, pos.max(0)-maxdist))
        is_center = (pos>bounds[0]).min(1) & (pos<bounds[1]).min(1) & is_in[outslab]
        nb = np.sum(is_center)
        nbdens[t] = nb/np.prod(bounds[1]-bounds[0])
        if l==6 and hasattr(time_step, 'boo12'):
            qlms = getattr(time_step, 'boo12').readCoordinates(outslab, 'q6m')[:]
            Qlms = getattr(time_step, 'boo12').readCoordinates(outslab, 'Q6m')[:]
        else:
            qlms = getattr(time_step, 'boo12_l%d'%l).readCoordinates(outslab, 'qlm')[:]
            Qlms = getattr(time_step, 'boo12_l%d'%l).readCoordinates(outslab, 'Qlm')[:]
        G6_table[t] = np.column_stack(boo.gG_l(pos, qlms, Qlms, is_center, Nbins, maxdist))
        if nb>0:
            G6_table[t]/=nb
    G6_array = h5file.createArray(
        sample_group.ROI, 'gG_l%d_extended'%l,
        G6_table, 'Spatial correlation of the q%dm and the Q%dm using as centers the particles inside the ROI and within maxdist from the boundary. Take correlation with particles outside the ROI.'%(l,l)
        )
    G6_array._v_attrs.bins = np.linspace(0, maxdist, Nbins)
    G6_array._v_attrs.nbdens = nbdens
    
def fill_steihardtgl_extended(h5file, sample_group, maxdist=75.0, Nbins =200, l=6):
    maxsq = float(maxdist**2)
    G6_table = np.zeros([sample_group._v_attrs.n_time_steps, Nbins, 2])
    nbdens =  np.zeros(sample_group._v_attrs.n_time_steps)
    roi = sample_group.ROI
    for t, inside in enumerate(roi.inside.iterrows()):
        time_step = getattr(sample_group, 't%03d'%t)
        #select all particles within maxdist of the ROI
        outslab = time_step.positions.getWhereList(
            '(z<=%f) & (z>=%f)'%(roi._v_attrs.zmax+maxdist, roi._v_attrs.zmin-maxdist)
            )
        #remove particles on the edges of the experimental window
        outslab = outslab[time_step.neighbours12.readCoordinates(outslab, 'inside')]
        is_outslab = np.zeros(len(time_step.positions), bool)
        is_outslab[outslab] = True
        #select bonds
        ngbs = time_step.neighbours12.readCoordinates(outslab, 'ngbs')
        gr = graph()
        gr.add_nodes(outslab)
        for i, n in zip(outslab, ngbs):
            for j in n:
                if is_outslab[j] and not gr.has_edge((i,j)):
                    gr.add_edge((i,j))
        bonds = np.asarray([[i,j] for i,j in gr.edges() if i<j])
        #load all the coordinates
        pos = np.column_stack([getattr(time_step.positions.cols, c)[:] for c in 'xyz'])
        #select particles in the ROI
        is_in = np.zeros(len(time_step.positions), bool)
        is_in[inside] = True
        #and within maxdist from the boundaries of outslab
        bounds = np.vstack((pos[outslab].min(0)+maxdist, pos[outslab].max(0)-maxdist))
        is_center = (pos>bounds[0]).min(1) & (pos<bounds[1]).min(1) & is_in
        nb = np.sum(is_center)
        nbdens[t] = nb/np.prod(bounds[1]-bounds[0])
        G6_table[t] = np.column_stack(boo.steinhardt_g_l(pos, bonds, is_center, Nbins, maxdist, l=l))
        #return G6_table[t]
        if nb>0:
            G6_table[t]/=nb
    G6_array = h5file.createArray(
        sample_group.ROI, 'steihardtg_l%d_extended'%l,
        G6_table, 'Spatial correlation of the spherical harmonics Y%dm using as centers the bonds inside the ROI and within maxdist from the boundary. Take correlation with bonds outside the ROI.'%(l)
        )
    G6_array._v_attrs.bins = np.linspace(0, maxdist, Nbins)
    G6_array._v_attrs.nbdens = nbdens
    
    
def fill_Sl_2D(h5file, sample_group, l=6, coarse=False, margin=3.0):
    if coarse:
        letter = 'Q'
    else:
        letter = 'q'
    roi = sample_group.ROI
    #find the edges
    edges = np.zeros([len(roi.inside), 2, 3])
    for t, inside in enumerate(roi.inside.iterrows()):
        #if t>roi._v_attrs.tmin: break
        time_step = getattr(sample_group, 't%03d'%t)
        if l==6 and hasattr(time_step, 'boo12'):
            qls = time_step.boo12.readCoordinates(inside, '%s6'%letter)[:]
        else:
            qls = getattr(time_step, 'boo12_l%d'%l).readCoordinates(inside, '%sl'%letter)[:]
        isin = (qls>0)#np.ones(len(inside), bool) #time_step.neighbours12.readCoordinates(inside, 'inside') #& (qls>0)
        if coarse:
            ngbs = time_step.neighbours12.cols.ngbs[:]
            isin[ngbs[inside[isin],-1]==-1] = False
            isin[np.min(ngbs[ngbs[inside[isin]]][:,:,-1] == -1, -1)] = False
        pos = np.column_stack([time_step.positions.readCoordinates(inside[isin], c) for c in 'xyz'])
        edges[t,0] = pos.min(0)
        edges[t,1] = pos.max(0)
    bounds = np.vstack([edges[:,0].max(0)+margin, edges[:,1].min(0)-margin])
    bounds[0,-1] = np.maximum(bounds[0,-1], roi._v_attrs.zmin)
    bounds[1,-1] = np.minimum(bounds[1,-1], roi._v_attrs.zmax)
    shape = np.array(np.floor(bounds[1]-bounds[0]), int)+1
    factors = np.ones(3)
    im = np.zeros(shape, np.complex64)
    spectrum = np.zeros(shape[:-1])
    #mask of wavenumbers
    qx = np.fft.fftfreq(shape[0], d=1/factors[0])
    qy = np.fft.fftfreq(shape[1], d=1/factors[1])
    qz = np.fft.fftfreq(shape[2], d=1/factors[2])
    dists = numexpr.evaluate('sqrt(qx**2+qy**2)', {
        'qx':qx[:,None],
        'qy':qy[None,:]
        })
    #bin the wavenumbers
    nbq, qs = np.histogram(dists.ravel(), qx[:len(qx)/2])
    Sl = np.zeros(nbq.shape)
    #tot = np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin+1, l+1], np.complex128)
    tot = np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin+1])
    nbs = np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin+1], int)
    #window function
    ws = [np.hamming(s) for s in shape]
    for t, inside in enumerate(roi.inside.iterrows()):
        if t<roi._v_attrs.tmin or t>roi._v_attrs.tmax:
            continue
        time_step = getattr(sample_group, 't%03d'%t)
        #select the particles with a defined ql
        if l==6 and hasattr(time_step, 'boo12'):
            qls = getattr(time_step, 'boo12').readCoordinates(inside, '%s6'%letter)[:]
        else:
            qls = getattr(time_step, 'boo12_l%d'%l).readCoordinates(inside, '%sl'%letter)[:]
        hasql = (qls>0)#np.ones(len(inside), bool) #time_step.neighbours12.readCoordinates(inside, 'inside') #& (qls>0)
        if coarse:
            ngbs = time_step.neighbours12.cols.ngbs[:]
            hasql[ngbs[inside[hasql],-1]==-1] = False
            hasql[np.min(ngbs[ngbs[inside[hasql]]][:,:,-1] == -1, -1)] = False
        #load coordinates
        pos = np.column_stack([time_step.positions.readCoordinates(inside[hasql],c)[:] for c in 'xyz'])
        isin = np.min(pos>bounds[0], 1) & np.min(pos<bounds[1], 1)
        pos = pos[isin]
        #load BOO
        #if l==6 and hasattr(time_step, 'boo12'):
         #   qlms = getattr(time_step, 'boo12').readCoordinates(inside[hasql][isin], '%6m'%letter)[:]
        #else:
         #   qlms = getattr(time_step, 'boo12_l%d'%l).readCoordinates(inside[hasql][isin], '%slm'%letter)[:]
        tot[t-roi._v_attrs.tmin] = qls[hasql][isin].sum(0)
        nbs[t-roi._v_attrs.tmin] = len(pos)
        spectrum.fill(0)
        #compute the correlation for each m
        #for m, qlm in enumerate(qlms.T):
        im.fill(0)
        #draw a valued pixel at the position of each particle
        for (x, y, z), qm in zip((pos-bounds[0])*factors, qls[hasql][isin]):#qlm):
            im[x,y,z] = qm
        #if m>0: return im;
        #remove offset
        im -= im.mean()
        #windowing
        for d, w in enumerate(ws):
            im *= w[tuple([None]*d + [slice(None)] + [None]*(2-d))]
        #do the (half)Fourier transform
        #spectrum += (2-(m==0))*np.abs(anfft.fftn(im, 3, measure=True)[:,:,0])**2
        spectrum = np.abs(anfft.fftn(im, 3, measure=True)[:,:,0])**2
        #break
        #radial average (sum)
        Sl += np.histogram(dists.ravel(), qs, weights=spectrum.ravel())[0]/spectrum.mean()*len(pos)
        #return qx, Sl[1:]/nbq[1:]/nbs[0] * 4.*np.pi/(2*l+1)
    #normalize by the total number of particles
    Sl /= nbs.sum() #*(2*l+1)/(4.*np.pi)
    #radial average (division)
    Sl[nbq>0] /= nbq[nbq>0]
    name = 'S_l%d'%l
    if coarse:
        name += "_cg"
    Sl_array = h5file.createArray(
        sample_group.ROI, name,
        Sl, 'Bond order structure factor for %s%dm (qz=0 correlations only)'%(letter,l)
        )
    Sl_array._v_attrs.bins = qs
    Sl_array._v_attrs.tot = tot
    Sl_array._v_attrs.nbs = nbs
    
def fill_Sl_2D_bin(h5file, sample_group, l=6, margin=3.0, thr=None, ratio=None):
    roi = sample_group.ROI
    #find the edges
    edges = np.zeros([len(roi.inside), 2, 3])
    isins = []
    poss = []
    qls = []
    for t, inside in enumerate(roi.inside.iterrows()):
        #if t>roi._v_attrs.tmin: break
        time_step = getattr(sample_group, 't%03d'%t)
        if l==6 and hasattr(time_step, 'boo12'):
            ql = time_step.boo12.readCoordinates(inside, 'Q6')[:]
        else:
            ql = getattr(time_step, 'boo12_l%d'%l).readCoordinates(inside, 'Ql')[:]
        isin = (ql>0)
        ngbs = time_step.neighbours12.cols.ngbs[:]
        isin[ngbs[inside[isin],-1]==-1] = False
        isin[np.min(ngbs[ngbs[inside[isin]]][:,:,-1] == -1, -1)] = False
        isins.append(isin)
        qls.append(ql[isin])
        pos = np.column_stack([time_step.positions.readCoordinates(inside[isin], c) for c in 'xyz'])
        poss.append(pos)
        edges[t,0] = pos.min(0)
        edges[t,1] = pos.max(0)
    bounds = np.vstack([edges[:,0].max(0)+margin, edges[:,1].min(0)-margin])
    bounds[0,-1] = np.maximum(bounds[0,-1], roi._v_attrs.zmin)
    bounds[1,-1] = np.minimum(bounds[1,-1], roi._v_attrs.zmax)
    #discard particles that are not inside
    for t, pos in enumerate(poss):
        isin = np.min(pos>bounds[0], 1) & np.min(pos<bounds[1], 1)
        poss[t] = pos[isin]
        isins[t] = isins[t][isin]
        qls[t] = qls[t][isin]
    #distribution of order parameter
    thrQ6 = np.array([
        np.histogram(ql, np.linspace(0, 0.6, 101))[0] 
        for t, ql in enumerate(qls) 
        if t>=roi._v_attrs.tmin and t<=roi._v_attrs.tmax
        ])
    if ratio is not None:
        #set the threshold in order to have a fixed ratio of ordered particles
        thr = np.linspace(0, 0.6, 101)[np.where(
            np.cumsum(thrQ6.sum(0)[::-1])[::-1] > thrQ6.sum()*ratio
            )[0][-1]]
    if thr is None:
        #threshold maximizing the susceptibility
        thr = np.linspace(0, 0.6, 101)[np.argmax(np.var(np.cumsum(thrQ6[:,::-1],1)[:,::-1],0))]
    #best shape for the FFT
    shape = np.array(np.floor(bounds[1]-bounds[0]), int)+1
    #initialise the structure factor calculator
    Sl = sta.StructureFactor2D(shape)
    nbs = np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin+1], int)
    #window function
    ws = [np.hamming(s) for s in shape]
    for t, inside in enumerate(roi.inside.iterrows()):
        if t<roi._v_attrs.tmin or t>roi._v_attrs.tmax:
            continue
        #compute the structure factor
        Sl(poss[t][qls[t]>thr] - bounds[0])
        #return qx, Sl[1:]/nbq[1:]/nbtot * 4.*np.pi/(2*l+1)
    Sl_array = h5file.createArray(
        sample_group.ROI, 'S_l%d_2D_bin'%l,
        Sl.get_S() / thrQ6.sum(), #normalize by the total number of particles
        'Binary bond order structure factor for Q%d (qz=0 correlations only)'%l
        )
    Sl_array._v_attrs.bins = Sl.qs
    Sl_array._v_attrs.threshold = thr
    Sl_array._v_attrs.distributions = np.var(np.cumsum(thrQ6[:,::-1],1)[:,::-1],0) / thrQ6.sum(1).mean()
    
def get_w6Q6_hist(sample_group):
    bins = [np.linspace(-0.052, 0.052, 101), np.linspace(0, 0.6, 101)]
    h = np.zeros([100, 100], int)
    for t, inside in enumerate(sample_group.ROI.inside.iterrows()):
        time_step = getattr(sample_group, 't%03d'%t)
        Q6 = time_step.boo12.readCoordinates(inside, 'Q6')
        w6 = time_step.boo12.readCoordinates(inside, 'w6')
        h += np.histogram2d(w6, Q6, bins=bins)[0]
    return h/(0.06*0.0052)/h.max()

def gen_spanning(sample_group, dt):
    """Position insides at t0 and t0+dt of the trajectories in the ROI at t0; for each t0 in the ROI"""
    roi = sample_group.ROI
    #load all trajectories in memory
    alltrajs = sample_group.trajectories[:]
    #indexing starting time and ending time
    tr_start_stop = np.asarray([[tr[0], len(tr)] for tr in alltrajs])
    tr_start_stop[:,1] = numexpr.evaluate(
        'l-1-t0',
        {'l':tr_start_stop[:,1], 't0': tr_start_stop[:,0]}
        )
    for t in range(roi._v_attrs.tmin, roi._v_attrs.tmax+1-dt):
        #all the trajectories spanning [t, t+dt]
        spanning = numexpr.evaluate(
            '(start<=t0) & (stop>=t1)',
            {
                't0':t, 't1':t+dt+1,
                'start':tr_start_stop[:,0],
                'stop':tr_start_stop[:,1]
                })
        trajs = np.vstack([
            tr[[t-tr[0]+1, t+dt-tr[0]+1]]
            for tr, sp in zip(alltrajs, spanning)
            if sp
            ])
        #restrict to the trajectories that are in the ROI at t
        inside = np.zeros(len(getattr(sample_group, 't%03d'%t).positions), bool)
        inside[roi.inside[t]] = True
        trajs = trajs[inside[trajs[:,0]]]
        yield (t,trajs)
    
def fill_sd_boo(h5file, sample_group, dt):
    roi = sample_group.ROI
    inside = roi.inside
    bins = {
        'Q6': np.linspace(0, 0.6, 101),
        'w6': np.linspace(-0.052, 0.052, 101),
        'q6': np.linspace(0, 0.7, 101)
        }
    h = np.zeros((3, 100))
    h2 = np.zeros((3, 100))
    hn = np.zeros((3, 100), int)
    for t, good in gen_spanning(sample_group, dt):
        t0 = getattr(sample_group, 't%03d'%t)
        t1 = getattr(sample_group, 't%03d'%(t+dt))
        #Displacement between the two time steps
        d = np.column_stack([numexpr.evaluate(
            'a-b',
            {
                'a': t1.positions.readCoordinates(good[:,1], c),
                'b': t0.positions.readCoordinates(good[:,0], c)}
            ) for c in 'xyz'])
        #Square displacement with drift removed
        sd = np.sum(numexpr.evaluate(
            '(d-drift)**2', {'d':d, 'drift':d.mean(0)}
            ),1)
        #histograms
        for i, (boo, bi) in enumerate(bins.iteritems()):
            values = t0.boo12.readCoordinates(good[:,0], boo)
            h[i] += np.histogram(values, bi, weights=sd)[0]
            h2 += np.histogram(values, bi, weights=sd**2)[0]
            hn[i] += np.histogram(values, bi)[0]
    for i, (boo, bi) in enumerate(bins.iteritems()):
        nonz = hn[i]>0
        msd = h[i].sum()/hn[i].sum()
        sd_boo_array = h5file.createArray(
            roi, 'sd_%s'%boo,
            np.column_stack((
                bi[:-1][nonz],
                h[i][nonz]/hn[i][nonz]/msd,
                hn[i][nonz],
                np.sqrt(h2[i][nonz]/hn[i][nonz] - (h[i][nonz]/hn[i][nonz])**2)/msd
                )),
            'Normalised square displacement function of %s'%boo
            )
        sd_boo_array._v_attrs.dt = dt

def fill_sd_w6Q6(h5file, sample_group, dt):
    roi = sample_group.ROI
    inside = roi.inside
    bins = {
        'Q6': np.linspace(0, 0.6, 101),
        'w6': np.linspace(-0.052, 0.052, 101),
        }
    h = np.zeros((100, 100))
    h2 = np.zeros((100, 100))
    hn = np.zeros((100, 100), int)
    for t, good in gen_spanning(sample_group, dt):
        t0 = getattr(sample_group, 't%03d'%t)
        t1 = getattr(sample_group, 't%03d'%(t+dt))
        #Displacement between the two time steps
        d = np.column_stack([numexpr.evaluate(
            'a-b',
            {
                'a': t1.positions.readCoordinates(good[:,1], c),
                'b': t0.positions.readCoordinates(good[:,0], c)}
            ) for c in 'xyz'])
        #Square displacement with drift removed
        sd = np.sum(numexpr.evaluate(
            '(d-drift)**2', {'d':d, 'drift':d.mean(0)}
            ),1)
        #histograms
        values = dict([
            (boo, t0.boo12.readCoordinates(good[:,0], boo))
            for (boo, bi) in bins.iteritems()
            ])
        h += np.histogram2d(
            values['w6'], values['Q6'],
            (bins['w6'], bins['Q6']),
            weights=sd)[0]
        h2 += np.histogram2d(
            values['w6'], values['Q6'],
            (bins['w6'], bins['Q6']),
            weights=sd**2)[0]
        hn += np.histogram2d(
            values['w6'], values['Q6'],
            (bins['w6'], bins['Q6']))[0]
    nonz = hn>0
    msd = h.sum()/hn.sum()
    h[nonz] /= hn[nonz]
    h2[nonz] /= hn[nonz]
    sd_w6Q6_array = h5file.createArray(
        roi, 'sd_w6Q6',
        np.dstack((
            h/msd,
            hn,
            np.sqrt(h2 - h**2)/msd)),
        'Normalised square displacement maped over the (w6, Q6)-plane'
        )
    sd_w6Q6_array._v_attrs.dt = dt
    sd_w6Q6_array._v_attrs.bins = bins
  
def get_overlap(a, b, thr):
    assert a.shape == b.shape
    overlap = np.zeros(len(a), bool)
    weave.blitz('overlap = sum(pow2(a-b), secondIndex())<thr*thr')
    return overlap
    
def fill_zoverlap(h5file, sample_group, dt, Nbins=10, over_thr=4.0):
    roi = sample_group.ROI
    zbins = np.linspace(roi._v_attrs.zmin, roi._v_attrs.zmax, Nbins+1)
    zoverlapping = np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, Nbins], int)
    zNb = np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, Nbins], int)
    #load all trajectories in memory
    alltrajs = sample_group.trajectories[:]
    #indexing starting time and ending time SLOW
    tr_start_stop = np.column_stack(([tr[0] for tr in alltrajs], map(len, alltrajs)))
    weave.blitz('tr_start_stop[:,1] += tr_start_stop[:,0]-1')
    for t, inside in enumerate(roi.inside.iterrows()):
        if t<roi._v_attrs.tmin or t+dt>roi._v_attrs.tmax:
            continue
        time_step = getattr(sample_group, 't%03d'%t)
        time_step1 = getattr(sample_group, 't%03d'%(t+dt))
        #restrict to the particles having a future of dt
        trnumber = time_step.positions.readCoordinates(inside, 'trnumber')
        haveafuture = (t + dt < tr_start_stop[trnumber,1])
        outslab = inside[haveafuture]
        trnumber = trnumber[haveafuture]
        #select future positions
        future = np.asarray([alltrajs[tr][1+t-tr_start_stop[tr,0]+dt] for tr in trnumber])
        #load the coordinates
        pos = np.column_stack([
            time_step.positions.readCoordinates(outslab, c)
            for c in 'xyz'])
        pos1 = np.column_stack([
            time_step1.positions.readCoordinates(future, c)
            for c in 'xyz'])
        #remove drift
        pos1 -= np.mean(pos1-pos, 0)
        zoverlapping[t-roi._v_attrs.tmin] = np.histogram(
            pos[:,-1], zbins, 
            weights=np.array(numexpr.evaluate(
                'sum((a-b)**2, -1)', 
                {'a': pos, 'b':pos1}
                )<over_thr**2, int)
            )[0]
        zNb[t-roi._v_attrs.tmin] = np.histogram(pos[:,-1], zbins)[0]
    overlap_array = h5file.createArray(
        sample_group.ROI, 'overlap_dt%d'%dt,
        zoverlapping, 'Number of ROI particles overlapping between t and t+dt (dt=%d) (z subdivisions).'%(dt)
        )
    overlap_array._v_attrs.zbins = zbins
    overlap_array._v_attrs.zNbs = zNb
     
def fill_G_overlap_extended(h5file, sample_group, dt, maxdist=75.0, Nbins =200, over_thr=4.0):
    maxsq = float(maxdist**2)
    roi = sample_group.ROI
    zbins = np.linspace(roi._v_attrs.zmin-maxdist, roi._v_attrs.zmax+maxdist, 11)
    Goverlap_table = np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, Nbins], int)
    nbdens =  np.zeros(roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1)
    overlapping =  np.zeros(roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1)
    novZ =  np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, 10], int)
    ncenters =  np.zeros(roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, int)
    ncentersZ = np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, 10], int)
    #load all trajectories in memory
    alltrajs = sample_group.trajectories[:]
    #indexing starting time and ending time SLOW
    tr_start_stop = np.column_stack(([tr[0] for tr in alltrajs], map(len, alltrajs)))
    weave.blitz('tr_start_stop[:,1] += tr_start_stop[:,0]-1')
    for t, inside in enumerate(roi.inside.iterrows()):
        if t<roi._v_attrs.tmin or t+dt>roi._v_attrs.tmax:
            continue
        time_step = getattr(sample_group, 't%03d'%t)
        time_step1 = getattr(sample_group, 't%03d'%(t+dt))
        #select all particles within maxdist of the ROI
        outslab = time_step.positions.getWhereList(
            '(z<=%f) & (z>=%f)'%(roi._v_attrs.zmax+maxdist, roi._v_attrs.zmin-maxdist)
            )
        #remove the particles that are on the edge of the experimental window
        #outslab = outslab[time_step.neighbours12.readCoordinates(outslab, 'inside')]
        #restrict to the particles having a future of dt
        trnumber = time_step.positions.readCoordinates(outslab, 'trnumber')
        haveafuture = (t + dt < tr_start_stop[trnumber,1])
        outslab = outslab[haveafuture]
        trnumber = trnumber[haveafuture]
        #select future positions
        goodtrajs = np.zeros(len(alltrajs), bool)
        goodtrajs[trnumber] = True
        future = np.asarray([alltrajs[tr][1+t-tr_start_stop[tr,0]+dt] for tr in trnumber])
        #load the coordinates
        pos = np.column_stack([
            time_step.positions.readCoordinates(outslab, c)
            for c in 'xyz'])
        pos1 = np.column_stack([
            time_step1.positions.readCoordinates(future, c)
            for c in 'xyz'])
        #remove drift
        pos1 -= np.mean(pos1-pos, 0)
        ncentersZ[t-roi._v_attrs.tmin] = np.histogram(pos[:,2], bins=zbins)[0]
        #compute the overlap
        overlap = numexpr.evaluate('sum((a-b)**2, -1)', {'a': pos, 'b':pos1})<over_thr**2
        overlapping[t-roi._v_attrs.tmin] = overlap.sum() / float(len(overlap))
        if overlap.sum()<2:
            continue
        #restrict to slow particles
        pos = pos[overlap]
        outslab = outslab[overlap]
        novZ[t-roi._v_attrs.tmin] = np.histogram(pos[:,2], bins=zbins)[0]
        #select particles in the ROI
        is_in = np.zeros(len(time_step.positions), bool)
        is_in[inside] = True
        #and within maxdist from the boundaries of outslab
        bounds = np.vstack((pos.min(0)+maxdist, pos.max(0)-maxdist))
        is_center = (pos>bounds[0]).min(1) & (pos<bounds[1]).min(1) & is_in[outslab]
        ncenters[t-roi._v_attrs.tmin] = is_center.sum()
        nbdens[t-roi._v_attrs.tmin] = (len(pos)-1)/np.prod(pos.ptp(0))
        Goverlap_table[t-roi._v_attrs.tmin] = colloids.particles.get_rdf(pos, is_center, Nbins, maxdist)
    Goverlap_array = h5file.createArray(
        sample_group.ROI, 'Goverlap_dt%d'%dt,
        Goverlap_table, 'Radial distribution function of the particles overlapping between t and t+dt (dt=%d) using as centers the particles inside the ROI and within maxdist from the boundary. Take correlation with particles outside the ROI.'%(dt)
        )
    Goverlap_array._v_attrs.bins = np.linspace(0, maxdist, Nbins+1)
    Goverlap_array._v_attrs.nbdens = nbdens
    Goverlap_array._v_attrs.overlapping = overlapping
    Goverlap_array._v_attrs.ncenters = ncenters
    Goverlap_array._v_attrs.novZ = novZ
    Goverlap_array._v_attrs.ncentersZ = ncentersZ
    
def fill_G_overlap_planar(h5file, sample_group, dt, maxdist=75.0, Nbins =200, over_thr=4.0, maxangle=np.pi/3):
    maxsq = float(maxdist**2)
    roi = sample_group.ROI
    Zmaxdist = np.sin(maxangle) * maxdist
    zbins = np.linspace(roi._v_attrs.zmin-Zmaxdist, roi._v_attrs.zmax+Zmaxdist, 11)
    Goverlap_table = np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, Nbins], int)
    nbdens =  np.zeros(roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1)
    overlapping =  np.zeros(roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1)
    novZ =  np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, 10], int)
    ncenters =  np.zeros(roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, int)
    ncentersZ = np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, 10], int)
    #load all trajectories in memory
    alltrajs = sample_group.trajectories[:]
    #indexing starting time and ending time SLOW
    tr_start_stop = np.column_stack(([tr[0] for tr in alltrajs], map(len, alltrajs)))
    weave.blitz('tr_start_stop[:,1] += tr_start_stop[:,0]-1')
    for t, inside in enumerate(roi.inside.iterrows()):
        if t<roi._v_attrs.tmin or t+dt>roi._v_attrs.tmax:
            continue
        time_step = getattr(sample_group, 't%03d'%t)
        time_step1 = getattr(sample_group, 't%03d'%(t+dt))
        #select all particles within Zmaxdist of the ROI
        outslab = time_step.positions.getWhereList(
            '(z<=%f) & (z>=%f)'%(roi._v_attrs.zmax+Zmaxdist, roi._v_attrs.zmin-Zmaxdist)
            )
        #restrict to the particles having a future of dt
        trnumber = time_step.positions.readCoordinates(outslab, 'trnumber')
        haveafuture = (t + dt < tr_start_stop[trnumber,1])
        outslab = outslab[haveafuture]
        trnumber = trnumber[haveafuture]
        #select future positions
        goodtrajs = np.zeros(len(alltrajs), bool)
        goodtrajs[trnumber] = True
        future = np.asarray([alltrajs[tr][1+t-tr_start_stop[tr,0]+dt] for tr in trnumber])
        #load the coordinates
        pos = np.column_stack([
            time_step.positions.readCoordinates(outslab, c)
            for c in 'xyz'])
        pos1 = np.column_stack([
            time_step1.positions.readCoordinates(future, c)
            for c in 'xyz'])
        #remove drift
        pos1 -= np.mean(pos1-pos, 0)
        ncentersZ[t-roi._v_attrs.tmin] = np.histogram(pos[:,2], bins=zbins)[0]
        #compute the overlap
        overlap = numexpr.evaluate('sum((a-b)**2, -1)', {'a': pos, 'b':pos1})<over_thr**2
        overlapping[t-roi._v_attrs.tmin] = overlap.sum() / float(len(overlap))
        if overlap.sum()<2:
            continue
        #restrict to slow particles
        pos = pos[overlap]
        outslab = outslab[overlap]
        novZ[t-roi._v_attrs.tmin] = np.histogram(pos[:,2], bins=zbins)[0]
        #select particles in the ROI
        is_in = np.zeros(len(time_step.positions), bool)
        is_in[inside] = True
        #and within maxdist from the boundaries of outslab
        bounds = np.vstack((pos.min(0)+[maxdist, maxdist, Zmaxdist], pos.max(0)-[maxdist, maxdist, Zmaxdist]))
        is_center = (pos>bounds[0]).min(1) & (pos<bounds[1]).min(1) & is_in[outslab]
        ncenters[t-roi._v_attrs.tmin] = is_center.sum()
        nbdens[t-roi._v_attrs.tmin] = (len(pos)-1)/np.prod(pos.ptp(0))
        Goverlap_table[t-roi._v_attrs.tmin] = colloids.particles.get_planar_rdf(pos, is_center, Nbins, maxdist, maxangle)
    Goverlap_array = h5file.createArray(
        sample_group.ROI, 'GoverlapXY_dt%d'%dt,
        Goverlap_table, 'Planar radial distribution function of the particles overlapping between t and t+dt (dt=%d) using as centers the particles inside the ROI and within maxdist from the boundary. Take correlation with particles outside the ROI.'%(dt)
        )
    Goverlap_array._v_attrs.bins = np.linspace(0, maxdist, Nbins+1)
    Goverlap_array._v_attrs.nbdens = nbdens
    Goverlap_array._v_attrs.overlapping = overlapping
    Goverlap_array._v_attrs.ncenters = ncenters
    Goverlap_array._v_attrs.novZ = novZ
    Goverlap_array._v_attrs.ncentersZ = ncentersZ
    Goverlap_array._v_attrs.maxangle = maxangle
    
def fill_G_overlap_slice(h5file, sample_group, dt, maxdist=75.0, Nbins =200, over_thr=4.0, width=10.0, margin=10.0):
    maxsq = float(maxdist**2)
    roi = sample_group.ROI
    Zmaxdist = 0.5*width
    zbins = np.linspace(roi._v_attrs.zmin-Zmaxdist, roi._v_attrs.zmax+Zmaxdist, 11)
    Goverlap_table = np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, Nbins], int)
    nbdens =  np.zeros(roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1)
    overlapping =  np.zeros(roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1)
    novZ =  np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, 10], int)
    ncenters =  np.zeros(roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, int)
    ncentersZ = np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, 10], int)
    #load all trajectories in memory
    alltrajs = sample_group.trajectories[:]
    #indexing starting time and ending time SLOW
    tr_start_stop = np.column_stack(([tr[0] for tr in alltrajs], map(len, alltrajs)))
    weave.blitz('tr_start_stop[:,1] += tr_start_stop[:,0]-1')
    for t, inside in enumerate(roi.inside.iterrows()):
        if t<roi._v_attrs.tmin or t+dt>roi._v_attrs.tmax:
            continue
        time_step = getattr(sample_group, 't%03d'%t)
        time_step1 = getattr(sample_group, 't%03d'%(t+dt))
        #select all particles within Zmaxdist of the ROI
        outslab = time_step.positions.getWhereList(
            '(z<=%f) & (z>=%f)'%(roi._v_attrs.zmax+Zmaxdist, roi._v_attrs.zmin-Zmaxdist)
            )
        #remove the particles that are on the edge of the experimental window
        #outslab = outslab[time_step.neighbours12.readCoordinates(outslab, 'inside')]
        #restrict to the particles having a future of dt
        trnumber = time_step.positions.readCoordinates(outslab, 'trnumber')
        haveafuture = (t + dt < tr_start_stop[trnumber,1])
        outslab = outslab[haveafuture]
        trnumber = trnumber[haveafuture]
        #select future positions
        goodtrajs = np.zeros(len(alltrajs), bool)
        goodtrajs[trnumber] = True
        future = np.asarray([alltrajs[tr][1+t-tr_start_stop[tr,0]+dt] for tr in trnumber])
        #load the coordinates
        pos = np.column_stack([
            time_step.positions.readCoordinates(outslab, c)
            for c in 'xyz'])
        pos1 = np.column_stack([
            time_step1.positions.readCoordinates(future, c)
            for c in 'xyz'])
        #discard particles in an outer margin
        inmargin = (pos[:,:2] > pos[:,:2].min(0)+margin).min(1) & (pos[:,:2] < pos[:,:2].max(0)-margin).min(1)
        pos = pos[inmargin]
        pos1 = pos1[inmargin]
        outslab = outslab[inmargin]
        #remove drift
        pos1 -= np.mean(pos1-pos, 0)
        ncentersZ[t-roi._v_attrs.tmin] = np.histogram(pos[:,2], bins=zbins)[0]
        #compute the overlap
        overlap = numexpr.evaluate('sum((a-b)**2, -1)', {'a': pos, 'b':pos1})<over_thr**2
        overlapping[t-roi._v_attrs.tmin] = overlap.sum() / float(len(overlap))
        if overlap.sum()<2:
            continue
        #restrict to slow particles
        pos = pos[overlap]
        outslab = outslab[overlap]
        novZ[t-roi._v_attrs.tmin] = np.histogram(pos[:,2], bins=zbins)[0]
        #select particles in the ROI
        is_in = np.zeros(len(time_step.positions), bool)
        is_in[inside] = True
        #and within maxdist from the boundaries of outslab
        bounds = np.vstack((pos.min(0)+[maxdist, maxdist, Zmaxdist], pos.max(0)-[maxdist, maxdist, Zmaxdist]))
        is_center = (pos>bounds[0]).min(1) & (pos<bounds[1]).min(1) & is_in[outslab]
        ncenters[t-roi._v_attrs.tmin] = is_center.sum()
        nbdens[t-roi._v_attrs.tmin] = (len(pos)-1)/np.prod(pos.ptp(0))
        Goverlap_table[t-roi._v_attrs.tmin] = colloids.particles.get_slice_rdf(pos, is_center, Nbins, maxdist, width)
    Goverlap_array = h5file.createArray(
        sample_group.ROI, 'Goverlap_slice_dt%d'%dt,
        Goverlap_table, 'Radial distribution function of the particles overlapping between t and t+dt (dt=%d) using as centers the particles inside the ROI and within maxdist from the boundary. Take correlation with particles in the same plane (a XY slice of thickness %f).'%(dt,width)
        )
    Goverlap_array._v_attrs.bins = np.linspace(0, maxdist, Nbins+1)
    Goverlap_array._v_attrs.nbdens = nbdens
    Goverlap_array._v_attrs.overlapping = overlapping
    Goverlap_array._v_attrs.ncenters = ncenters
    Goverlap_array._v_attrs.novZ = novZ
    Goverlap_array._v_attrs.ncentersZ = ncentersZ
    Goverlap_array._v_attrs.width = width
    
def fill_S_overlap(h5file, sample_group, dt, over_thr=4.0):
    roi = sample_group.ROI
    #find the edges
    edges = np.zeros([len(roi.inside), 2, 2])
    for t, inside in enumerate(roi.inside.iterrows()):
        pos = np.column_stack([getattr(sample_group, 't%03d'%t).positions.readCoordinates(inside, c) for c in 'xy'])
        edges[t,0] = pos.min(0)
        edges[t,1] = pos.max(0)
    bounds = np.column_stack((
        np.vstack([edges[:,0].min(0), edges[:,1].max(0)]), 
        [roi._v_attrs.zmin, roi._v_attrs.zmax]
        ))
    #bounds = np.vstack([edges[:,0].min(0), edges[:,1].max(0)])
    #select the best shape for FFT
    #shape = 2**(np.array(np.log(bounds[1]-bounds[0]+1)/np.log(2), int)+1)
    #factors = shape/(bounds[1]-bounds[0]+1)
    shape = np.array(np.floor(bounds[1]-bounds[0]), int)+1
    factors = np.ones(3)
    im = np.zeros(shape)
    #mask of wavenumbers
    qx = np.fft.fftfreq(shape[0], d=1/factors[0])
    qy = np.fft.fftfreq(shape[1], d=1/factors[1])
    qz = np.fft.fftfreq(shape[2], d=1/factors[2])
    dists = numexpr.evaluate('sqrt(qx**2+qy**2+qz**2)', {
        'qx':qx[:,None,None],
        'qy':qy[None,:,None],
        'qz':qz[None,None,:shape[2]/2+1]
        })
    #bin the wavenumbers
    nbq, qs = np.histogram(dists.ravel(), qx[:len(qx)/2])
    S4 = np.zeros(nbq.shape)
    #window function
    ws = [np.hamming(s) for s in shape]
    #load all trajectories in memory
    alltrajs = sample_group.trajectories[:]
    #indexing starting time and ending time
    tr_start_stop = np.column_stack(([tr[0] for tr in alltrajs], map(len, alltrajs)))
    weave.blitz('tr_start_stop[:,1] += tr_start_stop[:,0]-1')
    nbtot = 0
    for t, inside in enumerate(roi.inside.iterrows()):
        if t<roi._v_attrs.tmin or t+dt>roi._v_attrs.tmax:
            continue
        pos = np.column_stack([getattr(sample_group, 't%03d'%t).positions.readCoordinates(inside,c)[:] for c in 'xyz'])
        #pos = np.column_stack([getattr(sample_group, 't%03d'%t).positions.col(c)[:] for c in 'xyz'])
        nbtot += len(pos)
        trnumber = getattr(sample_group, 't%03d'%t).positions.readCoordinates(inside, 'trnumber')[:] 
        haveafuture = (t + dt < tr_start_stop[trnumber,1])
        future = np.asarray([alltrajs[tr][1+t-tr_start_stop[tr,0]+dt] for tr in trnumber[haveafuture]])
        pos1 = np.column_stack([getattr(sample_group, 't%03d'%(t+dt)).positions.readCoordinates(future, c)[:] for c in 'xyz'])
        pos1 -= np.mean(pos1-pos[haveafuture], 0)
        overlap = numexpr.evaluate('sum((a-b)**2, -1)', {'a': pos[haveafuture], 'b':pos1})<(over_thr)**2
        pos0 = pos[haveafuture][overlap]
        #draw a white pixel at the position of each slow particle
        im.fill(0)
        #wsum = 0.0
        for x, y, z in (pos0-bounds[0])*factors:
            im[x,y,z] = 1
            #wsum += ws[0][x]*ws[1][y]*ws[2][z]
        #wsum /= len(pos)
        #remove offset
        im -= im.mean()
        #windowing
        for d, w in enumerate(ws):
            im *= w[tuple([None]*d + [slice(None)] + [None]*(2-d))]
        #do the (half)Fourier transform
        spectrum = numexpr.evaluate('real(abs(a))**2', {'a':anfft.rfftn(im, 3, measure=True)})
        #return spectrum, dists
        #radial average (sum)
        S4 += np.histogram(dists.ravel(), qs, weights=spectrum.ravel())[0]/spectrum.mean()*len(pos0)
        #S4[nbq>0] /= nbq[nbq>0]
        #S4 /= nbtot
        #return S4, qs
    #normalize by the total number of particles and NOT by the number of slow particles
    S4 /= nbtot
    #radial average (division)
    S4[nbq>0] /= nbq[nbq>0]
    Soverlap_array = h5file.createArray(
        sample_group.ROI, 'Soverlap_dt%d'%dt,
        S4, 'Structure factor of the particles overlapping between t and t+dt (dt=%d)'%dt
        )
    Soverlap_array._v_attrs.bins = qs
    
def fill_S_overlap2D(h5file, sample_group, dt, over_thr=4.0):
    roi = sample_group.ROI
    #find the edges
    edges = np.zeros([len(roi.inside), 2, 2])
    for t, inside in enumerate(roi.inside.iterrows()):
        pos = np.column_stack([getattr(sample_group, 't%03d'%t).positions.readCoordinates(inside, c) for c in 'xy'])
        edges[t,0] = pos.min(0)
        edges[t,1] = pos.max(0)
    bounds = np.column_stack((
        np.vstack([edges[:,0].min(0), edges[:,1].max(0)]), 
        [roi._v_attrs.zmin, roi._v_attrs.zmax]
        ))
    shape = np.array(np.floor(bounds[1]-bounds[0]), int)+1
    #initialise the structure factor calculator
    S4 = sta.StructureFactor2D(shape)
    #load all trajectories in memory
    alltrajs = sample_group.trajectories[:]
    #indexing starting time and ending time
    tr_start_stop = np.column_stack(([tr[0] for tr in alltrajs], map(len, alltrajs)))
    weave.blitz('tr_start_stop[:,1] += tr_start_stop[:,0]-1')
    nbtot = 0
    for t, inside in enumerate(roi.inside.iterrows()):
        if t<roi._v_attrs.tmin or t+dt>roi._v_attrs.tmax:
            continue
        pos = np.column_stack([getattr(sample_group, 't%03d'%t).positions.readCoordinates(inside,c)[:] for c in 'xyz'])
        #pos = np.column_stack([getattr(sample_group, 't%03d'%t).positions.col(c)[:] for c in 'xyz'])
        nbtot += len(pos)
        trnumber = getattr(sample_group, 't%03d'%t).positions.readCoordinates(inside, 'trnumber')[:] 
        haveafuture = (t + dt < tr_start_stop[trnumber,1])
        future = np.asarray([alltrajs[tr][1+t-tr_start_stop[tr,0]+dt] for tr in trnumber[haveafuture]])
        pos1 = np.column_stack([getattr(sample_group, 't%03d'%(t+dt)).positions.readCoordinates(future, c)[:] for c in 'xyz'])
        pos1 -= np.mean(pos1-pos[haveafuture], 0)
        overlap = numexpr.evaluate('sum((a-b)**2, -1)', {'a': pos[haveafuture], 'b':pos1})<(over_thr)**2
        pos0 = pos[haveafuture][overlap]
        #compute the structure factor
        S4(pos0-bounds[0])
    Soverlap_array = h5file.createArray(
        sample_group.ROI, 'Soverlap2D_dt%d'%dt,
        S4.get_S() / nbtot, #normalize by the total number of particles and NOT by the number of slow particles
        'Structure factor of the particles overlapping between t and t+dt (dt=%d) (qz=0 correlations only)'%dt
        )
    Soverlap_array._v_attrs.bins = S4.qs
    
def fill_S_overlap_slice(h5file, sample_group, dt, shape=[256]*3, over_thr=4.0, width=10):
    #declare c++ code for weave
    histogram_code = """
    #pragma omp parallel for
    for(int z=0; z<Nspectrum[0]; ++z)
        for(int y=0; y<Nspectrum[1]; ++y)
            for(int x=0; x<Nspectrum[2]; ++x)
            {
                const int q = dists(y,x);
                if (q<0 || q>= NS4[1]) continue;
                S4(z, q) += std::norm(spectrum(z,y,x));
            }
    """
    nbtot_code = """
    #pragma omp parallel for
    for(int p=0; p<Npos[0]; ++p)
    {
        blitz::Range r(
            std::max(0, (int)(pos(p,2)-radii(p))),
            std::min(Nnbtot[0], (int)(pos(p,2)+radii(p))+1)-1
            );
        #pragma omp critical
        {
        nbtot(r) += 1;
        }
    }
    """
    im = np.zeros(shape)
    #mask of all the "distances" in (half)fourier space
    dists = np.fft.fftshift(np.sqrt(
        (np.arange(-shape[1]/2, shape[1]/2, dtype=np.float32)**2)[:,None] + 
        (np.arange(0, shape[2]/2+1, dtype=np.float32)**2)[None,:]
        ), [0])
    #do not take pure directions into account to avoid window harmonics
    dists[0]=-1
    dists[:,0]=-1
    roi = sample_group.ROI
    #load all trajectories in memory
    alltrajs = sample_group.trajectories[:]
    #indexing starting time and ending time
    tr_start_stop = np.column_stack(([tr[0] for tr in alltrajs], map(len, alltrajs)))
    weave.blitz('tr_start_stop[:,1] += tr_start_stop[:,0]-1')
    nbtot = np.zeros(shape[0], int)
    S4 = np.zeros([shape[0], min(shape[1:])/2])
    for t in range(roi._v_attrs.tmin, roi._v_attrs.tmax-dt+1):
        time_step = getattr(sample_group, 't%03d'%t)
        pos = np.column_stack([getattr(time_step.positions.cols, c)[:] for c in 'xyz'])
        if width is None:
            radii = time_step.positions.cols.r[:]
        else:
            radii = 0.5*width * np.ones(len(pos))
        weave.inline(
            nbtot_code, ['pos', 'radii', 'nbtot'],
            type_converters =converters.blitz,
            extra_compile_args =['-O3 -fopenmp'],
            extra_link_args=['-lgomp'],
            verbose=2, compiler='gcc')
        trnumber = time_step.positions.cols.trnumber[:] 
        haveafuture = (t + dt < tr_start_stop[trnumber,1])
        future = np.asarray([alltrajs[tr][1+t-tr_start_stop[tr,0]+dt] for tr in trnumber[haveafuture]])
        pos1 = np.column_stack([getattr(sample_group, 't%03d'%(t+dt)).positions.readCoordinates(future, c) for c in 'xyz'])
        pos1 -= np.mean(pos1-pos[haveafuture], 0)
        overlap = numexpr.evaluate('sum((a-b)**2, -1)', {'a': pos[haveafuture], 'b':pos1})<(over_thr)**2
        pos0 = pos[haveafuture][overlap]
        radii = radii[haveafuture][overlap]
        #draw a column of white pixels at the position of each slow particle
        im.fill(0)
        track.draw_rods(pos0, radii, im=im)
        #do the (half)Fourier transform along X and Y
        spectrum = anfft.rfftn(im, 2, measure=True)
        #radial average (sum)
        weave.inline(
            histogram_code,['spectrum', 'dists', 'S4'],
            type_converters =converters.blitz,
            extra_compile_args =['-O3 -fopenmp'],
            extra_link_args=['-lgomp'],
            verbose=2, compiler='gcc')
    #normalize by the total number of particles and NOT by the number of slow particles
    S4 /= nbtot[:,None]
    nb = np.histogram(dists.ravel(), np.arange(S4.shape[1]+1))[0]
    #radial average (division)
    S4[:, nb>0] /= nb[nb>0]
    Soverlap_array = h5file.createArray(
        sample_group.ROI, 'Soverlap_slice_dt%d'%dt,
        S4, 'Structure factor of the particles overlapping between t and t+dt (dt=%d) for each plane'%dt
        )
    Soverlap_array._v_attrs.bins = np.arange(len(S4))
    Soverlap_array._v_attrs.nbtot = nbtot
    
def fill_G_u_extended(h5file, sample_group, dt, maxdist=75.0, Nbins =200):
    maxsq = float(maxdist**2)
    roi = sample_group.ROI
    Gu_table = np.zeros([roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, 2, Nbins], int)
    nbdens =  np.zeros(roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1)
    ncenters =  np.zeros(roi._v_attrs.tmax - roi._v_attrs.tmin -dt+1, int)
    #load all trajectories in memory
    alltrajs = sample_group.trajectories[:]
    #indexing starting time and ending time SLOW
    tr_start_stop = np.add([tr[0]for tr in alltrajs], map(len, alltrajs))
    for t, inside in enumerate(roi.inside.iterrows()):
        if t<roi._v_attrs.tmin or t+dt>roi._v_attrs.tmax:
            continue
        time_step = getattr(sample_group, 't%03d'%t)
        time_step1 = getattr(sample_group, 't%03d'%(t+dt))
        #select all particles within maxdist of the ROI
        outslab = time_step.positions.getWhereList(
            '(z<=%f) & (z>=%f)'%(roi._v_attrs.zmax+maxdist, roi._v_attrs.zmin-maxdist)
            )
        #restrict to the particles having a future of dt
        trnumber = time_step.positions.readCoordinates(outslab, 'trnumber')
        haveafuture = (t + dt < tr_start_stop[trnumber,1])
        outslab = outslab[haveafuture]
        trnumber = trnumber[haveafuture]
        #select future positions
        goodtrajs = np.zeros(len(alltrajs), bool)
        goodtrajs[trnumber] = True
        future = np.asarray([alltrajs[tr][1+t-tr_start_stop[tr,0]+dt] for tr in trnumber])
        #load the coordinates
        pos = np.column_stack([
            time_step.positions.readCoordinates(outslab, c)
            for c in 'xyz'])
        pos1 = np.column_stack([
            time_step1.positions.readCoordinates(future, 'c')
            for c in 'xyz'])
        #remove drift
        pos1 -= (pos1-pos).mean(0)
        #compute the norm of the displacements
        displ = np.sqrt(numexpr.evaluate('sum((a-b)**2, -1)', {'a': pos, 'b':pos1}))
        #fluctuations
        u = displ - displ.mean()
        #select particles in the ROI
        is_in = np.zeros(len(time_step.positions), bool)
        is_in[inside] = True
        #and within maxdist from the boundaries of outslab
        bounds = np.vstack((pos.min(0)+maxdist, pos.max(0)-maxdist))
        is_center = (pos>bounds[0]).min(1) & (pos<bounds[1]).min(1) & is_in[outslab]
        ncenters[t-roi._v_attrs.tmin] = is_center.sum()
        nbdens[t-roi._v_attrs.tmin] = (len(pos)-1)/np.prod(pos.ptp(0))
        Gu_table[t-roi._v_attrs.tmin] = colloids.particles.get_rdf(pos, is_center, Nbins, maxdist)
    Gu_array = h5file.createArray(
        sample_group.ROI, 'Goverlap_dt%d'%dt,
        Gu_table, 'Radial distribution function of the particles overlapping between t and t+dt (dt=%d) using as centers the particles inside the ROI and within maxdist from the boundary. Take correlation with particles outside the ROI.'%(dt)
        )
    Gu_array._v_attrs.bins = np.linspace(0, maxdist, Nbins+1)
    Gu_array._v_attrs.nbdens = nbdens
    Gu_array._v_attrs.ncenters = ncenters

def fill_g4(h5file, sample_group, bins=np.linspace(0, 256, 1024)):
    sigmasq = sample_group._v_attrs.mean_diam**2
    roi = sample_group.ROI
    maxsq = bins[-1]**2
    A = np.zeros([
        roi._v_attrs.tmax+1 - roi._v_attrs.tmin,
        len(bins)-1])
    B = np.zeros(A.shape, int)
    C = np.zeros(A.shape[0])
    D = np.zeros(C.shape, int)
    for dt in range(1, len(A)):
        for t, trajs in gen_spanning(sample_group, dt):
            t0 = getattr(sample_group, 't%03d'%t)
            t1 = getattr(sample_group, 't%03d'%(t+dt))
            p0 = np.column_stack([t0.positions.readCoordinates(trajs[:,0], c) for c in 'xyz'])
            p1 = np.column_stack([t1.positions.readCoordinates(trajs[:,1], c) for c in 'xyz'])
            #Displacement between the two time steps
            d = p1 - p0
            #remove drift and take the norm
            nd = np.sum(numexpr.evaluate(
                '(d-drift)**2', {'d':d, 'drift':d.mean(0)}
                ),1)
            #mobility as defined by Berthier
            mobility = numexpr.evaluate('exp(-nd/sigmasq)')
            #binning of the spatial average
            C[dt] += mobility.sum()
            D[dt] += len(mobility)
            #binning of the spatial correlation
            for i, (p, u) in enumerate(zip(p0, mobility)):
                distsq = numexpr.evaluate(
                    '(p-q)**2',
                    {'p': p, 'q': p0}
                    ).sum(-1)
                inbound = distsq<maxsq
                dist = np.sqrt(distsq[inbound])
                products = u*mobility[inbound]
                A[dt] += np.histogram(dist, bins, weights = products)[0]
                B[dt] = np.histogram(dist, bins)[0]
    #prevent division by 0
    B[B==0] = 1
    D[D==0] = 1
    g4_array = h5file.createArray(
        sample_group.ROI, 'g4',
        numexpr.evaluate(
            'A/B - (C/D)**2',
            {'A':A, 'B':B, 'C':C[:,None], 'D':D[:,None]}
            ),
        'Four-point correlation function'
        )
    g4_array._v_attrs.bins = bins

def get_percol(time_step, centers):
    ngbs = time_step.neighbours12.readCoordinates(centers, 'ngbs')
    gr = graph()
    gr.add_nodes(np.union1d(ngbs.ravel(), centers))
    for i, n in zip(centers, ngbs):
        for j in n:
            if not gr.has_edge((i,j)):
                gr.add_edge((i,j))
    cc = connected_components(gr)
    ncl = max(cc.values())
    sizes = np.histogram(cc.values(), bins=range(1, ncl+2))[0]
    pos = np.column_stack([getattr(time_step.positions.cols, c)[:] for c in 'xyz'])
    clpos = [[] for cl in range(ncl)]
    for p, cl in cc.iteritems():
        clpos[cl-1].append(pos[p])
    rgs = np.asarray([np.var(cl,0).sum(0) for cl in clpos])
    return sizes, rgs
    
def get_percol_iter(time_step, centers, by=10):
    return [get_percol(time_step, centers[m:]) for m in range(0, len(centers), by)]
    
def find_percol(time_step, centers, thr=5000):
    ngbs = time_step.neighbours12.readCoordinates(centers, 'ngbs')
    gr = graph()
    gr.add_nodes(np.union1d(ngbs.ravel(), centers))
    for count, i, n in zip(range(len(centers)), centers, ngbs):
        for j in n:
            if not gr.has_edge((i,j)):
                gr.add_edge((i,j))
        cc = connected_components(gr)
        ncl = max(cc.values())
        sizes = np.histogram(cc.values(), bins=range(1, ncl+2))[0]
        if sizes.max()>thr:
            return count

def exportPOV(
    filename, sample_group, t=None, mrco_thr=0.25, ico_thr=-0.033,
    header='go1.inc'
    ):
    #input
    roi = sample_group.ROI
    if t is None:
        t = roi._v_attrs.tmin
    time_step = getattr(sample_group, 't%03d'%t)
    inside = roi.inside[t]
    #prepare output
    f = pov.File(filename, "colors.inc", header)
    #MRCO
    mrco_centers = inside[
        time_step.boo12.readCoordinates(inside, 'Q6')[:]>mrco_thr
        ]
    mrco_ngbs = np.unique(time_step.neighbours12.readCoordinates(
        mrco_centers, 'ngbs'
        ))
    mrco_all = np.union1d(mrco_centers, mrco_ngbs)
    pov_mrco = [
        pov.Sphere(tuple([row[c] for c in 'xyz']), row['r']*np.sqrt(2))
        for row in time_step.positions.readCoordinates(mrco_all)
        ]
    pov_mrco = pov.Union(*pov_mrco + [pov.Texture(pov.Pigment(color="Green"))])
    f.write(pov_mrco)
    #Icosahedral clusters
    ico_centers = inside[
        time_step.boo12.readCoordinates(inside, 'w6')[:]<ico_thr
        ]
    ico_ngbs = time_step.neighbours12.readCoordinates(ico_centers, 'ngbs')
    ico_all = np.union1d(ico_centers, np.unique(ico_ngbs))
    gr = graph()
    gr.add_nodes(ico_all)
    for i, ngbs in zip(ico_centers, ico_ngbs):
        for n in ngbs:
            edge = (i, n)
            if not gr.has_edge(edge):
                gr.add_edge(edge)
    try:
        cc = connected_components(gr)
    except RuntimeError:
        print "Graph is too large for ico_thr=%g, lower the threshold."%ico_thr
        return
    nb_clusters = max(cc.values())
    pov_ico = [[] for cl in range(nb_clusters)
        ]
    for row in time_step.positions.readCoordinates(ico_all):
        cl = cc[row['idnumber']] - 1
        pov_ico[cl].append(pov.Sphere(
            tuple([row[c] for c in 'xyz']),
            row['r']*np.sqrt(2)
            ))
    for icl, cl in enumerate(pov_ico):
        f.write(pov.Union(*cl + [pov.Texture(pov.Pigment(
            color="COLORSCALE(%f)"%(icl*120.0/nb_clusters)
            ))]))
    f.file.flush()

def fill_voro(h5file, sample_group):
    for t in range(sample_group._v_attrs.n_time_steps):
        time_step = getattr(sample_group, 't%03d'%t)
        pos = np.column_stack([
            time_step.positions.read(field=c) for c in 'xyzr'])
        pos[:,3] *= np.sqrt(2)
        mins = numexpr.evaluate(
            'p-r', {'p':pos[:,:3], 'r':pos[:,3,None]}).min(0)
        maxs = numexpr.evaluate(
            'p+r', {'p':pos[:,:3], 'r':pos[:,3,None]}).max(0)
        pr = subprocess.Popen(
            [
                'voro++', '-r', '-c', '"%i %v %n"', '%g'%(sample_group._v_attrs.mean_diam)
                ]+(' '.join(['%g %g'%b for b in zip(mins, maxs)])).split()+['-'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE
        )
        out = pr.communicate(''.join([
            '%d %f %f %f %f '%(i, p[0], p[1], p[2], p[3])
            for i, p in enumerate(pos)
            ]))[0]
        vol = np.zeros(len(pos))
        ngbs = [[] for i in range(len(pos))]
        for line in out.split("\n"):
            if len(line)==0:
                continue
            l = line[1:-1].split()
            i = int(l[0])
            vol[i] = float(l[1])
            ngbs[i] = map(int, l[2:])
        ngbs_array = h5file.createVLArray(
            time_step, 'voro_ngb',
            Int64Atom(shape=()),
            "Voronoi neighbours. Negative numbers are boundaries.")
        for ngb in ngbs:
            ngbs_array.append(ngb)
        vol_array = h5file.createArray(
            time_step, 'voro_vol',
            vol, 'Volume of the voronoi cell'
            )
def fill_phi_boo(h5file, sample_group):
    roi = sample_group.ROI
    bins = {
        'Q6': np.linspace(0, 0.6, 101),
        'w6': np.linspace(-0.052, 0.052, 101),
        'q6': np.linspace(0, 0.7, 101)
        }
    h = np.zeros((3, 100))
    hn = np.zeros((3, 100), int)
    for t, inside in enumerate(sample_group.ROI.inside.iterrows()):
        if t<roi._v_attrs.tmin or t>roi._v_attrs.tmax:
            continue
        time_step = getattr(sample_group, 't%03d'%t)
        vf = numexpr.evaluate(
            'pi/6.0*(2*sqrt(2)*r)**3/v',
            {
                 'pi':np.pi, 
                 'r':time_step.positions.readCoordinates(inside,'r'),
                 'v':time_step.voro_vol[:][inside]
                 })
        #histograms
        for i, (boo, bi) in enumerate(bins.iteritems()):
            values = time_step.boo12.readCoordinates(inside, boo)
            h[i] += np.histogram(values, bi, weights=vf)[0]
            hn[i] += np.histogram(values, bi)[0]
#    return [np.column_stack((bi[:-1], h[i], hn[i])) for i, (boo, bi) in enumerate(bins.iteritems())]
    for i, (boo, bi) in enumerate(bins.iteritems()):
        nonz = hn[i]>0
        phi_boo_array = h5file.createArray(
            roi, 'phi_%s'%boo,
            np.column_stack((
                bi[:-1][nonz],
                h[i][nonz]/hn[i][nonz],
                hn[i][nonz])),
            'Voronoi volume fraction function of %s'%boo
            )
