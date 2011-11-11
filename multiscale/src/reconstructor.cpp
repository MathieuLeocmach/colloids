/*
 * reconstructor.cpp
 *
 *  Created on: 3 août 2011
 *      Author: mathieu
 */

#include "reconstructor.hpp"
//#include <kdtree++/kdtree.hpp>
#include "multiscalefinder.hpp"
#include <boost/ptr_container/ptr_map.hpp>
#include <math.h>

using namespace std;

namespace Colloids {

Reconstructor::Reconstructor() {
	// TODO Auto-generated constructor stub

}

Reconstructor::~Reconstructor() {
	// TODO Auto-generated destructor stub
}

void Reconstructor::clear()
{
	this->clusters.clear();
	this->trajectories.reset();
	this->last_frame.clear();
}
/**
 * \param tolerance Fraction of the contact distance (sum of radii) accepted. For tolerance<=1 accept overlap only.
 */
void Reconstructor::push_back(const Frame &frame, const double &max_dist)
{
	//remove overlapping
	Frame fr = frame;
	std::auto_ptr<RTree> tree = removeOverlapping<2>(fr);
	if(this->empty())
	{
		this->trajectories.reset(new TrajIndex(fr.size()));
		for(Frame::const_iterator c=fr.begin(); c!=fr.end(); ++c)
			this->clusters.push_back(Cluster(1, Center3D(*c, 0)));
	}
	else
	{
		std::vector<double> distances;
		std::vector<size_t> from, to;
		//use RStarTree spatial indexing to create links
		this->links_by_RStarTree(fr, *tree, distances, from, to, max_dist);
		//brute force
		//this->links_by_brute_force(fr, distances, from, to, tolerance);
		//remember the time step and the number of previously existing trajectories
		const size_t
			t = this->size(),
			old_traj = this->trajectories->size();
		//link trajectories
		this->trajectories->add_Frame(fr.size(), distances, from, to);
		for(size_t p=0; p<fr.size(); ++p)
		{
			const size_t tr = this->trajectories->getTraj(t, p);
			if(tr<old_traj)
				//update formerly existing clusters
				this->clusters[tr].push_back(Center3D(fr[p], t));
			else
				//create new cluster
				this->clusters.push_back(Cluster(1, Center3D(fr[p], t)));
		}
	}
	//keep a copy of the inserted frame
	this->last_frame = fr;
}

void Reconstructor::split_clusters()
{
	const size_t cl_end =  this->clusters.size();
	for(size_t cl=0; cl<cl_end;++cl)
	{
		if(this->clusters[cl].size()<3)
			continue;
		std::vector<double> grad(this->clusters[cl].size()-2);
		Cluster::const_iterator p=this->clusters[cl].begin(), q=this->clusters[cl].begin();
		q++; q++;
		for(size_t i=0; i<this->clusters[cl].size()-2; ++i)
			grad[i] = pow((*p)[0]-(*q)[0], 2) + pow((*p)[1]-(*q)[1], 2);
		//look for local maxima of gradient
		std::list<size_t> blobs;
		for(size_t i=0; i+2<grad.size(); ++i)
			if(grad[i]<=grad[i+1] && grad[i+2]<=grad[i+1] && grad[i+1]>1)
				blobs.push_back(i+2);
		//split, from end to begin
		for(std::list<size_t>::reverse_iterator it=blobs.rbegin(); it!=blobs.rend(); ++it)
		{
			this->clusters.push_back(Cluster());
			Cluster::iterator u = this->clusters[cl].begin();
			std::advance(u, *it);
			this->clusters.back().splice(this->clusters.back().begin(), this->clusters[cl], u, this->clusters[cl].end());
		}
	}
	/*boost::ptr_map<size_t, MultiscaleFinder1D> finders;
	const size_t cl_end =  this->clusters.size();
	for(size_t cl=0; cl<cl_end;++cl)
	{
		if(this->clusters[cl].size()<6)
			continue;
		OctaveFinder::Image signal(1, this->clusters[cl].size());
		//fetch the needed 1D finder
		boost::ptr_map<size_t, MultiscaleFinder1D>::iterator f_it = finders.find(this->clusters[cl].size());
		if(f_it==finders.end())
		{
			//create the finder and cache it
			std::auto_ptr<MultiscaleFinder1D> finder(new MultiscaleFinder1D(this->clusters[cl].size()));
			f_it = finders.insert(this->clusters[cl].size(), finder).first;
		}
		MultiscaleFinder1D &finder = *f_it->second;
		Cluster::const_iterator p=this->clusters[cl].begin(), q=this->clusters[cl].begin();
		q++;
		OctaveFinder::PixelType * s= &signal(0, 0);
		for(size_t i=0; i<this->clusters[cl].size()-1; ++i)
		{
			*s++ = -(*p++ - *q++);
		}
		//where does the radius get smaller
		std::vector<Center2D> blobs;
		finder.get_centers(signal, blobs);
		if(blobs.empty())
			continue;

		//split, from end to begin
		for(size_t i=0; i<blobs.size(); ++i)
		{
			this->clusters.push_back(Cluster());
			Cluster::iterator u = this->clusters[cl].begin();
			std::advance(u, (size_t)blobs[blobs.size()-i-1][0]);
			this->clusters.back().splice(this->clusters.back().begin(), this->clusters[cl], u, this->clusters[cl].end());
		}
	}*/
}

void Reconstructor::get_blobs(std::deque<Center3D>& centers)
{
	//const size_t margin = 6;
	centers.clear();
	boost::ptr_map<size_t, MultiscaleFinder1D> finders;
	for(std::deque<Cluster>::const_iterator cl=this->clusters.begin(); cl!=this->clusters.end(); ++cl)
	{
		if(cl->size()*3<6)
			continue;
		const size_t margin = cl->size();
		//fetch the needed 1D finder
		boost::ptr_map<size_t, MultiscaleFinder1D>::iterator f_it = finders.find(cl->size());
		if(f_it==finders.end())
		{
			//create the finder and cache it
			std::auto_ptr<MultiscaleFinder1D> finder(new MultiscaleFinder1D(cl->size()+2*margin));
			f_it = finders.insert(cl->size(), finder).first;
		}
		MultiscaleFinder1D &finder = *f_it->second;

		//copy the radii adding margins on each size to allow blob tracking on short signals
		OctaveFinder::Image signal(1, cl->size()+2*margin);
		signal.setTo(0.9*cl->back().r);
		std::fill_n(signal.begin(), margin, 0.9*cl->front().r);
		Cluster::const_iterator c=cl->begin();
		OctaveFinder::PixelType * s= &signal(0, margin);
		for(size_t i=0; i<cl->size(); ++i)
		{
			*s++ = c->r;
			c++;
		}
		std::vector<Center1D> blobs;
		finder.get_centers(signal, blobs);
		removeOverlapping_brute_force(blobs);

		for(std::vector<Center1D>::const_iterator b=blobs.begin(); b!=blobs.end(); ++b)
		{
			//get in the cluster just before the blob
			const size_t pos = (*b)[0];
			if(pos<margin || pos>cl->size()+margin)
				continue;
			const double frac = (*b)[0] - pos;
			Cluster::const_iterator it = cl->begin();
			std::advance(it, pos-margin);
			//add the center
			centers.push_back(*it);
			//centers.back()[2] = (*b)[0] - margin;
			//modulate by the next position (that may be closer to the blob position)
			it++;
			if(it!=cl->end())
			{
				centers.back()[0] += frac * ((*it)[0] - centers.back()[0]);
				centers.back()[1] += frac * ((*it)[1] - centers.back()[1]);
				centers.back()[2] += frac * ((*it)[2] - centers.back()[2]) - 0.5;
				centers.back().r += frac * (it->r - centers.back().r);
				centers.back().intensity += frac * (it->intensity - centers.back().intensity);
			}
		}
		/*if(blobs.empty())
		{
			//the signal is probably too short to localize a blob. We just take the maximum of the signal.
			centers.push_back(*std::max_element(cl->begin(), cl->end(), compare_radii<3>()));
		}*/
	}
}

void Reconstructor::links_by_brute_force(const Frame& fr, std::vector<double> &distances, std::vector<size_t> &from, std::vector<size_t> &to, const double &tolerance) const
{
	const size_t n = fr.size() * this->last_frame.size();
	distances.reserve(n);
	from.reserve(n);
	to.reserve(n);
	std::back_insert_iterator<std::vector<double> > d(distances);
	std::back_insert_iterator<std::vector<size_t> > i(from), j(to);
	for(size_t f=0; f<this->last_frame.size(); ++f)
		for(size_t t=0; t<fr.size(); ++t)
		{
			const double dist = pow(this->last_frame[f][0] - fr[t][0], 2) + pow(this->last_frame[f][1] - fr[t][1], 2);
			if(dist < pow((this->last_frame[f].r + fr[t].r) * tolerance, 2))
			{
				*i++ = f;
				*j++ = t;
				*d++ = dist;
			}
		}
}


/**
 * \param tolerance Fraction of the contact distance (sum of radii) accepted. For tolerance<=1 accept overlap only.
 */
void Reconstructor::links_by_RStarTree(const Frame& fr, const RTree& tree, std::vector<double> &distances, std::vector<size_t> &from, std::vector<size_t> &to, const double &max_dist) const
{
	//(over)reserve memory
	const size_t n = 12 * max(fr.size(), this->last_frame.size());
	distances.clear();
	distances.reserve(n);
	from.clear();
	from.reserve(n);
	to.clear();
	to.reserve(n);

	//spatial index the new frame
	/*RTree tree;
	for(size_t p=0; p<fr.size(); ++p)
		tree.Insert(p, get_bb(fr[p]));*/

	//for each particle in previous frame, get all the particles in new frame that have an overlap with it
	const double max_distsq = max_dist * max_dist;
	for(size_t p=0; p<this->last_frame.size(); ++p)
	{
		std::list<size_t> ngb1;
		tree.Query(
				RTree::AcceptOverlapping(get_bb_margin(this->last_frame[p], max_dist)),
				Gatherer<2>(ngb1)
		);
		for(std::list<size_t>::const_iterator it= ngb1.begin(); it!=ngb1.end(); ++it)
		{
			const double dist = pow(this->last_frame[p][0] - fr[*it][0], 2) + pow(this->last_frame[p][1] - fr[*it][1], 2);
			if(dist < max_distsq)
			{
				distances.push_back(dist);
				from.push_back(p);
				to.push_back(*it);
			}
		}
	}

}

/**
 * \param tolerance Fraction of the contact distance (sum of radii) accepted. For tolerance<=1 accept overlap only.
 */
/*void Reconstructor::links_by_kdTree(const Frame& fr, std::vector<double> &distances, std::vector<size_t> &from, std::vector<size_t> &to, const double &tolerance) const
{
	typedef KDTree::KDTree<2, Center2D> tree_type;
	//(over)reserve memory
	const size_t n = 12 * max(fr.size(), this->last_frame.size());
	distances.clear();
	distances.reserve(n);
	from.clear();
	from.reserve(n);
	to.clear();
	to.reserve(n);

	//spatial index the new frame
	tree_type tree;
	for(size_t p=0; p<fr.size(); ++p)
		tree.insert(p);

	//for each particle in previous frame, get all the particles in new frame that have an overlap with it
	for(size_t p=0; p<this->last_frame.size(); ++p)
	{
		std::list<size_t> ngb1;
		tree.Query(
				RTree::AcceptOverlapping(get_bb(this->last_frame[p], tolerance)),
				Gatherer(ngb1)
		);
		for(std::list<size_t>::const_iterator it= ngb1.begin(); it!=ngb1.end(); ++it)
		{
			const double dist = pow(this->last_frame[p][0] - fr[*it][0], 2) + pow(this->last_frame[p][1] - fr[*it][1], 2);
			if(dist < pow((this->last_frame[p].r + fr[*it].r) * tolerance, 2))
			{
				distances.push_back(dist);
				from.push_back(p);
				to.push_back(*it);
			}
		}
	}

}*/

}//namespace Colloids
