/*
 * reconstructor.cpp
 *
 *  Created on: 3 août 2011
 *      Author: mathieu
 */

#include "reconstructor.h"
#include "RStarTree/RStarTree.h"
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

void Reconstructor::push_back(const Frame &fr)
{
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
		//use RStarTree spatial indexing
		this->links_by_RStarTree(fr, distances, from, to);
		//brute force
		//this->links_by_brute_force(fr, distances, from, to);
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

void Reconstructor::links_by_brute_force(const Frame& fr, std::vector<double> &distances, std::vector<size_t> &from, std::vector<size_t> &to) const
{
	const size_t n = fr.size() * this->last_frame.size();
	distances.resize(n);
	from.resize(n);
	to.resize(n);
	std::vector<double>::iterator d = distances.begin();
	std::vector<size_t>::iterator i = from.begin(), j = to.begin();
	for(size_t f=0; f<this->last_frame.size(); ++f)
		for(size_t t=0; t<fr.size(); ++t)
		{
			*i++ = f;
			*j++ = t;
			*d++ = pow(this->last_frame[f][0] - fr[t][0], 2) + pow(this->last_frame[f][1] - fr[t][1], 2);
		}
}

inline RStarBoundingBox<2,double> get_bb(const Center2D &c)
{
	RStarBoundingBox<2,double> bb;
	bb.edges[0].first = c[0] - c.r;
	bb.edges[1].first = c[1] - c.r;
	bb.edges[0].second = c[0] + c.r;
	bb.edges[1].second = c[1] + c.r;
	return bb;
}
typedef RStarTree<size_t, 2, 4, 32, double> 	RTree;

struct Gatherer {
	std::list<size_t> *gathered;
	bool ContinueVisiting;

	Gatherer(std::list<size_t> &result) : gathered(&result), ContinueVisiting(true) {};

	void operator()(const RTree::Leaf * const leaf)
	{
		gathered->push_back(leaf->leaf);
	}
};

void Reconstructor::links_by_RStarTree(const Frame& fr, std::vector<double> &distances, std::vector<size_t> &from, std::vector<size_t> &to) const
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
	RTree tree;
	for(size_t p=0; p<fr.size(); ++p)
		tree.Insert(p, get_bb(fr[p]));

	//for each particle in previous frame, get all the particles in new frame that have an overlap with it
	for(size_t p=0; p<this->last_frame.size(); ++p)
	{
		std::list<size_t> ngb1;
		tree.Query(
				RTree::AcceptEnclosing(get_bb(this->last_frame[p])),
				Gatherer(ngb1)
		);
		for(std::list<size_t>::const_iterator it= ngb1.begin(); it!=ngb1.end(); ++it)
		{
			const double dist = pow(this->last_frame[p][0] - fr[*it][0], 2) + pow(this->last_frame[p][1] - fr[*it][1], 2);
			if(dist < pow(this->last_frame[p].r + fr[*it].r, 2))
			{
				distances.push_back(dist);
				from.push_back(*it);
				to.push_back(p);
			}
		}
	}

}

}//namespace Colloids
