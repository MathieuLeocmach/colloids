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
		//brute force
		const size_t n = fr.size() * this->last_frame.size();
		std::vector<double> distances(n);
		std::vector<size_t> from(n), to(n);
		std::vector<double>::iterator d = distances.begin();
		std::vector<size_t>::iterator i = from.begin(), j = to.begin();
		for(size_t f=0; f<this->last_frame.size(); ++f)
			for(size_t t=0; t<fr.size(); ++t)
			{
				*i++ = f;
				*j++ = t;
				*d++ = pow(this->last_frame[f][0] - fr[t][0], 2) + pow(this->last_frame[f][1] - fr[t][1], 2);
			}
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

}//namespace Colloids
