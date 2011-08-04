/*
 * reconstructor.h
 *
 *  Created on: 3 août 2011
 *      Author: mathieu
 */

#ifndef RECONSTRUCTOR_H_
#define RECONSTRUCTOR_H_

#include "center.h"
#include "traj.hpp"
#include <list>
#include <memory>

namespace Colloids
{

	class Reconstructor :boost::noncopyable
	{
	public:
		typedef std::list<Center3D> Cluster;
		typedef std::vector<Center2D> Frame;

		Reconstructor();
		virtual ~Reconstructor();

		//accessors
		inline bool empty() const {return !trajectories.get();}
		inline const size_t size() const {return empty()?0:trajectories->nbFrames();}
		inline const size_t nb_cluster() const {return clusters.size();}
		inline const std::deque<Cluster>& get_clusters() const {return clusters;}
		inline const TrajIndex& get_trajectories() const {assert(!empty()); return *trajectories;}

		//processing
		void clear();
		void push_back(const Frame &fr);
		void split_clusters();
		void get_blobs(std::list<Center3D>& blobs);

	private:
		std::auto_ptr<TrajIndex> trajectories;
		std::deque<Cluster> clusters;
		Frame last_frame;

		void links_by_brute_force(const Frame& fr, std::vector<double> &distances, std::vector<size_t> &from, std::vector<size_t> &to) const;
		void links_by_RStarTree(const Frame& fr, std::vector<double> &distances, std::vector<size_t> &from, std::vector<size_t> &to) const;

	};

}

#endif /* RECONSTRUCTOR_H_ */
