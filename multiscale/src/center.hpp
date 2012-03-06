/*
 * Center2D.h
 *
 *  Created on: 3 août 2011
 *      Author: mathieu
 */

#ifndef CENTER_H_
#define CENTER_H_

#include "../RStarTree/RStarTree.h"
#include <boost/array.hpp>

namespace Colloids {

	struct Center_base
	{
		double r, intensity;
		explicit Center_base(const double &r=0, const double &i=0) : r(r), intensity(i){};
		virtual const double& operator[](size_t i) const=0;
		virtual double& operator[](size_t i)=0;
	};

	template<int D>
	struct Center : public Center_base
	{
		boost::array<double, D> coords;
		static const int dimension = D;

		explicit Center(const double &v=0.0, const double &r=0, const double &i=0) : Center_base(r, i)
		{
			std::fill(coords.begin(), coords.end(), v);
		};
		explicit Center(const Center<D-1> &c, const double &additional_coord): Center_base(c.r, c.intensity)
		{
			std::copy(c.coords.begin(), c.coords.end(), coords.begin());
			coords.back() = additional_coord;
		};
		virtual const double& operator[](size_t i) const {return coords[i];};
		virtual double& operator[](size_t i){return coords[i];};
		inline double operator-(Center<D> other) const
		{
			double d = 0;
			for(size_t i=0; i<coords.size(); ++i)
				d += pow(coords[i]-other[i], 2);
			return d;
		}
	};
	typedef Center<1> Center1D;
	typedef Center<2> Center2D;
	typedef Center<3> Center3D;

	template<int D>
	struct compare_radii : std::binary_function<bool, const Center<D>&, const Center<D>& >
	{
		bool operator()(const Center<D>& a, const Center<D> &b) const {return a.r < b.r;}
	};
	template<int D>
	struct compare_intensities : std::binary_function<bool, const Center<D>&, const Center<D>& >
	{
		bool operator()(const Center<D>& a, const Center<D> &b) const {return a.intensity < b.intensity;}
	};
	template<int d>
	struct compare_coord : std::binary_function<bool, const Center_base&, const Center_base& >
	{
		bool operator()(const Center_base& a, const Center_base &b) const {return a[d] < b[d];}
	};

	template<int D>
	inline RStarBoundingBox<D,double> get_bb(const Center<D> &c, const double &tolerance=1.0)
	{
		RStarBoundingBox<D,double> bb;
		for(size_t d=0; d<D; ++d)
		{
			bb.edges[d].first = c[d] - c.r * tolerance;
			bb.edges[d].second = c[d] + c.r * tolerance;
		}
		return bb;
	}
	template<int D>
	inline RStarBoundingBox<D,double> get_bb_margin(const Center<D> &c, const double &margin=0.0)
	{
		RStarBoundingBox<D,double> bb;
		for(size_t d=0; d<D; ++d)
		{
			bb.edges[d].first = c[d] - margin;
			bb.edges[d].second = c[d] + margin;
		}
		return bb;
	}

	template<int D, std::size_t min_child_items=4, std::size_t max_child_items=32>
	struct Gatherer {
		typedef RStarTree<size_t, D, min_child_items, max_child_items, double> RTree;
		std::list<size_t> *gathered;
		bool ContinueVisiting;

		Gatherer(std::list<size_t> &result) : gathered(&result), ContinueVisiting(true) {};

		void operator()(const typename RTree::Leaf * const leaf)
		{
			gathered->push_back(leaf->leaf);
		}
	};

	template<int D>
	std::auto_ptr< RStarTree<size_t, D, 4, 32, double> > removeOverlapping(std::vector<Center<D> > &centers, const double& tolerance=0.5)
	{
		typedef RStarTree<size_t, D, 4, 32, double> RTree;
		typedef std::vector<Center<D> > Centers;
		std::auto_ptr<RTree> tree(new RTree());
		Centers filtered;
		filtered.reserve(centers.size());
		//sort the centers by decreasing response (increasing negative intensity)
		std::sort(centers.begin(), centers.end(), compare_intensities<D>());
		//insert the centers one by one, if no overlap
		const double tolsq = tolerance * tolerance;
		for(size_t p=0; p<centers.size(); ++p)
		{
			//norm 1 overlapping
			typename RTree::BoundingBox bb = get_bb(centers[p], tolerance);
			std::list<size_t> overlapping;
			tree->Query(typename RTree::AcceptOverlapping(bb), Gatherer<D>(overlapping));
			bool is_overlapping = false;
			//norm 2 overlapping
			for(std::list<size_t>::const_iterator q= overlapping.begin(); q!=overlapping.end(); ++q)
				if(tolsq*(centers[p]-filtered[*q]) < pow(centers[p].r + filtered[*q].r ,2))
				{
					is_overlapping = true;
					break;
				}
			if(!is_overlapping)
			{
				tree->Insert(filtered.size(), bb);
				filtered.push_back(centers[p]);
			}
		}
		centers.swap(filtered);
		return tree;
	}
	template<int D>
	void removeOverlapping_brute_force(std::vector<Center<D> > &centers)
	{
		typedef std::vector<Center<D> > Centers;
		Centers filtered;
		filtered.reserve(centers.size());
		//sort the centers by decreasing response (increasing negative intensity)
		std::sort(centers.begin(), centers.end(), compare_intensities<D>());
		//insert the centers one by one, if no overlap
		for(size_t p=0; p<centers.size(); ++p)
		{
			bool is_overlapping = false;
			for(size_t q=0; q<filtered.size(); ++q)
				if(centers[p]-filtered[q] < pow(centers[p].r + filtered[q].r ,2))
				{
					is_overlapping = true;
					break;
				}
			if(!is_overlapping)
				filtered.push_back(centers[p]);
		}
		centers.swap(filtered);
	}
	template<int D>
	std::auto_ptr< RStarTree<size_t, D, 4, 32, double> > removeHalfOverlapping(std::vector<Center<D> > &centers)
	{
		typedef RStarTree<size_t, D, 4, 32, double> RTree;
		typedef std::vector<Center<D> > Centers;
		std::auto_ptr<RTree> tree(new RTree());
		Centers filtered;
		filtered.reserve(centers.size());
		//sort the centers by decreasing response (increasing negative intensity)
		std::sort(centers.begin(), centers.end(), compare_intensities<D>());
		//insert the centers one by one, if no overlap
		for(size_t p=0; p<centers.size(); ++p)
		{
			//norm 1 overlapping
			typename RTree::BoundingBox bb = get_bb(centers[p], 1.0);
			std::list<size_t> overlapping;
			tree->Query(typename RTree::AcceptOverlapping(bb), Gatherer<D>(overlapping));
			bool is_overlapping = false;
			//norm 2 overlapping
			for(std::list<size_t>::const_iterator q= overlapping.begin(); q!=overlapping.end(); ++q)
			{
				double distsq = centers[p]-filtered[*q];
				if(distsq < pow(std::max(centers[p].r, filtered[*q].r), 2))
				{
					is_overlapping = true;
					break;
				}
			}
			if(!is_overlapping)
			{
				tree->Insert(filtered.size(), bb);
				filtered.push_back(centers[p]);
			}
		}
		centers.swap(filtered);
		return tree;
	}

}

#endif /* CENTER_H_ */
