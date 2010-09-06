/*
 *  Copyright (c) 2008 Dustin Spicuzza <dustin@virtualroadside.com>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of version 2 of the GNU General Public License
 *  as published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */


#ifndef RStarBoundingBox_H
#define RStarBoundingBox_H

#include <limits>
#include <utility>
#include <cstddef>
#include <string>
#include <sstream>

/**
    \struct RstarBoundingBox RstarBoundingBox.h
    \brief A N-dim rectangle bounding objects
    \tparam dimensions Dimension of the working space
    \tparam grain Type of the bounding coordinates (int,float...)
*/
template <std::size_t dimensions,typename grain>
struct RStarBoundingBox {

	typedef grain grain_type;

	// edges[x].first is low value, edges[x].second is high value
	std::pair<grain, grain> edges[dimensions];

	// forces all edges to their extremes so we can stretch() it
	void reset()
	{
		for (std::size_t axis = 0; axis < dimensions; axis++)
		{
			edges[axis].first = std::numeric_limits<grain>::max();
			edges[axis].second = std::numeric_limits<grain>::min();
		}
	}

	// returns a new bounding box that has the maximum boundaries
	static RStarBoundingBox MaximumBounds()
	{
		RStarBoundingBox<dimensions,grain> bound;
		bound.reset();
		return bound;
	}


	// fits another box inside of this box, returns true if a stretch occured
	bool stretch(const RStarBoundingBox<dimensions,grain> &bb)
	{
		bool ret = false;

		for (std::size_t axis = 0; axis < dimensions; axis++)
		{

			if (edges[axis].first > bb.edges[axis].first)
			{
				edges[axis].first = bb.edges[axis].first;
				ret = true;
			}

			if (edges[axis].second < bb.edges[axis].second)
			{
				edges[axis].second = bb.edges[axis].second;
				ret = true;
			}
		}

		return ret;
	}

	// the sum of all deltas between edges
	inline grain edgeDeltas() const
	{
		grain distance = 0;
		for (std::size_t axis = 0; axis < dimensions; axis++)
			distance += edges[axis].second - edges[axis].first;

		return distance;
	}

	// calculates the area of a bounding box
	inline double area() const
	{
		double area = 1;
		for (std::size_t axis = 0; axis < dimensions; axis++)
			area *= (double)(edges[axis].second - edges[axis].first);

		return area;
	}

	// this determines if a bounding box is fully contained within this bounding box
	inline bool encloses(const RStarBoundingBox<dimensions,grain>& bb) const
	{
		// if (y1 < x1 || x2 < y2)
		for (std::size_t axis = 0; axis < dimensions; axis++)
			if (bb.edges[axis].first < edges[axis].first || edges[axis].second < bb.edges[axis].second)
				return false;

		return true;
	}

	// a quicker way to determine if two bounding boxes overlap
	inline bool overlaps(const RStarBoundingBox<dimensions,grain>& bb) const
	{
		// do it this way so theres no equal signs (in case of doubles)
		// if (!(x1 < y2) && !(x2 > y1))
		for (std::size_t axis = 0; axis < dimensions; axis++)
		{
			if (!(edges[axis].first < bb.edges[axis].second) || !(bb.edges[axis].first < edges[axis].second))
				return false;
		}

		return true;
	}

	// calculates the total overlapping area of two boxes
	double overlap(const RStarBoundingBox<dimensions,grain>& bb) const
	{
		double area = 1.0;
		for (std::size_t axis = 0; area && axis < dimensions; axis++)
		{
			// this makes it easier to understand
			const grain x1 = edges[axis].first;
			const grain x2 = edges[axis].second;
			const grain y1 = bb.edges[axis].first;
			const grain y2 = bb.edges[axis].second;

			// left edge outside left edge
			if (x1 < y1)
			{
				// and right edge inside left edge
				if (y1 < x2)
				{
					// right edge outside right edge
					if (y2 < x2)
						area *= (double)( y2 - y1 );
					else
						area *= (double)( x2 - y1 );

					continue;
				}
			}
			// right edge inside left edge
			else if (x1 < y2)
			{
				// right edge outside right edge
				if (x2 < y2)
					area *= (double)( x2 - x1 );
				else
					area *= (double)( y2 - x1 );

				continue;
			}

			// if we get here, there is no overlap
			return 0.0;
		}

		return area;
	}

	// sums the total distances from the center of another bounding box
	double distanceFromCenter(const RStarBoundingBox<dimensions,grain>& bb) const
	{
		double distance = 0, t;
		for (std::size_t axis = 0; axis < dimensions; axis++)
		{
			t = ((double)edges[axis].first + (double)edges[axis].second +
			     (double)bb.edges[axis].first + (double)bb.edges[axis].second)
				 /2.0;
			distance += t*t;
		}

		return distance;
	}

	// determines if two bounding boxes are identical
	bool operator==(const RStarBoundingBox<dimensions,grain>& bb)
	{
		for (std::size_t axis = 0; axis < dimensions; axis++)
			if (edges[axis].first != bb.edges[axis].first || edges[axis].second != bb.edges[axis].second)
				return false;

		return true;
	}

	/** \brief intersection assignement operator */
	RStarBoundingBox<dimensions,grain>& operator&=(const RStarBoundingBox<dimensions,grain> &bb)
	{
	    for (std::size_t axis = 0; axis < dimensions; axis++)
	    {
	        edges[axis].first = max(bb.edges[axis].first,edges[axis].first);
	        edges[axis].second = min(bb.edges[axis].second,edges[axis].second);
	    }
	    return *this;
	}

    /** \brief intersection operator */
	RStarBoundingBox<dimensions,grain> operator&&(RStarBoundingBox<dimensions,grain> bb)
	{
	    return bb&=(*this);
	}

	// very slow, use for debugging only
	std::string ToString() const
	{
		std::stringstream name("");
		name << "[";
		for (std::size_t axis = 0; axis < dimensions; axis++)
		{
			name << "(" << edges[axis].first << "," << edges[axis].second << ")";
			if (axis != dimensions -1)
				name << ",";
		}
		name << "]";

		return name.str();
	}

	/** \brief translation of the box  */
	template<class Coordinates>
	RStarBoundingBox<dimensions,grain>& operator+=(const Coordinates &v)
	{
		for (std::size_t axis = 0; axis < dimensions; axis++)
		{
			edges[axis].first+=v[axis];
			edges[axis].second+=v[axis];
		}
		return *this;
	}
};



template <std::size_t dimensions,typename grain>
struct RStarBoundedItem {
	typedef RStarBoundingBox<dimensions,grain> BoundingBox;

	BoundingBox bound;
};


/**********************************************************
 * Functor used to iterate over a set and stretch a
 * bounding box
 **********************************************************/

// for_each(items.begin(), items.end(), StretchBoundedItem::BoundingBox(bound));
template <typename BoundedItem>
struct StretchBoundingBox :
	public std::unary_function< const BoundedItem * const, void >
{
	typename BoundedItem::BoundingBox * m_bound;
	explicit StretchBoundingBox(typename BoundedItem::BoundingBox * bound) : m_bound(bound) {}

	void operator() (const BoundedItem * const item)
	{
		m_bound->stretch(item->bound);
	}
};


/**********************************************************
 * R* Tree related functors used for sorting BoundedItems
 *
 * TODO: Take advantage of type traits
 **********************************************************/

template <typename BoundedItem>
struct SortBoundedItemsByFirstEdge :
	public std::binary_function< const BoundedItem * const, const BoundedItem * const, bool >
{
	const std::size_t m_axis;
	explicit SortBoundedItemsByFirstEdge (const std::size_t axis) : m_axis(axis) {}

	bool operator() (const BoundedItem * const bi1, const BoundedItem * const bi2) const
	{
		return bi1->bound.edges[m_axis].first < bi2->bound.edges[m_axis].first;
	}
};

template <typename BoundedItem>
struct SortBoundedItemsBySecondEdge :
	public std::binary_function< const BoundedItem * const, const BoundedItem * const, bool >
{
	const std::size_t m_axis;
	explicit SortBoundedItemsBySecondEdge (const std::size_t axis) : m_axis(axis) {}

	bool operator() (const BoundedItem * const bi1, const BoundedItem * const bi2) const
	{
		return bi1->bound.edges[m_axis].second < bi2->bound.edges[m_axis].second;
	}
};


template <typename BoundedItem>
struct SortBoundedItemsByDistanceFromCenter :
	public std::binary_function< const BoundedItem * const, const BoundedItem * const, bool >
{
	const typename BoundedItem::BoundingBox * const m_center;
	explicit SortBoundedItemsByDistanceFromCenter(const typename BoundedItem::BoundingBox * const center) : m_center(center) {}

	bool operator() (const BoundedItem * const bi1, const BoundedItem * const bi2) const
	{
		return bi1->bound.distanceFromCenter(*m_center) < bi2->bound.distanceFromCenter(*m_center);
	}
};

template <typename BoundedItem>
struct SortBoundedItemsByAreaEnlargement :
	public std::binary_function< const BoundedItem * const, const BoundedItem * const, bool >
{
	const double area;
	explicit SortBoundedItemsByAreaEnlargement(const typename BoundedItem::BoundingBox * center) : area(center->area()) {}

	bool operator() (const BoundedItem * const bi1, const BoundedItem * const bi2) const
	{
		return area - bi1->bound.area() < area - bi2->bound.area();
	}
};

template <typename BoundedItem>
struct SortBoundedItemsByOverlapEnlargement :
	public std::binary_function< const BoundedItem * const, const BoundedItem * const, bool >
{
	const typename BoundedItem::BoundingBox * const m_center;
	explicit SortBoundedItemsByOverlapEnlargement(const typename BoundedItem::BoundingBox * const center) : m_center(center) {}

	bool operator() (const BoundedItem * const bi1, const BoundedItem * const bi2) const
	{
		return bi1->bound.overlap(*m_center) < bi2->bound.overlap(*m_center);
	}
};

/**********************************************************
 * Functor used to translate the bounding boxes of a set
 * of BoundedItems by the same vector
 **********************************************************/

template <typename BoundedItem, typename Coordinates>
struct TranslateBoundingBox :
	public std::unary_function< BoundedItem *, void >
{
	const Coordinates * trans;
	const bool ContinueVisiting;
	explicit TranslateBoundingBox(const Coordinates * v) : trans(v), ContinueVisiting(true) {}

	void operator() (BoundedItem * item)
	{
		item->bound+=(*trans);
	}
};


#endif
