/**
    Copyright 2008,2009 Mathieu Leocmach

    This file is part of Colloids.

    Colloids is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Colloids is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Colloids.  If not, see <http://www.gnu.org/licenses/>.

 * \file particles.hpp
 * \brief Defines classes for spatial indexing
 * \author Mathieu Leocmach
 * \date 21 Novemeber 2009
 *
 */

#ifndef index_H
#define index_H

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <tvmet/Vector.h>
#include <boost/ptr_container/ptr_container.hpp>

#include "RStarTree/RStarTree.h"

namespace Colloids
{
    typedef tvmet::Vector<double, 3>            Coord;
    typedef RStarBoundingBox<3,double>          BoundingBox;
    typedef std::pair<size_t,size_t>            Interval;
    typedef std::pair<Interval, BoundingBox>    TimeBox;

    /** \brief A template class defining the common interface of index classes   */
    template<class BoundingItem>
    struct BasicIndex
    {
        virtual ~BasicIndex(){return;};
        virtual void insert(const size_t &i, const BoundingItem &b){};
        explicit BasicIndex(const std::vector<BoundingItem> &items)
        {
            for(size_t i=0;i<items.size();++i)
                this->insert(i,items[i]);
            return;
        };
        virtual std::set<size_t> operator()(const BoundingItem &b) const = 0;
    };

    /** \brief A virtual class defining the common interface of spatial index classes   */
    class SpatialIndex : public BasicIndex<BoundingBox>
    {
        public:
            explicit SpatialIndex(const std::vector<BoundingBox> &items) : BasicIndex<BoundingBox>(items){return;};
            virtual std::set<size_t> getInside(const double &margin) const;
            /** @brief Translate index */
            virtual void operator+=(const Coord &v) = 0;
            virtual BoundingBox getOverallBox() const = 0;
    };

    /** \brief A virtual class defining the common interface of temporal index classes   */
    class TemporalIndex : public BasicIndex<Interval>
    {
        public:
            explicit TemporalIndex(const std::vector<Interval> &items) : BasicIndex<Interval>(items){return;};
            virtual Interval getOverallInterval() const = 0;
    };

    /** \brief A virtual class defining the common interface of spatio-temporal index classes   */
    class SpatioTemporalIndex : public BasicIndex<TimeBox>
    {
        public:
            explicit SpatioTemporalIndex(const std::vector<TimeBox> &items) : BasicIndex<TimeBox>(items){return;};
            virtual std::set<size_t> operator()(const TimeBox &b) const = 0;
            virtual std::set<size_t> operator()(const BoundingBox &b) const = 0;
            virtual std::set<size_t> operator()(const Interval &in) const = 0;
            virtual std::set<size_t> getSpanningInside(const Interval &in,const double &margin) const;
            virtual std::set<size_t> getInside(const double &margin) const;
            virtual BoundingBox getOverallBox() const = 0;
            virtual Interval getOverallInterval() const = 0;
    };

    /** \brief Brute force implementation of spatial index. Utterly inefficient. */
    class BruteSpatialIndex : public SpatialIndex
    {
        BoundingBox *overallBox;
        std::multimap<size_t, BoundingBox> items;

        public:
            BruteSpatialIndex(const std::vector<BoundingBox> &items) : SpatialIndex(items){return;};
            void insert(const size_t &i, const BoundingBox &b);
            std::set<size_t> operator()(const BoundingBox &b) const;
            void operator+=(const Coord &v);
            BoundingBox getOverallBox() const{return *overallBox;};
    };

    /** \brief R*Tree implementation of spatial index. */
    class RStarIndex_S : public SpatialIndex
    {
        public:
            typedef RStarTree<size_t, 3, 2, 16,double> 	RTree;
            RTree tree;

        /** \brief Visitor gathering particles indices */
        struct Gatherer {
            std::set<size_t> gathered;
            bool ContinueVisiting;

            Gatherer() : gathered(), ContinueVisiting(true) {};

            void operator()(const RTree::Leaf * const leaf)
            {
                gathered.insert(leaf->leaf);
            }
        };
            RStarIndex_S(const std::vector<BoundingBox> &items) : SpatialIndex(items){return;};
            void insert(const size_t &i, const BoundingBox &b);
            std::set<size_t> operator()(const BoundingBox &b) const;
            void operator+=(const Coord &v) {tree+=v;};
            BoundingBox getOverallBox() const{return tree.getOverallBox();};
    };

    /** \brief simple tree implementation of temporal index */
    class TreeIndex_T : public TemporalIndex
    {
        /** \brief  tree[start][end] is the set of objects spanning exactly the interval [start,end] */
        boost::ptr_vector< boost::ptr_vector< std::set<size_t> > > tree;
        Interval *overallInterval;

        public:
            TreeIndex_T(const std::vector<Interval> &items) : TemporalIndex(items){return;};
            void insert(const size_t &i, const Interval &in);
            std::set<size_t> operator()(const Interval &in) const;
            Interval getOverallInterval() const{return *overallInterval;};
    };

    /** \brief A slice of a SpatioTemporalIndex at time t is a SpatialIndex */
    template<class ParentSTIndex, class Converter>
    class SpatioTemporalIndex_slice : public SpatialIndex
    {
        const ParentSTIndex *parent;
        const size_t time;
        Converter &conv;

        public:
            SpatioTemporalIndex_slice(ParentSTIndex &par, const size_t &t, Converter &converter)
             : parent(&par), time(t), conv(converter) {return;};

            std::set<size_t> operator()(const BoundingBox &b) const
            {
                std::set<size_t> ret, trs = (*parent)(std::make_pair(std::make_pair(time,time),b));
                std::transform(trs.begin(), trs.end(), std::inserter(ret,ret.end()), this->conv);
                return ret;
            }
            BoundingBox getOverallBox() const{return parent->getOverallBox();};
            /** \brief A slice is read only */
            void insert(const size_t &i, const BoundingBox &b){return;};
            void operator+=(const Coord &v) const{return;};

    };


};
#endif
