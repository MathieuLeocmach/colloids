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

*/

#include "index.hpp"
using namespace std;
using namespace Colloids;

/** @brief Get the indices of the objects contained inside a reduction of the maximum bounding box  */
set<size_t> SpatialIndex::getInside(const double &margin) const
{
    BoundingBox insideBox = getOverallBox();
    for(size_t i=0;i<3;++i)
        if(insideBox.edges[i].second-insideBox.edges[i].first>2*margin)
        {
            insideBox.edges[i].first  += margin;
            insideBox.edges[i].second -= margin;
        }
    //gets the indexes of the particles totally contained inside this volume
    return (*this)(insideBox);
}

/** @brief Get objects spanning the whole interval inside the query box  */
set<size_t> SpatioTemporalIndex::operator()(const BoundingBox &b) const
{
    return (*this)(TimeBox(getOverallInterval(),b));
}

/** @brief Get all the objects spanning the query interval  */
set<size_t> SpatioTemporalIndex::operator()(const Interval &in) const
{
    return (*this)(TimeBox(in,getOverallBox()));
}

/** @brief Get the indices of the objects spanning the query interval inside a reduction of the maximum bounding box  */
set<size_t> SpatioTemporalIndex::getSpanningInside(const Interval &in,const double &margin) const
{
    BoundingBox insideBox = getOverallBox();
    for(size_t i=0;i<3;++i)
        if(insideBox.edges[i].second-insideBox.edges[i].first>2*margin)
        {
            insideBox.edges[i].first  += margin;
            insideBox.edges[i].second -= margin;
        }
    //gets the indexes of the particles totally contained inside this volume
    return (*this)(TimeBox(in,insideBox));
}

/** @brief Get the indices of the objects spanning the whole interval inside a reduction of the maximum bounding box
  */
set<size_t> SpatioTemporalIndex::getInside(const double &margin) const
{
    return getSpanningInside(getOverallInterval(),margin);
}



/** @brief insert  */
void BruteSpatialIndex::insert(const size_t &i, const BoundingBox &b)
{
    if(!overallBox)
        overallBox = new BoundingBox(b);
    else
        overallBox->stretch(b);
    items.insert(make_pair(i,b));
}

/** @brief get all items included in the query box
    complexity ~N
*/
set<size_t> BruteSpatialIndex::operator()(const BoundingBox &b) const
{
    set<size_t> ret;
    for(multimap<size_t, BoundingBox>::const_iterator it = items.begin(); it!= items.end();++it)
        if(b.encloses(it->second))
            ret.insert(ret.end(), it->first);
    return ret;
}

/** @brief Translate all bounding boxes  */
void BruteSpatialIndex::operator+=(const Coord &v)
{
    for(multimap<size_t, BoundingBox>::iterator it = items.begin(); it!= items.end();++it)
        it->second += v;
}

void RStarIndex_S::insert(const size_t &i, const BoundingBox &b)
{
    tree.Insert(i,b);
}

/** @brief Get the indices of the objects whose bounding boxes are contained inside the query box */
set<size_t> RStarIndex_S::operator()(const BoundingBox &b) const
{
    return tree.Query(RTree::AcceptEnclosing(b), Gatherer()).gathered;
}

/** @brief insertion  */
void TreeIndex_T::insert(const size_t &i, const Interval &in)
{
    if(!overallInterval)
        overallInterval = new Interval(in);
    else
    {
        overallInterval->first = min(overallInterval->first, in.first);
        overallInterval->second = max(overallInterval->second, in.second);
    }
    for(size_t t = tree.size();t<=in.first;++t)
        tree.push_back(new boost::ptr_vector< std::set<size_t> >);
    for(size_t t = tree[in.first].size();t<=in.second-in.first;++t)
        tree[in.first].push_back(new std::set<size_t>);
    tree[in.first][in.second-in.first].insert(tree[in.first][in.second-in.first].end(),i);
}

/** @brief Get the objects spanning at least the query interval */
set<size_t> TreeIndex_T::operator()(const Interval &in) const
{
    set<size_t> sel;
    for(size_t t0=0;t0<=min(in.first,tree.size());++t0)
        for(size_t t1=in.second-in.first;t1<tree[t0].size();++t1)
            copy(tree[t0][t1].begin(), tree[t0][t1].end(), inserter(sel,sel.end()));
    return sel;
}


