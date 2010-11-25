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
**/

#include "traj.hpp"
#include <boost/tokenizer.hpp>

using namespace std;
using namespace Colloids;


typedef multimap<double, size_t>   FolMap;

/** @brief Constructor from the number of particles in the first frame

Initial trajectories have the same indicies as initial positions
*/
TrajMap::TrajMap(const size_t &firstFrameSize)
{
    bm.assign(1,Frame());
    for(size_t i=0;i<firstFrameSize;++i)
        bm.back().insert(Link(i,i));
    nbTraj = firstFrameSize;
}

/** @brief Create a new frame by attributing (if possible) a follower to all the trajectories existing in previous frame.
Create new trajectories in the new frame if needed.
  * \param followersByDist The possible followers of each position in the last frame, sorted by the distance (last posistion, new position)
  * \param frameSize The total number of positions in the new frame
  * Any distence in the topological meaning is ok : d(a,b)=0 <=> a=b and for all a!=b d(a,b)>0.
  */
void Colloids::TrajMap::push_back(const vector< multimap<double, size_t> > &followersByDist, const size_t &frameSize)
{
    //convert the input into a list of links (pos,tr) sorted by distence between pos and the last position of the trajectory
    multimap<double, Link> potential_links;
    for(Frame::map_by<Coord>::const_iterator p_tr=bm.back().by<Coord>().begin();p_tr!=bm.back().by<Coord>().end();++p_tr)
        for(FolMap::const_iterator dist_fol=followersByDist[p_tr->first].begin();dist_fol!=followersByDist[p_tr->first].end();++dist_fol)
            potential_links.insert(make_pair(dist_fol->first, Link(dist_fol->second, p_tr->second)));

    //create the new frame
    bm.push_back(Frame());

    //insert the links into the new frame, statring by the shortest link
    //The links pointing to the same trajectory or to the same object as an already inserted link are automatically ignored
    for(multimap<double, Link>::const_iterator l=potential_links.begin();l!=potential_links.end();++l)
        bm.back().insert(l->second);

    //some elements of the new frame may not have found a trajectory to be linked to. Creates a new trajectory for each of them.
    vector<bool> linked(frameSize, false);
    for(Frame::map_by<Coord>::const_iterator p=bm.back().by<Coord>().begin(); p!=bm.back().by<Coord>().end(); ++p)
        linked[p->first]=true;
    for(size_t p =0;p<linked.size();++p)
        if(!linked[p])
            bm.back().insert(Link(p, nbTraj++));
}

/** @brief get the size of each Frame  */
vector<size_t> TrajMap::getFrameSizes() const
{
	vector<size_t> frameSizes(bm.size());
    std::transform(
		bm.begin(), bm.end(),
		frameSizes.begin(), mem_fun_ref(&Frame::size)
		);
	return frameSizes;
}




/** \brief feed an initialized trajectory with the position indices contained by an input stream */
istream& Colloids::operator>> (istream& is, Traj& tr )
{
    size_t index;
    while(!is.eof())
    {
        is >> index;
        tr.push_back(index);
    }
    return is;
}

/** \brief output a trajectory to an output stream */
ostream& Colloids::operator<< (ostream& os, const Traj& tr )
{
    os << tr.start_time << "\n";
    os << tr.steps[0];
    for(size_t i=1;i<tr.steps.size();++i)
        os << "\t" << tr.steps[i];
    return os;
}

/** @brief Get the subtrajectory starting at t0 and ending at t1  */
Traj Traj::subtraj(const size_t &t0, const size_t &t1) const
{
    Traj tr(t0);
    for(size_t t=t0; t<=t1; ++t)
        tr.push_back((*this)[t]);
    return tr;
}



/** \brief explains the trajectory error */
const char* TrajError::what() const throw()
{
    ostringstream os(ios_base::in | ios_base::out);
    os<< "This trajectory doesn't exist at time step "<<asked<<". It spans ["<<start<<","<<stop<<")";
    return os.str().c_str();
}
/** \brief explains the trajectory error */
const char* IdTrajError::what() const throw()
{
    ostringstream os(ios_base::in | ios_base::out);
    os<< "Trajectory "<<id<<" doesn't exist at time step "<<asked<<". It spans ["<<start<<","<<stop<<")";
    return os.str().c_str();
}

/** @brief Constructor from a TrajMap  */
TrajIndex::TrajIndex(const TrajMap &tm)
{
    typedef TrajMap::Frame::map_by<Traj>::const_iterator TrajP_it;
    for(size_t t=0; t<tm.size();++t)
    {
        //Delimitate already existing trajectories (tr<this.size()) for the trajectories that we have to make start now
        TrajP_it last_tr;
        if(this->empty())
            last_tr = tm[t].by<Traj>().begin();
        else
            last_tr = tm[t].by<Traj>().upper_bound(this->size()-1);
        //prolongating existing trajectories
        for(TrajP_it tr_p = tm[t].by<Traj>().begin(); tr_p!=last_tr; ++tr_p)
            (*this)[tr_p->first].push_back(tr_p->second);
        //starting new trajectories
        for(TrajP_it tr_p = last_tr; tr_p!=tm[t].by<Traj>().end(); ++tr_p)
            this->push_back(Traj(t, tr_p->second));
    }
    makeInverse(tm.getFrameSizes());
}



/** @brief get the latest time spaned by a trajectory  */
size_t TrajIndex::getMaxTime() const
{
    size_t max=0;
    for(const_iterator tr=begin();tr!=end();++tr)
        if(max<(*tr).last_time())
            max = (*tr).last_time();
    return max;
}

/** @brief longest_span return the number of time steps that spans the longest trajectory  */
size_t TrajIndex::longest_span() const
{
    vector<size_t> lengths(size());
    transform(
        this->begin(), this->end(),
        lengths.begin(),
        mem_fun_ref(&Traj::size)
        );
    return *max_element(lengths.begin(), lengths.end());
}



/** @brief count the number of trajectories in each time step  */
vector<size_t> TrajIndex::getFrameSizes(const size_t &length) const
{
	vector<size_t> frameSizes(length?length:(getMaxTime()+1), 0);
	for(const_iterator tr=begin(); tr!=end(); ++tr)
		for(size_t t=tr->start_time; t<=tr->last_time();++t)
			frameSizes[t]++;
	return frameSizes;
}



/** @brief make/reset the inverse index  */
void TrajIndex::makeInverse(const std::vector<size_t> &frameSizes)
{
    inverse.resize(0);
    inverse.reserve(frameSizes.size());
    for(size_t t=0; t<frameSizes.size(); ++t)
		inverse.push_back(new vector<size_t>(frameSizes[t]));
    for(size_t tr=0;tr<size();++tr)
        for(size_t t = (*this)[tr].start_time;t<=(*this)[tr].last_time();++t)
            inverse[t][(*this)[tr][t]]=tr;
}

/** @brief fill the TrajIndex with the content of a file stream  */
istream & Colloids::operator>>(std::istream& is, TrajIndex& tri)
{
    size_t start;
    string indexString;
    is>>start; //sliding the reading in order to avoid reading final empty line
    while(is.good())
    {
        is.get(); //escape the endl
        Traj tr(start);

        getline(is,indexString);
        istringstream indexStream(indexString,istringstream::in);
        indexStream>>tr;
        tri.push_back(tr);
        is>>start;
    }
    return is;
}

/** @brief operator<<
  *
  * @todo: document this function
  */
ostream & Colloids::operator<<(std::ostream& os, const TrajIndex& tri)
{
    copy(
        tri.begin(), tri.end(),
        ostream_iterator<Traj>(os, "\n")
        );
    return os;
}


