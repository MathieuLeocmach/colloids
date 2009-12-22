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

/** \brief feed an initialized trajectory with the position indices contained by an input stream */
istream& operator>> (istream& is, Traj& tr )
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
ostream& operator<< (ostream& os, const Traj& tr )
{
    os << tr.start_time << endl;
    os << tr.steps[0];
    for(size_t i=1;i<tr.steps.size();++i)
        os << "\t" << tr.steps[i];
    return os;
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

/** \brief constructor from file */
TrajIndex::TrajIndex(const string &filename)
{
    //extract the path from the filename
    const size_t endfolder = filename.find_last_of("/\\");
    const string folder = filename.substr(0,endfolder+1);

    ifstream input(filename.c_str(), ios::in);
    if(input)
    {
        //header
        input >> radius >> dt;
        input.get(); //escape the endl

        //data of the file serie containing the positions
        string base_name,token;
        getline(input,base_name);
        getline(input,token);
        vector<string> tokens(1,token);
        tt = TokenTree(tokens,folder+base_name);
        //cout << tt <<endl;

        input >> t_offset >> t_size;

        size_t start;
        string indexString;
        input>>start; //sliding the reading in oder to avoid reading final empty line
        while(!input.eof())
        {
            input.get(); //escape the endl
            Traj tr(start);

            getline(input,indexString);
            istringstream indexStream(indexString,istringstream::in);
            indexStream>>tr;
            push_back(tr);
            input>>start;
        }
        input.close();
        makeInverse();
    }
    else
        throw invalid_argument("No such file as "+filename);

	return;
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

/** @brief make/reset the inverse index
  */
void TrajIndex::makeInverse()
{
    inverse.assign(getMaxTime()+1,vector<size_t>(size()));
    for(size_t tr=0;tr<size();++tr)
        for(size_t t = (*this)[tr].start_time;t<=(*this)[tr].last_time();++t)
            inverse[t][(*this)[tr][t]]=tr;
}


/** @brief Interrupt or attribute a follower to all the trajectories existing at lastTime. Create new trajectories at lastTime+1 if needed.
  * \param lastTime The time step to grow from
  * \param frameSize The total number of elements existing at lastTime+1
  * \param followersByDist The possible followers of each trajectory, sorted by their distence to the last position of the trajectory
  * Any distence in the topological meaning is ok : d(a,b)=0 <=> a=b and for all a!=b d(a,b)>0.
  */
void TrajIndex::addTimeStep(const size_t &lastTime, const size_t &frameSize, std::vector< std::multimap<double, size_t> > &followersByDist)
{
    //First we treat the simple cases : the trajectories not sharing their nearest possible follower
    const size_t notLinked = size();
    vector<size_t> linkedTo(frameSize,notLinked);
    vector<bool> conflict(size(),false), linked(frameSize,false);
    size_t nN;
    for(size_t tr=0;tr<size();++tr)
        if(at(tr).exist(lastTime) && !followersByDist[tr].empty())
        {
            nN = (*followersByDist[tr].begin()).second;
            if(linkedTo[nN]==notLinked)
                linkedTo[nN]=tr;
            else
            {
                conflict[tr]=true;
                conflict[linkedTo[nN]]=true;
            }
        }
    //simple and difficult cases told appart
    //now adding the non conflicting particles to the non conflicting trajectories
    for(size_t p = 0;p<frameSize;++p)
        if(linkedTo[p]!=notLinked && !conflict[linkedTo[p]])
        {
            at(linkedTo[p]).push_back(p);
            linked[p]=true;
        }

    //clean up followersByDist :
    // - Removing the follower ranking for the used trajectories
    // - Removing the new frame elements that where linked to a trajectory
    for(size_t tr=0;tr<size();++tr)
    {
        if(!conflict[tr])
            followersByDist[tr].clear();
        else
        {
            deque< multimap<double, size_t>::iterator > toBeRemoved;
            for(multimap<double, size_t>::iterator f=followersByDist[tr].begin();f!=followersByDist[tr].end();++f)
                if(linked[(*f).second])
                    toBeRemoved.push_back(f);
            for(deque< multimap<double, size_t>::iterator >::iterator r=toBeRemoved.begin();r!=toBeRemoved.end();++r)
                followersByDist[tr].erase(*r);
        }
    }
    //difficult case : two or more trajectories have the same nearest possible follower

    //index of difficult trajectories sorted by square distance to their nearest follower
    multimap<double,size_t> difficultTraj;
    for(size_t tr = 0;tr<size();++tr)
        if(!followersByDist[tr].empty())
            difficultTraj.insert(make_pair((*followersByDist[tr].begin()).first,tr));

    //cout << difficultTraj.size() << " difficult cases ... ";

    //link the difficult trajectories one by one, updating the difficult followers
    for(multimap<double,size_t>::iterator tr= difficultTraj.begin();tr!=difficultTraj.end();++tr)
        if(!followersByDist[(*tr).second].empty())
        {
            //index of the particle closest to the trajectory's last position
            size_t p = (*followersByDist[(*tr).second].begin()).second;
            //link the particle to the trajectory
            at((*tr).second).push_back(p);
            linked[p]=true;

            //remove p as possible followers of any (difficult) trajectory
            multimap<double,size_t>::iterator nexttr = tr;
            nexttr++;
            deque< multimap<double,size_t>::iterator > had_p_asFollower;
            while(nexttr!=difficultTraj.end())
            {
                for(multimap<double,size_t>::iterator placeOfp = followersByDist[(*nexttr).second].begin();placeOfp!=followersByDist[(*nexttr).second].end();++placeOfp)
                    if((*placeOfp).second==p)
                    {
                        followersByDist[(*nexttr).second].erase(placeOfp);
                        had_p_asFollower.push_back(nexttr);
                        break;
                    }
                nexttr++;
            }

            //re-insert in the ranking all the modified trajectories
            multimap<double,size_t> fol;
            for(deque<multimap<double,size_t>::iterator>::iterator f=had_p_asFollower.begin();f!=had_p_asFollower.end();++f)
            {
                const size_t actualTraj = (*(*f)).second;
                fol = followersByDist[actualTraj];
                difficultTraj.erase(*f);
                //re-insert if the follower list is not empty
                if(!fol.empty())
                    difficultTraj.insert(make_pair((*fol.begin()).first,actualTraj));
            }
        }

    //some elements of the new frame may not have found a trajectory to be linked to. Creates new a new trajectory for each of them.
    size_t nbnew=0;
    for(size_t p =0;p<frameSize;++p)
        if(!linked[p])
        {
            push_back(Traj(lastTime+1,p));
            nbnew++;
        }
    //cout << nbnew << " new trajectories" << endl;
}

