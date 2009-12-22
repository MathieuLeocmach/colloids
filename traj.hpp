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


 * \file traj.hpp
 * \brief Defines a trajectory class
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 3 March 2009
 *
 */

#ifndef traj_H
#define traj_H

#include <deque>
#include <map>
#include <istream>
#include <stdexcept>
#include "files_series.hpp"

class TrajError : std::exception
{
    public:
        size_t asked,start,stop;

        TrajError(const size_t &ask,const size_t &sta,const size_t &sto){asked=ask;start=sta;stop=sto;return;}
        virtual const char* what() const throw();
};
class IdTrajError : TrajError
{
    public :
        size_t id;
        IdTrajError(const TrajError &tre,const size_t &i) : TrajError(tre){id=i;return;};
        const char* what() const throw();
};

/**
    \brief Trajectory object
    A trajectory is an indexation linking the time step to a position index (typically the index of a position inside a "Particles" object).
    The trajectory is defined by a starting time step and a list of position indexes
*/
class Traj
{
    public :
        size_t start_time;
        std::deque<size_t> steps;

        Traj(const size_t &start);
        Traj(const size_t &start,const size_t &first_step);
        //Traj(Traj &tr,size_t t0,size_t tmax);

        size_t last_time() const;
        bool exist(const size_t &t) const;
        bool span(const size_t &t0,const size_t &t1) const;
        size_t &operator[](const size_t &t) throw (TrajError);
        const size_t &operator[](const size_t &t) const throw (TrajError);
        void push_back(const size_t &pos);
};

std::istream& operator>> (std::istream& is, Traj& tr );
std::ostream& operator<< (std::ostream& os, const Traj& tr );

class TrajIndex : public std::deque<Traj>
{
    public:
        TokenTree tt;
        size_t t_offset, t_size;
        double radius,dt;
        std::vector< std::vector<size_t> >inverse;

        TrajIndex():std::deque<Traj>(){return;};
        TrajIndex(const std::string &filename);

        void addTimeStep(const size_t &lastTime, const size_t &frameSize, std::vector< std::multimap<double, size_t> > &followersByDist);

        size_t & getTraj(const size_t &t,const size_t &p);
        void makeInverse();
        size_t getMaxTime() const;

        template<typename T>
        std::map<size_t,T>& trajToPos(const size_t &t,const std::map<size_t,T> &trajMap,std::map<size_t,T> &posMap) const;
        template<typename T>
        std::map<size_t,T>& posToTraj(const size_t &t,const std::map<size_t,T> &posMap,std::map<size_t,T> &trajMap) const;


		struct Converter : public std::unary_function<const size_t&,size_t>
		{
			const size_t t;
			const TrajIndex *const index ;
			Converter(const size_t &time, const TrajIndex &i) : t(time), index(&i){};

			/** gives the position index of a trajectory in the frame t*/
			size_t operator()(const size_t &in) const
			{
				return (*index)[in][t];
			}
		};
		struct Inverser : public std::unary_function<const size_t&,size_t>
		{
			const size_t t;
			const TrajIndex *const index ;
			Inverser(const size_t &time, const TrajIndex &i) : t(time), index(&i){};

			/** gives the trajectory index of a position in the frame t*/
			size_t operator()(const size_t &in) const
			{
				return index->inverse[t][in];
			}
		};
};
/** \brief easy constructor */
inline Traj::Traj(const size_t &start)
{
    start_time=start;
    steps = std::deque<size_t>(0);
    return;
}

/** \brief constructor from first data */
inline Traj::Traj(const size_t &start,const size_t &first_step)
{
    start_time=start;
    steps = std::deque<size_t>(1,first_step);
    return;
}

/** \brief access to the particle index at a given time step */
inline size_t Traj::last_time() const
{
    return start_time+steps.size()-1;
}

/** \brief access to the particle index at a given time step */
inline size_t& Traj::operator[](const size_t &t) throw (TrajError)
{
    if (!exist(t)) throw TrajError(t,start_time,last_time());
    return steps[t-start_time];
}

/** \brief access to the particle index at a given time step */
inline const size_t& Traj::operator[](const size_t &t) const throw (TrajError)
{
    if (!exist(t)) throw TrajError(t,start_time,last_time());
    return steps[t-start_time];
}

/** \brief is the trajectory existing at time t ? */
inline bool Traj::exist(const size_t &t) const
{
    if(t<start_time || t>last_time()) return false;
    else return true;
}

/** \brief does the trajectory span at least from t0 to t1 ? */
inline bool Traj::span(const size_t &t0,const size_t &t1) const
{
    if(t0<start_time || t1>last_time() ) return false;
    else return true;
}

/** \brief add a new position index at the end of the trajectory */
inline void Traj::push_back(const size_t &pos)
{
    steps.push_back(pos);
}

/** @brief convert a mapping on positions in frame t to a mapping to trajectories  */
template<typename T>
std::map<size_t,T> & TrajIndex::posToTraj(const size_t &t,const std::map<size_t,T> &posMap,std::map<size_t,T> &trajMap) const
{
	for(typename std::map<size_t,T>::const_iterator it = posMap.begin();it!=posMap.end();++it)
		trajMap.insert(trajMap.end(),make_pair(inverse[t][(*it).first],(*it).second));
	return trajMap;
}

/** @brief convert a mapping on trajectories to a mapping on positions in frame t  */
template<typename T>
std::map<size_t,T> & TrajIndex::trajToPos(const size_t &t,const std::map<size_t,T> &trajMap,std::map<size_t,T> &posMap) const
{
	for(typename std::map<size_t,T>::const_iterator it = trajMap.begin();it!=trajMap.end();++it)
	{
		//std::cout<<(*it).first;
		if((*this)[(*it).first].exist(t))
		{
			//std::cout<<"\texist at t="<<t<<"\t"<<*(double*)&(*it).second;
			posMap.insert(posMap.end(),make_pair(at((*it).first)[t],(*it).second));
		}
		//std::cout<<std::endl;
	}
	return posMap;
}



#endif
