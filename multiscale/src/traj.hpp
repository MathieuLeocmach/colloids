/**
    Copyright 2008-2011 Mathieu Leocmach

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
 * \date 3 March 2009
 *
 */

#ifndef traj_H
#define traj_H

#include <boost/ptr_container/ptr_vector.hpp>
#include <map>
#include <vector>
#include <deque>
#include <list>
#include <iterator>
#include <iostream>
#include <stdexcept>

namespace Colloids
{
    /**
     * \brief Trajectory object
     * A trajectory is an index linking the time step to a position index (typically the index of a position inside a "Particles" object).
     * The trajectory is defined by a starting time step and a list of position indexes
    */
    class Traj
    {
    	int start_time;
		std::deque<size_t> steps;
        public :
            explicit Traj(const int &start_time, const size_t &start_position):
				start_time(start_time), steps(1, start_position){};

            const int& get_start() const {return this->start_time;}
            int get_finish() const {return start_time + steps.size();}
            bool exist(const int &t) const {return start_time<=t && t<get_finish();}
            bool span(const int &t0,const int &t1) const {return start_time<=t0 && t1<=get_finish();}
            size_t size() const {return steps.size();};
            const size_t &operator[](const int &t) const {assert(exist(t));return steps[t-start_time];};
            void update_last(const size_t& p){steps.back()=p;};
            void push_back(const size_t &pos){steps.push_back(pos);};
    };

    /**
	 * \brief Trajectory index links "frames" of instantaneous objects (ex: positions) into persisting objects (ex:trajectories).
	 * Frames can be 2D slices of a 3D data set and then a given trajectory links the positions of slices of the same 3D object.
	*/
    class TrajIndex : boost::noncopyable
    {
        private:
    	boost::ptr_vector<Traj> tr2pos;
        boost::ptr_vector< std::vector<size_t> >pos2tr;

        public:
        TrajIndex(const size_t& nb_initial_positions);

		//accessors
        size_t size(void) const {return tr2pos.size();};
        size_t nbFrames(void) const {return pos2tr.size();};
		const Traj& operator[](const size_t &tr) const {return tr2pos[tr];}
		const size_t & getTraj(const size_t &t,const size_t &p) const {return pos2tr[t][p];};
		const std::vector<size_t>& getInverse(const size_t &t) const {return pos2tr[t];};


		//processing
		void add_Frame(const size_t &frame_size, const std::vector<double> &distances, const std::vector<size_t> &p_from, const std::vector<size_t> &p_to);
		template<class InputIterator>
		void add_Frame(InputIterator first, InputIterator last);

    };

    template<class InputIterator>
    void TrajIndex::add_Frame(InputIterator first, InputIterator last)
    {
    	//check for uniqueness
    	std::list<size_t> input(first, last);
    	const size_t original_size = input.size();
    	input.sort();
    	input.unique();
    	if(input.size() < original_size)
    		throw std::invalid_argument("TrajIndex::add_Frame: a trajectory index appear more than once.");
    	//limit between old and new trajectories
    	const size_t nbtrajs = this->size();
    	std::list<size_t>::const_iterator
			old_end = std::upper_bound(input.begin(), input.end(), nbtrajs-1);
    	//check that old trajectories can be continued
    	for(std::list<size_t>::const_iterator it = input.begin(); it!=old_end; ++it)
			if(this->tr2pos[*it].get_finish() < this->pos2tr.size())
				throw std::invalid_argument("TrajIndex::add_Frame: cannot continue a trajectory that do not exist in the previous frame");
    	//check that the new trajectories are consecutive
    	if((old_end != input.end()) && (input.back() - *old_end +1 != std::distance(old_end, (std::list<size_t>::const_iterator)input.end())))
    		throw std::invalid_argument("TrajIndex::add_Frame: indices of new trajectories are not consecutive");

    	//create all necessary new trajectories
		for(std::list<size_t>::const_iterator it = old_end; it!=input.end(); ++it)
			this->tr2pos.push_back(new Traj(this->pos2tr.size(), 0));
    	//fill pos2tr
    	this->pos2tr.push_back(new std::vector<size_t>(input.size()));
    	std::copy(first, last, this->pos2tr.back().begin());

    	//update tr2pos
    	for(size_t p=0; p<this->pos2tr.back().size(); ++p)
    	{
    		const size_t tr = this->pos2tr.back()[p];
			if(tr < nbtrajs)
				this->tr2pos[tr].push_back(p);
			else
				this->tr2pos[tr].update_last(p);
    	}

    }

};


#endif
