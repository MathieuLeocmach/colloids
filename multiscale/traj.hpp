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

//#include "index.hpp"
//#include "files_series.hpp"

//#include <boost/bimap.hpp>
//#include <boost/format.hpp>

#include <map>
#include <vector>
#include <deque>
#include <boost/ptr_container/ptr_vector.hpp>
//#include <istream>
//#include <stdexcept>

namespace Colloids
{
    /*class TrajError : public std::exception
    {
        public:
            size_t asked,start,stop;

            explicit TrajError(const size_t &ask,const size_t &sta,const size_t &sto){asked=ask;start=sta;stop=sto;return;}
            virtual const char* what() const throw();
    };
    class IdTrajError : public TrajError
    {
        public :
            size_t id;
            explicit IdTrajError(const TrajError &tre,const size_t &i) : TrajError(tre){id=i;return;};
            const char* what() const throw();
    };*/

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
            void push_back(const size_t &pos){steps.push_back(pos);};
    };

    class TrajIndex
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



		size_t getMaxTime() const;
		size_t longest_span() const;
		std::vector<size_t> getFrameSizes(const size_t &length=0) const;
    };

};


#endif
