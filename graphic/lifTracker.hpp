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

 * \file lifTracker.hpp
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 22 december 2009
 *
 *
 */

#ifndef lif_tracker_H
#define lif_tracker_file_H

#include "tracker.hpp"
#include "lifFile.hpp"
namespace Colloids{
/** \brief glue between Tracker and LifFile
Can be used as a stadard InputIterator
*/
class LifTracker : public TrackerIterator
{
 	LifSerie *serie;
 	std::istreambuf_iterator<char> iterator;

 	public:
		explicit LifTracker(LifSerie &serie, const size_t ch=0, const unsigned fs=FFTW_ESTIMATE);
		/** \brief default constructor that should be used only to get an "end" iterator*/
		LifTracker end(){return LifTracker(getLif().getNbTimeSteps());};

		LifSerie& getLif() const {return *serie;};

		void chooseChannel();
		void setTimeStep(size_t t);
		double getZXratio(){return serie->getZXratio();};

		LifTracker& operator++();
		bool operator==(const LifTracker& rhs) {return time_step==rhs.time_step && serie==rhs.serie;}
		//bool operator!=(const LifTracker& rhs) {return time_step!=rhs.time_step;}


	private:
        LifTracker(const size_t end_step) : iterator(){time_step = end_step;};
		boost::array<size_t,3> getTrackerDims() const;
		std::streampos tellg(){return serie->tellg();}

};
}
#endif
