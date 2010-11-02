/**
    Copyright 2010 Mathieu Leocmach

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

 * \file serieTracker.hpp
 * \brief Define a class to track series of 2D image files
 * \author Mathieu Leocmach
 * \date 11 January 2010
 *
 *
 */
#ifndef serie_tracker_H
#define serie_tracker_file_H

#include "tracker.hpp"
#include <boost/format.hpp>
namespace Colloids{
/** \brief glue between Tracker and image files
    Can be used as a stadard InputIterator
*/
class SerieTracker : public TrackerIterator
{
 	public:
		explicit SerieTracker(
            const std::string &namePattern, boost::array<size_t, 4> &xyzt,
            const double Zratio=1.0,
            const size_t ch=0,
            const unsigned fs=FFTW_ESTIMATE);
		SerieTracker end(){return SerieTracker(length);};
		bool reachedEnd(){return getTimeStep()>=length;};

		void setTimeStep(size_t t);
		double getZXratio(){return ZXratio;};

		SerieTracker& operator++();
		bool operator==(const SerieTracker& rhs) {return (getPattern()==rhs.getPattern()) && (time_step==rhs.time_step);}
		//bool operator!=(const SerieTracker& rhs) {return (getPattern()!=rhs.getPattern()) || (time_step!=rhs.time_step);}

		std::string getPattern() const;

    private:
        boost::format serie;
        size_t length;
        double ZXratio;
        bool hasTime, hasDepth;

        /** \brief default constructor that should be used only to get an "end" iterator*/
        SerieTracker(const size_t end_step) {time_step = end_step;};

};
}
#endif
