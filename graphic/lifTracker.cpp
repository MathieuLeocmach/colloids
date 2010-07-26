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
 * \brief glue between Tracker and LifFile
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 22 december 2009
 *
 *
 */

#include <algorithm>
#include "lifTracker.hpp"

using namespace std;
using namespace Colloids;

/** @brief Constructor from an existant LifFile object  */
LifTracker::LifTracker(LifSerie &serie, const size_t ch, const unsigned fs)
{
    this->serie = &serie;
    setChannel(ch);
    tracker = new Tracker(getTrackerDims(), fs);
	tracker->fortran_order=true;
	this->centers=0;
	setTimeStep(0);
	setThreshold(0);
	return;
}

/** @brief chooseChannel  */
void LifTracker::chooseChannel()
{
	channel = getLif().chooseChannel();
}

/** @brief set the current time step  */
void LifTracker::setTimeStep(size_t t)
{
    this->iterator = serie->begin(t);
    this->time_step = t;
    if(centers)
    {
        Particles* old_centers = this->centers;
        this->centers = 0;
        delete old_centers;
    }
    tracker->fillImage_charToUchar(this->iterator);
}




/** @brief Fill the tracker's image with the next time step  */
LifTracker & LifTracker::operator++()
{
    this->iterator = tracker->fillImage_charToUchar(this->iterator);
    if(centers)
    {
        Particles* old_centers = this->centers;
        this->centers = 0;
        delete old_centers;
    }
    this->time_step++;
    return *this;
}




/** @brief get the dimensions according to the content of the lif serie
	Convert the dimension order from row major (Leica files) to column major (c order)
	If less than 3 dimensions, the first(s) dimension(s) is/are set to 1.
	That way (last dim != 1), real to complex FFT is efficient.
*/
boost::array<size_t,3> LifTracker::getTrackerDims() const
{
	boost::array<size_t,3> dims = {1,1,1};
	vector<size_t> fortran_order_dims = getLif().getSpatialDimensions();
	copy(
		fortran_order_dims.begin(),
		fortran_order_dims.end(),
		dims.rbegin()
		);
	return dims;
}


/** @brief get the dimensions according to the content of the lif serie
	Convert the dimension order from row major (Leica files) to column major (c order)
	If less than 3 dimensions, the first(s) dimension(s) is/are set to 1.
	That way (last dim != 1), real to complex FFT is efficient.
*/
boost::array<size_t,3> LifTracker2D::getTrackerDims() const
{
	boost::array<size_t,3> dims = {1,1,1};
	vector<size_t> fortran_order_dims = getLif().getSpatialDimensions();
	fortran_order_dims.back()=1;
	copy(
		fortran_order_dims.begin(),
		fortran_order_dims.end(),
		dims.rbegin()
		);
	return dims;
}


