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

 * \file serieTracker.cpp
 * \brief Implement SerieTracker and the display routines of Tracker
 * \author Mathieu Leocmach
 * \date 11 January 2010
 *
 * Cimg is heavy to compile. That is why all functions needing it are set apart.
 *
 */
#include "serieTracker.hpp"
#include <CImg.h>

using namespace std;
using namespace cimg_library;
using namespace Colloids;
typedef boost::multi_array_types::index_range range;

/** \brief display the 3D image as it is.

    Large data may cause memory trouble
*/
void Tracker::display(const std::string &windowName) const
{
    if(!padded)
        //construct a cimg image sharing the tracker image data
        //imediately display it
        CImg<float>(
            data,
            centersMap.shape()[fortran_order?2:0],
            centersMap.shape()[1],
            centersMap.shape()[fortran_order?0:2],
            1,true
            ).display(windowName.c_str());
    else
        //padding is also displayed
        CImg<float>(
            data,
            fortran_order ? FFTmask.shape()[2]*2 : centersMap.shape()[0],
            centersMap.shape()[1],
            fortran_order ? centersMap.shape()[0] : FFTmask.shape()[2]*2,
            1,true
            ).display(windowName.c_str());
}

/** \brief Display only a slice of the image in the given dimension */
void Tracker::displaySlice(const size_t dim, const size_t position, const std::string &windowName) const
{
    boost::array<size_t,3> shape;
    copy(centersMap.shape(), centersMap.shape()+3, shape.begin());
    shape.back() = 2*(shape.back()/2 +1);
    //copy the slice's data contigiously
    boost::multi_array<float,3> thin_image
    (
        boost::const_multi_array_ref<float,3>(data, shape)
        [
            boost::indices
				[(dim==0) ? range(position, position+1) : range()]
				[(dim==1) ? range(position, position+1) : range()]
				[(dim==2) ? range(position, position+1) : range()]
		]
    );
    remove_copy(
        thin_image.shape(), thin_image.shape()+3,
        shape.begin(), 1u
        );
    CImg<float>(thin_image.origin(), shape[1], shape[0], 1, 1, true).display(windowName.c_str());
}

struct spectrum_norm : public unary_function<const complex<float>&, float>
{
    float operator()(const complex<float>&x) const {return log(1+abs(x));}
};

/** \brief display projections of the 3D spectrum.  */
void Tracker::displaySpectrum(const std::string &windowName) const
{
    using namespace boost;
    array<size_t,3> shape;
    copy(FFTmask.shape(), FFTmask.shape()+3, shape.begin());
    const_multi_array_ref<complex<float> ,3> c_spectrum((complex<float>*)data, shape);

    array<size_t,2> proj_shape;
    proj_shape[0] = FFTmask.shape()[0] + FFTmask.shape()[1];
    proj_shape[1] = FFTmask.shape()[1] + FFTmask.shape()[2];
    multi_array<complex<float> ,2> c_proj(proj_shape);
    multi_array<float,2> projections(proj_shape);

    c_proj[indices[range(0,FFTmask.shape()[0])][range(0,FFTmask.shape()[1])]]
        = c_spectrum[indices[range()][range()][0]];
    c_proj[indices[range().start(FFTmask.shape()[0])][range().start(FFTmask.shape()[1])]]
        = c_spectrum[indices[0][range()][range()]];
    c_proj[indices[range(0,FFTmask.shape()[0])][range().start(FFTmask.shape()[1])]]
        = c_spectrum[indices[range()][0][range()]];

    transform(
        c_proj.origin(), c_proj.origin()+c_proj.num_elements(),
        projections.origin(), spectrum_norm()
        );

    CImg<float>(
            projections.origin(),
            projections.shape()[1],
            projections.shape()[2],
            1,1,true
            ).display(windowName.c_str());
}

/** \brief display projections of the 3D mask.*/
void Tracker::displayMask(const std::string &windowName) const
{
    using namespace boost;
    array<size_t,2> proj_shape;
    proj_shape[0] = FFTmask.shape()[0] + FFTmask.shape()[1];
    proj_shape[1] = FFTmask.shape()[1] + FFTmask.shape()[2];
    multi_array<bool,2> projections(proj_shape);
    projections[indices[range(0,FFTmask.shape()[0])][range(0,FFTmask.shape()[1])]]
        = FFTmask[indices[range()][range()][0]];
    projections[indices[range().start(FFTmask.shape()[0])][range().start(FFTmask.shape()[1])]]
        = FFTmask[indices[0][range()][range()]];
    projections[indices[range(0,FFTmask.shape()[0])][range().start(FFTmask.shape()[1])]]
        = FFTmask[indices[range()][0][range()]];

    CImg<bool>(
            projections.origin(),
            projections.shape()[1],
            projections.shape()[2],
            1,1,true
            ).display(windowName.c_str());
}

/** \brief display the tracked centers.*/
void Tracker::displayCenters(const std::string &windowName) const
{
    CImg<bool>(
            centersMap.origin(),
            centersMap.shape()[fortran_order?2:0],
            centersMap.shape()[1],
            centersMap.shape()[fortran_order?0:2],
            1,true
            ).display(windowName.c_str());
}

/** @brief Constructor from an existant LifFile object  */
SerieTracker::SerieTracker(const std::string &namePattern, boost::array<size_t, 4> &xyzt, const double Zratio, const size_t ch, const unsigned fs)
{
    boost::array<size_t,3> dims;
    copy(xyzt.begin(), xyzt.begin()+3, dims.rbegin());
    tracker = new Tracker(dims, fs);
	tracker->fortran_order=true;
	setChannel(ch);
	setThreshold(0);
	this->ZXratio = Zratio;

	//parse the name pattern
	string head = namePattern, tail;
	size_t tdigits=0, zdigits=0,
        tpos = namePattern.rfind("_t"),
        zpos = namePattern.rfind("_z");

    if(tpos!=string::npos)
    {
        if(xyzt[3]>1)
            throw invalid_argument("Name pattern doesn't accept time dependence");
        this->hasTime=true;
        head.resize(tpos);
        tdigits = namePattern.find_first_not_of("0123456789", tpos+2);
        tail = namePattern.substr(tdigits);
        tdigits -= tpos;
    }
    else
        this->hasTime=false;
    if(zpos!=string::npos)
    {
        if(xyzt[2]>1)
            throw invalid_argument("Name pattern doesn't accept z dependence");
        this->hasDepth=true;
        head.resize(zpos);
        zdigits = namePattern.find_first_not_of("0123456789", zpos+2);
        //get the tail only if no time
        if(tail.empty())
            tail = namePattern.substr(zdigits);
        zdigits -= zpos;
    }
    else
        this->hasDepth=false;

    //reconstruct the format
    ostringstream os;
    os<<head;
    if(this->hasDepth)
        os<<"_z%|0"<<zdigits;
    if(this->hasTime)
        os<<"_t%|0"<<tdigits;
    os<<tail;
    this->serie.parse(os.str());

	this->centers=0;
	setTimeStep(0);

	return;
}

/** @brief set the current time step  */
void SerieTracker::setTimeStep(size_t t)
{
    this->time_step = t;
    if(centers)
    {
        Particles* old_centers = this->centers;
        this->centers = 0;
        delete old_centers;
    }
    CImg<unsigned char> buffer(tracker->centersMap.shape()[2],tracker->centersMap.shape()[1]);
    if(!hasTime)
    {
        if(!hasDepth)
        {
            buffer.load(serie.str().c_str());
            tracker->fillImage(buffer.data+buffer.offset(0,0,0,channel));
        }
        else
            for(size_t z=0; z<tracker->centersMap.shape()[0];++z)
            {
                buffer.load((serie%z).str().c_str());
                tracker->fillSlice(z, buffer.data+buffer.offset(0,0,0,channel));
            }
    }
    else
    {
        if(!hasDepth)
        {
            buffer.load((serie%t).str().c_str());
            tracker->fillImage(buffer.data+buffer.offset(0,0,0,channel));
        }
        else
            for(size_t z=0; z<tracker->centersMap.shape()[0];++z)
            {
                buffer.load((serie%z%t).str().c_str());
                tracker->fillSlice(z, buffer.data+buffer.offset(0,0,0,channel));
            }
    }
}



/** @brief Fill the tracker's image with the next time step  */
SerieTracker & SerieTracker::operator++()
{
    setTimeStep(this->time_step+1);
    return *this;
}
