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

 * \file radiiTracker.cpp
 * \brief Implement classes to extract particle radii from 3D image data
 * \author Mathieu Leocmach
 *
 */

#include "radiiTracker.hpp"
#include "particles.hpp"

using namespace std;
using namespace Colloids;
typedef boost::multi_array_types::index_range range;

/** @brief get a 3D Neighbourhood around a point and the coordinates of this point in the neighbourhood.  */
Coord RadiiTracker::getNeighbourhood(const Coord &point, array_type &ngb) const
{
    //initialize the neighbourhood
    fill_n(ngb.origin(), ngb.num_elements(), 0);
    boost::array<size_t, 3> low, high;
    Coord center(3);
    for(int d=0; d<3; ++d)
    {
        low[d] = (size_t)max(0.0, point[d] - ngb.shape()[d]/2);
        high[d] = (size_t)min((double)(data.shape()[d]), point[d] + (ngb.shape()[d]/2)+1);
        center[d] = point[d]-low[d];
    }
    ngb[boost::indices
				[range(0, high[0]-low[0])]
				[range(0, high[1]-low[1])]
				[range(0, high[2]-low[2])]
			] = data[boost::indices
				[range(low[0], high[0])]
				[range(low[1], high[1])]
				[range(low[2], high[2])]
			];
    return center;
}

/** @brief get the radius of a particle  */
double RadiiTracker::getRadius(const array_type &ngb, const Coord &center, const double &rmax)
{
    vector<double> val(N);
    vector<size_t> nbs(N);
    //bin linearly
    const double dmax = pow(rmax, 2);
    const unsigned char *px = ngb.origin();

    if(fortran_order)
        for(size_t i=0; i<ngb.shape()[0]; ++i)
            for(size_t j=0; j<ngb.shape()[1]; ++j)
                for(size_t k=0; k<ngb.shape()[2]; ++k)
                {
                    const double d = pow((i-center[0])*ZXratio, 2) + pow(j-center[1], 2) + pow(k-center[2], 2);
                    if(d<dmax)
                    {
                        size_t l = sqrt(d)/precision;
                        val[l] += *px;
                        nbs[l]++;
                    }
                    px++;
                }
    else
        for(size_t i=0; i<ngb.shape()[0]; ++i)
            for(size_t j=0; j<ngb.shape()[1]; ++j)
                for(size_t k=0; k<ngb.shape()[2]; ++k)
                {
                    const double d = pow(i-center[0], 2) + pow(j-center[1], 2) + pow((k-center[2])*ZXratio, 2);
                    if(d<dmax)
                    {
                        size_t l = sqrt(d)/precision;
                        val[l] += *px;
                        nbs[l]++;
                    }
                    px++;
                }
    //prevent division by 0
    nbs[0]=1;

    //integrates val and nbs (in place)
    partial_sum(val.begin(), val.end(), val.begin());
    partial_sum(nbs.begin(), nbs.end(), nbs.begin());

    //divide sum_intensity(r) by nb_pixels(r) => mean_intensity(r)
    transform(
        val.begin(), val.end(),
        nbs.begin(), val.begin(),
        divides<double>()
        );

    //smooth curve by convolving a gaussian via fft
    copy(val.begin(), val.end(), in);
    fftwf_execute(forward);

    for(ssize_t i=0; i<N/2+1; ++i)
        *reinterpret_cast<complex<float>*>(out+i) *= exp(-i/13.5);

    fftwf_execute(backward);
    copy(in,in+N, val.begin());


    //take derivative
    adjacent_difference(val.begin(), val.end(), val.begin());

    //find global minimum further than rmin
    const size_t margin = rmin/precision;
    return precision*distance(val.begin(), min_element(val.begin()+margin, val.end()-margin));
}

/** @brief get the radius of each particle  */
vector<double> RadiiTracker::getRadii()
{
    vector<double> radii(this->centers.size());

    boost::array<size_t, 3> ngb_sizes;
    for(size_t d=0; d<3; ++d)
        ngb_sizes[d] = (data.shape()[d]>2*rmax)?(2*rmax+1):1;

    //#pragma omp parallel for
    for(size_t p=0; p<this->centers.size(); ++p)
    {
        array_type ngb(ngb_sizes);
        radii[p] = getRadius(ngb, getNeighbourhood(centers[p], ngb), rmax);
    }

    return radii;
}

/** \brief helper functor to transform zyx into xyz*/
template<class Coordinates>
struct revert : public unary_function<Coordinates& ,void>
{
	void operator()(Coordinates& v){swap(v[0], v[2]);}
};

/** @brief Translate the input coordinates to voxel sizes  */
void RadiiTracker::setCenters(Particles c)
{
    centers = c;
    if(fortran_order)
	{
		//If the image was filled in row major ordering of the dimensions, the coordinates are given as z,y,x by the tracker.
		//Here we convert the human-readable x,y,z coordinates to the machine-readable z,y,x
		for_each(centers.begin(), centers.end(), revert<valarray<double> >());
	}

}

/** @brief resize  */
void RadiiTracker::resize()
{
    #pragma omp critical (fftw_create_plan)
    {
        N = rmax/precision;
        in = (float*) fftwf_malloc(sizeof(float) * N);
        out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * (N/2+1));
        forward = fftwf_plan_dft_r2c_1d(N, in, out, flags);
        backward = fftwf_plan_dft_c2r_1d(N, out, in, flags);
    }
}



/** @brief Free the memory and leave the iterator as the default constructed one. Tracker get deleted  */
void RadiiTrackerIterator::close()
{
    if(tracker) delete tracker;
    if(radii) delete radii;
}

/** @brief set radii tracker's centers by scaling z coordinate back to pixel number */
void RadiiTrackerIterator::setCenters(Particles c)
{
    Coord v(1.0,3);
    v[2] = 1/getZXratio();
    c *= v;
    if(!quiet()) cout <<"z divided by "<<getZXratio()<<endl;
    tracker->setCenters(c);
}



/** @brief find radii and return the result.
  *
  * Contrary to the dereferencement operator for standard iterator, this is a heavy operation.
  * The result is cached until the next incrementation or setTimeStep
  */
vector<double>& RadiiTrackerIterator::operator*()
{
    if(!radii)
    {
        if(!quiet()) cout<<"radii at t="<< this->time_step <<endl;
        radii = new vector<double>(tracker->getRadii());
        if(!quiet()) cout<<"done !"<<endl;
    }
    return *radii;
}

/** @brief Constructor from an existant LifFile object  */
LifRadiiTracker::LifRadiiTracker(LifSerie &serie, const std::string &centerSeriePrefix, const size_t ch, const unsigned fs)
{
    this->serie = &serie;
    setChannel(ch);
    tracker = new RadiiTracker(getTrackerDims(), fs);
	tracker->fortran_order=true;
	this->radii=0;
	setCentersSerie(FileSerie(
        FileSerie::get0th(centerSeriePrefix, getLif().getNbTimeSteps())+".dat",
        "_t", getLif().getNbTimeSteps(), 0
        ));
    cout<<"reading centers from "<<FileSerie::get0th(centerSeriePrefix, getLif().getNbTimeSteps())<<".dat"<<endl;
	setTimeStep(0);
	return;
}

/** @brief chooseChannel  */
void LifRadiiTracker::chooseChannel()
{
	channel = getLif().chooseChannel();
}

/** @brief set the current time step  */
void LifRadiiTracker::setTimeStep(size_t t)
{
    this->iterator = serie->begin(t);
    this->time_step = t;
    if(radii)
    {
        vector<double>* old_radii = this->radii;
        this->radii = 0;
        delete old_radii;
    }
    tracker->fillImage_charToUchar(this->iterator);
    assert(centersSerie);
    if(this->time_step < getLif().getNbTimeSteps())
        setCenters(Particles(*centersSerie % this->time_step));
}




/** @brief Fill the tracker's image with the next time step  */
LifRadiiTracker & LifRadiiTracker::operator++()
{
    this->iterator = tracker->fillImage_charToUchar(this->iterator);
    if(radii)
    {
        vector<double>* old_radii = this->radii;
        this->radii = 0;
        delete old_radii;
    }
    this->time_step++;
    assert(centersSerie);
    if(this->time_step < getLif().getNbTimeSteps())
        setCenters(Particles(*centersSerie % this->time_step));
    return *this;
}




/** @brief get the dimensions according to the content of the lif serie
	Convert the dimension order from row major (Leica files) to column major (c order)
	If less than 3 dimensions, the first(s) dimension(s) is/are set to 1.
	That way (last dim != 1), real to complex FFT is efficient.
*/
boost::array<size_t,3> LifRadiiTracker::getTrackerDims() const
{
	boost::array<size_t,3> dims = {{1,1,1}};
	vector<size_t> fortran_order_dims = getLif().getSpatialDimensions();
	copy(
		fortran_order_dims.begin(),
		fortran_order_dims.end(),
		dims.rbegin()
		);
	return dims;
}
