/*
 * liffinder.cpp
 *
 *  Created on: 11 août 2011
 *      Author: mathieu
 */

#include "locatorfromlif.h"

namespace Colloids {

LocatorFromLif::LocatorFromLif(LifSerie * serie):
		serie(serie), t(0), input(serie->begin())
{
	this->dims = serie->getSpatialDimensions();
	this->total_t = serie->getNbTimeSteps();
	this->finder.reset(new MultiscaleFinder2D(dims[0], dims[1]));
	this->slice.create(dims[0], dims[1]);
}

LocatorFromLif::~LocatorFromLif() {
	// TODO Auto-generated destructor stub
}

    void LocatorFromLif::fill_next_slice()
    {
    	std::vector<Center2D> centers;
    	if(this->t >= this->total_t || this->get_z() >= this->dims[2])
    		throw std::out_of_range("End of file reached");
    	this->serie->fill2DBuffer(static_cast<void*>(this->slice.data), this->t, this->get_z());
    	this->finder->get_centers(this->slice, centers);
    	this->rec.push_back(centers);
    }

}
