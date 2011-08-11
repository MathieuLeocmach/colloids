/*
 * liffinder.cpp
 *
 *  Created on: 11 août 2011
 *      Author: mathieu
 */

#include "locatorfromlif.h"

namespace Colloids {

LocatorFromLif::LocatorFromLif(LifSerie * serie):
		serie(serie), input(serie->begin()), t(0), z(0)
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
    	if(this->t >= this->total_t || this->z >= this->dims[2])
    		throw std::out_of_range("End of file reached");

    }

}
