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
	this->total_z = this->dims.size()>2 ? this->dims[2] : 1;
	this->total_t = serie->getNbTimeSteps();
	this->finder.reset(new MultiscaleFinder2D(dims[0], dims[1]));
	this->slice.create(dims[0], dims[1]);
	this->slice.setTo(0);
}

LocatorFromLif::~LocatorFromLif() {
	// TODO Auto-generated destructor stub
}

	void LocatorFromLif::clear()
	{
		this->slice.setTo(0);
		this->rec.clear();
	}

	/** \brief Read the next 2D slice from the lif file, track it and add it to the reconstructor */
    void LocatorFromLif::fill_next_slice()
    {
    	std::vector<Center2D> centers;
    	if(this->t >= this->total_t || this->get_z() >= this->total_z)
    		throw std::out_of_range("End of file reached");
    	this->serie->fill2DBuffer(static_cast<void*>(this->slice.data), 0, this->get_z());//this->t, this->get_z());
    	this->finder->get_centers(this->slice, centers);
    	this->rec.push_back(centers);
    }

    void LocatorFromLif::fill_time_step()
    {
    	//ensure that we start at the beginning of the stack
    	this->clear();
    	//track each slice
    	for(size_t z=0; z<this->total_z; ++z)
    		this->fill_next_slice();
    }
    void LocatorFromLif::get_centers(Centers & centers)
	{

		//3D reconstruction
		//this->rec.split_clusters();
		this->rec.get_blobs(centers);
		//convert the z coordinate according to the z-spacing of acquisition
		for(Centers::iterator c=centers.begin(); c!=centers.end(); ++c)
			(*c)[2] *= this->serie->getZXratio();
		this->t++;
	}

}
