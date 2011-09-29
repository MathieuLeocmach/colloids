/*
 * multiscalefinder.cpp
 *
 *  Created on: 11 juil. 2011
 *      Author: mathieu
 */

#include "multiscalefinder.hpp"
#include <stdexcept>
#include <iostream>

using namespace std;

namespace Colloids {

	MultiscaleFinder2D::MultiscaleFinder2D(const int nrows, const int ncols, const int nbLayers, const double &preblur_radius)
	{
		//must not fail if the image is too small, construct the 0th octave anyway
		const int s = (nrows<12 || ncols<12)?1:(log(min(nrows, ncols)/6)/log(2));
		this->octaves.reserve((size_t)s);
		this->octaves.push_back(new OctaveFinder(2*nrows, 2*ncols, nbLayers, preblur_radius));
		int ocr = nrows, occ = ncols;
		while(ocr >= 12 && occ >= 12)
		{
			this->octaves.push_back(new OctaveFinder(ocr, occ, nbLayers, preblur_radius));
			ocr /= 2;
			occ /= 2;
		}
		this->small = Image(nrows, ncols, 0.0);
		this->upscaled = Image(2*nrows, 2*ncols, 0.0);
	}
	MultiscaleFinder1D::MultiscaleFinder1D(const int ncols, const int nbLayers, const double &preblur_radius)
	{
		//must not fail if the signal is too short, construct the 0th octave anyway
		const int s = (ncols<12)?1:(log(ncols/6)/log(2));
		this->octaves.reserve((size_t)s);
		this->octaves.push_back(new OctaveFinder1D(2*ncols, nbLayers, preblur_radius));
		int occ = ncols;
		while(occ >= 12)
		{
			this->octaves.push_back(new OctaveFinder1D(occ, nbLayers, preblur_radius));
			occ /= 2;
		}
		this->small = Image(1, ncols, 0.0);
		this->upscaled = Image(1, 2*ncols, 0.0);
	}

	MultiscaleFinder::~MultiscaleFinder() {
		while(!this->octaves.empty())
		{
			delete octaves.back();
			octaves.pop_back();
		}
	}

    void MultiscaleFinder::set_radius_preblur(const double & k)
    {
    	for(size_t o=0; o<this->octaves.size(); ++o)
    		this->octaves[o]->set_radius_preblur(k);
    }

    /**
     * \brief fill all the octaves
     */
    void MultiscaleFinder::fill(const cv::Mat & input)
    {
    	if(input.rows != (int)this->get_width())
    	{
    		std::cerr<< "input.rows="<<input.rows<<"\twidth="<<this->get_width()<<std::endl;
    		throw std::invalid_argument("MultiscaleFinder::fill : the input's rows must match the width of the finder");
    	}
    	if(input.cols != (int)this->get_height())
    	{
    		std::cerr<< "input.cols="<<input.cols<<"\theight="<<this->get_height()<<std::endl;
    	    throw std::invalid_argument("MultiscaleFinder::fill : the input's cols must match the height of the finder");
    	}

    	input.convertTo(small, this->small.type());
    	//upscale the input to fill the first octave
    	//cv::resize does not work with double input, so we do it by hand
    	this->upscale();
    	this->octaves[0]->preblur_and_fill(this->upscaled);
    	if(this->octaves.size()>1)
			//Octave 1 corresponds to the size of the input image.
			//To avoid errors in the upsampling+downsampling process, we use the input directly
			this->octaves[1]->preblur_and_fill(input);
    	//For higher octaves we use the second to last layer of the previous octave, downsampled
    	for(size_t o=2; o<this->octaves.size(); ++o)
    		this->octaves[o]->fill(this->downscale(o));
    }
    /**
     * \brief Locate centers with pixel resolutions and scale resolution
     */
	void MultiscaleFinder::initialize_binary()
	{
		//initialize binary for each octave
    	for(size_t o=0; o<this->octaves.size(); ++o)
			this->octaves[o]->initialize_binary();
    	//Remove pixel centers that exist in consecutive octaves
    	for(size_t o=0; o<this->octaves.size()-1; ++o)
    		this->octaves[o]->seam_binary(*this->octaves[o+1]);
	}

    const MultiscaleFinder::Image MultiscaleFinder2D::downscale(const size_t &o)
	{
    	//second to last Gaussian layer of octave o-1 has a blurring radius two time larger than the original
    	Image roi2 = small(
				cv::Range(0, this->octaves[o]->get_width()),
				cv::Range(0, this->octaves[o]->get_height())
		);
		const Image & a = this->octaves[o-1]->get_layersG(this->octaves[o-1]->get_n_layers());
		for(int j=0; j<roi2.cols && 2*j+1<a.cols; ++j)
			for(int i=0; i<roi2.rows && 2*i+1<a.rows; ++i)
				roi2(i,j) = (a(2*i, 2*j) + a(2*i+1, 2*j) + a(2*i, 2*j+1) + a(2*i+1, 2*j+1))/4.0;
		return roi2;
	}
    const MultiscaleFinder::Image MultiscaleFinder1D::downscale(const size_t &o)
	{
		//second to last Gaussian layer of octave o-1 has a blurring radius two time larger than the original
    	Image roi2 = small(
				cv::Range(0, this->octaves[o]->get_width()),
				cv::Range(0, this->octaves[o]->get_height())
		);
		const Image & a = this->octaves[o-1]->get_layersG(this->octaves[o-1]->get_n_layers());
		for(int i=0; i<roi2.cols; ++i)
			roi2(0, i) = (a(0, 2*i) + a(0, 2*i+1))/2.0;
		return roi2;
	}
    void MultiscaleFinder::upscale()
    {
    	for(int j=0; 2*j<this->upscaled.rows; ++j)
		{
			OctaveFinder::PixelType * u = &this->upscaled(2*j, 0);
			const OctaveFinder::PixelType * s = &this->small(j,0);
			for(int i=0; 2*i<this->upscaled.cols; ++i)
			{
				*u++ = *s++;
				u++;
			}
			u = &this->upscaled(2*j, 1);
			s = &this->small(j,0);
			for(int i=0; 2*i+1<this->upscaled.cols && i+1<small.cols; ++i)
			{
				*u = 0.5 * *s++;
				*u++ += 0.5 * *s;
				u++;
			}
		}
		for(int j=0; 2*j+1<this->upscaled.rows && j+1<small.rows; ++j)
		{
			OctaveFinder::PixelType * u = &this->upscaled(2*j+1, 0);
			const OctaveFinder::PixelType * s1 = &this->small(j,0);
			const OctaveFinder::PixelType * s2 = &this->small(j+1,0);
			for(int i=0; 2*i<this->upscaled.cols; ++i)
			{
				*u++ = 0.5 * (*s1++ + *s2++);
				u++;
			}
			u = &this->upscaled(2*j+1, 1);
			s1 = &this->small(j,0);
			s2 = &this->small(j+1,0);
			for(int i=0; 2*i+1<this->upscaled.cols && i+1<small.cols; ++i)
			{
				*u = 0.25 * (*s1++ + *s2++);
				*u++ += 0.25 * (*s1 + *s2);
				u++;
			}
		}
    }
}
