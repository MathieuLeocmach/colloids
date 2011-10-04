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
	MultiscaleFinder3D::MultiscaleFinder3D(const int nplanes, const int nrows, const int ncols, const int nbLayers, const double &preblur_radius)
	{
		//must not fail if the image is too small, construct the 0th octave anyway
		const int s = (nplanes<10 || nrows<10 || ncols<10)?1:(log(min(nplanes, min(nrows, ncols))/5)/log(2));
		this->octaves.reserve((size_t)s);
		this->octaves.push_back(new OctaveFinder3D(2*nplanes, 2*nrows, 2*ncols, nbLayers, preblur_radius));
		int opl = nplanes, ocr = nrows, occ = ncols;
		while(opl >=10 && ocr >= 10 && occ >= 10)
		{
			this->octaves.push_back(new OctaveFinder3D(opl, ocr, occ, nbLayers, preblur_radius));
			opl /= 2;
			ocr /= 2;
			occ /= 2;
		}
		//random file name in the working directory
		this->path.reserve(30);
		this->path.push_back('_');
		this->path.push_back('_');
		for(int i=0; i<28;++i)
			this->path.push_back('a'+rand()%('Z'-'a'));
		//create a memory mapped file to contain the images data
		boost::iostreams::mapped_file_params params(this->path);
		const size_t nbpixels =  nplanes * nrows * ncols;
		params.new_file_size = nbpixels * sizeof(OctaveFinder::PixelType) * 9;
		params.flags = boost::iostreams::mapped_file::readwrite;
		this->file.open(params);
		//create the images inside the memory mapped file
		int dims[3] = {nplanes, nrows, ncols};
		this->small = cv::Mat(3, dims, this->octaves[0]->get_layers(0).type(), (void*)this->file.data());
		this->small.setTo(0);
		int ldims[3] = {2*nplanes, 2*nrows, 2*ncols};
		this->upscaled = cv::Mat(3, ldims, small.type(), (void*)(this->file.data() + nbpixels * sizeof(OctaveFinder::PixelType)));
		this->upscaled.setTo(0);

	}

	MultiscaleFinder::~MultiscaleFinder() {
		while(!this->octaves.empty())
		{
			delete octaves.back();
			octaves.pop_back();
		}
	}
	MultiscaleFinder3D::~MultiscaleFinder3D() {
		remove(this->path.c_str());
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
    	if(input.size[input.dims-2] != (int)this->get_width())
    	{
    		std::cerr<< "input.rows="<<input.rows<<"\twidth="<<this->get_width()<<std::endl;
    		throw std::invalid_argument("MultiscaleFinder::fill : the input's rows must match the width of the finder");
    	}
    	if(input.size[input.dims-1] != (int)this->get_height())
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
    const MultiscaleFinder::Image MultiscaleFinder3D::downscale(const size_t &o)
	{
		//second to last Gaussian layer of octave o-1 has a blurring radius two time larger than the original
    	int dims[3] = {
    			dynamic_cast<OctaveFinder3D*>(this->octaves[o])->get_depth(),
    			this->octaves[o]->get_width(),
    			this->octaves[o]->get_height()};
		Image roi2 = cv::Mat(3, dims, small.type(), small.data);
		const Image & a = this->octaves[o-1]->get_layersG(this->octaves[o-1]->get_n_layers());
		for(int k=0; k<roi2.size[0] && 2*k+1<a.size[0]; ++k)
			for(int j=0; j<roi2.size[1] && 2*j+1<a.size[1]; ++j)
				for(int i=0; i<roi2.size[0] && 2*i+1<a.size[0]; ++i)
					roi2(k,j,i) = 0.125 * (
							a(2*k, 2*j, 2*i) + a(2*k, 2*j, 2*i+1) +
							a(2*k, 2*j+1, 2*i) + a(2*k, 2*j+1, 2*i+1) +
							a(2*k+1, 2*j, 2*i) + a(2*k+1, 2*j, 2*i+1) +
							a(2*k+1, 2*j+1, 2*i) + a(2*k+1, 2*j+1, 2*i+1)
							);
		return roi2;
	}
    void MultiscaleFinder::upscale()
    {
    	//fill the even lines
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
    	//fill the odd lines
    	for(int j=0; 2*j+2<this->upscaled.size[0]; ++j)
		{
			OctaveFinder::PixelType * u = &this->upscaled(2*j, 0),
					* v = &this->upscaled(2*j+1, 0),
					* w = &this->upscaled(2*j+2, 0);
			for(int i=0; i<this->upscaled.size[1]; ++i)
				*v++ = 0.5 * (*u++ + *w++);
		}
    }
    void MultiscaleFinder3D::upscale()
	{
    	//fill the even lines of the even planes
    	for(int k=0; k<this->small.size[0]; ++k)
    		for(int j=0; j<this->small.size[1]; ++j)
    		{
    			OctaveFinder::PixelType * u = &this->upscaled(2*k, 2*j, 0);
				const OctaveFinder::PixelType * s = &this->small(k, j, 0);
				for(int i=0; i<this->small.size[2]; ++i)
				{
					*u++ = *s++;
					u++;
				}
				u = &this->upscaled(2*k, 2*j, 1);
				s = &this->small(k, j, 0);
				for(int i=0; 2*i+2<this->upscaled.size[2]; ++i)
				{
					*u = 0.5 * *s++;
					*u++ += 0.5 * *s;
					u++;
				}
    		}
    	//fill the odd lines of the even planes
    	for(int k=0; 2*k<this->upscaled.size[0]; ++k)
    		for(int j=0; 2*j+2<this->upscaled.size[1]; ++j)
    		{
    			OctaveFinder::PixelType * u = &this->upscaled(2*k, 2*j, 0),
    					* v = &this->upscaled(2*k, 2*j+1, 0),
    					* w = &this->upscaled(2*k, 2*j+2, 0);
    			for(int i=0; i<this->upscaled.size[2]; ++i)
    				*v++ = 0.5 * (*u++ + *w++);
    		}
    	//fill the odd planes
    	for(int k=0; 2*k+2<this->upscaled.size[0]; ++k)
    	{
    		OctaveFinder::PixelType * u = &this->upscaled(2*k, 0, 0),
						* v = &this->upscaled(2*k+1, 0, 0),
						* w = &this->upscaled(2*k+2, 0, 0);
			for(int i=0; i<this->upscaled.size[1]*this->upscaled.size[2]; ++i)
				*v++ = 0.5 * (*u++ + *w++);
    	}
	}
}
