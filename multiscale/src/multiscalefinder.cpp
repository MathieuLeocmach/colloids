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
	}
	MultiscaleFinder3D::MultiscaleFinder3D(const int nplanes, const int nrows, const int ncols, const int nbLayers, const double &preblur_radius, bool incore):
			deconv(false)
	{
		//must not fail if the image is too small, construct the 0th octave anyway
		const int s = (nplanes<10 || nrows<10 || ncols<10)?1:(log(min(nplanes, min(nrows, ncols))/5)/log(2));
		this->octaves.reserve((size_t)s);
		this->octaves.push_back(new OctaveFinder3D(2*nplanes, 2*nrows, 2*ncols, nbLayers, preblur_radius, false));
		int opl = nplanes, ocr = nrows, occ = ncols;
		while(opl >=10 && ocr >= 10 && occ >= 10)
		{
			this->octaves.push_back(new OctaveFinder3D(opl, ocr, occ, nbLayers, preblur_radius, incore));
			opl /= 2;
			ocr /= 2;
			occ /= 2;
		}
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
    void MultiscaleFinder3D::set_ZXratio(const double &ratio)
    {
    	for(size_t o=0; o<this->octaves.size(); ++o)
    		dynamic_cast<OctaveFinder3D*>(this->octaves[o])->set_ZXratio(ratio);
    }
    void MultiscaleFinder3D::set_halfZpreblur(bool value)
    {
    	if(this->octaves.size()>1)
    		dynamic_cast<OctaveFinder3D*>(this->octaves[1])->set_halfZpreblur(value);
    }
    void MultiscaleFinder3D::load_deconv_kernel(const std::vector<PixelType> &kernel)
    {
    	if(kernel.size() != this->get_height()/2+1)
    		throw std::invalid_argument("The kernel size must match the last dimension of the finder");
    	this->deconvKernel = kernel;
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
    	if(this->use_Octave0())
    	{
			//half preblur and upscale the input to fill the first octave
			Image upscaled = this->upscale(input);
			this->octaves[0]->fill(upscaled);
    	}
    	if(this->octaves.size()>1)
			//Octave 1 corresponds to the size of the input image.
			//To avoid errors in the upsampling+downsampling process, we use the input directly
			this->octaves[1]->preblur_and_fill(input);
    	//For higher octaves we use the second to last layer of the previous octave, downsampled
    	for(size_t o=2; o<this->octaves.size(); ++o)
    	{
    		Image small = this->downscale(o);
    		this->octaves[o]->fill(small);
    	}
    }
    /**
     * \brief Locate centers with pixel resolutions and scale resolution
     */
	void MultiscaleFinder::initialize_binary()
	{
		const size_t o0 = this->use_Octave0()?0:1;
		//initialize binary for each octave
    	for(size_t o=o0; o<this->octaves.size(); ++o)
			this->octaves[o]->initialize_binary();
    	//Remove pixel centers that exist in consecutive octaves
    	/*for(size_t o=0; o<this->octaves.size()-1; ++o)
    		this->octaves[o]->seam_binary(*this->octaves[o+1]);*/
	}

    MultiscaleFinder::Image MultiscaleFinder2D::downscale(const size_t &o) const
	{
    	//second to last Gaussian layer of octave o-1 has a blurring radius two time larger than the original
    	Image roi2(this->octaves[o]->get_width(), this->octaves[o]->get_height());
		const Image & a = this->octaves[o-1]->get_layersG(this->octaves[o-1]->get_n_layers());
		for(int j=0; j<roi2.cols && 2*j+1<a.cols; ++j)
			for(int i=0; i<roi2.rows && 2*i+1<a.rows; ++i)
				roi2(i,j) = (a(2*i, 2*j) + a(2*i+1, 2*j) + a(2*i, 2*j+1) + a(2*i+1, 2*j+1))/4.0;
		return roi2;
	}
    MultiscaleFinder::Image MultiscaleFinder1D::downscale(const size_t &o) const
	{
		//second to last Gaussian layer of octave o-1 has a blurring radius two time larger than the original
    	Image roi2(this->octaves[o]->get_width(), this->octaves[o]->get_height());
		const Image & a = this->octaves[o-1]->get_layersG(this->octaves[o-1]->get_n_layers());
		for(int i=0; i<roi2.cols; ++i)
			roi2(0, i) = (a(0, 2*i) + a(0, 2*i+1))/2.0;
		return roi2;
	}
    MultiscaleFinder::Image MultiscaleFinder3D::downscale(const size_t &o) const
	{
		//second to last Gaussian layer of octave o-1 has a blurring radius two time larger than the original
    	int dims[3] = {
    			dynamic_cast<OctaveFinder3D*>(this->octaves[o])->get_depth(),
    			this->octaves[o]->get_width(),
    			this->octaves[o]->get_height()};
		Image roi2(3, dims, (OctaveFinder::PixelType)0);
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
    MultiscaleFinder::Image MultiscaleFinder2D::upscale(const cv::Mat &input) const
    {
    	Image halfblured(input);
    	cv::GaussianBlur(halfblured, halfblured, cv::Size(0,0), this->get_radius_preblur()/2.0);
    	Image upscaled(
    			input.size[0]*2,
    			input.size[1]*2,
    			(OctaveFinder::PixelType)0);
    	//fill the even lines
    	for(int j=0; 2*j<upscaled.size[0]; ++j)
		{
    		//convert the unknown input type to PixelType
    		Image input_row;
    		halfblured.row(j).convertTo(input_row, input_row.type());
    		//copy to the even pixels
    		PixelType * u = &upscaled(2*j, 0);
			const PixelType * s = &input_row(0,0);
			for(int i=0; 2*i<upscaled.size[1]; ++i)
			{
				*u++ = *s++;
				u++;
			}
			//interpolate the odd pixels
			u = &upscaled(2*j, 1);
			for(int i=0; 2*i+2<upscaled.size[1]; ++i)
			{
				*u = 0.5 * (*(u-1) + *(u+1));
				u++;
				u++;
			}
		}
    	//fill the odd lines
    	for(int j=0; 2*j+2<upscaled.size[0]; ++j)
		{
			OctaveFinder::PixelType * u = &upscaled(2*j, 0),
					* v = &upscaled(2*j+1, 0),
					* w = &upscaled(2*j+2, 0);
			for(int i=0; i<upscaled.size[1]; ++i)
				*v++ = 0.5 * (*u++ + *w++);
		}
    	return upscaled;
    }
    MultiscaleFinder::Image MultiscaleFinder1D::upscale(const cv::Mat &input) const
	{
		Image upscaled(1, input.size[1]*2, (OctaveFinder::PixelType)0);
		//convert the unknown input type to PixelType
		Image input_row;
		input.convertTo(input_row, input_row.type());
		cv::GaussianBlur(input_row, input_row, cv::Size(0,0), this->get_radius_preblur()/2.0);
		//copy to the even pixels
		PixelType * u = &upscaled(0, 0);
		const PixelType * s = &input_row(0,0);
		for(int i=0; i<input_row.size[1]; ++i)
		{
			*u++ = *s++;
			u++;
		}
		//interpolate the odd pixels
		u = &upscaled(0, 1);
		for(int i=0; i+1<input_row.size[1]; ++i)
		{
			*u = 0.5 * (*(u-1) + *(u+1));
			u++;
			u++;
		}
		return upscaled;
	}
    MultiscaleFinder::Image MultiscaleFinder3D::upscale(const cv::Mat &input) const
	{
    	int dims[3] = {
			dynamic_cast<OctaveFinder3D*>(this->octaves[0])->get_depth(),
			this->octaves[0]->get_width(),
			this->octaves[0]->get_height()};
    	Image halfblurred;
    	input.convertTo(halfblurred, halfblurred.type());
    	inplace_blur3D(halfblurred, this->get_radius_preblur()/2.0, this->get_ZXratio());
    	Image upscaled(3, dims, (OctaveFinder::PixelType)0);
    	cv::Mat input2D(input.size[0]*input.size[1], input.size[2], halfblurred.type(), halfblurred.data);
    	//fill the even lines of the even planes
    	for(int k=0; k<input.size[0]; ++k)
    		for(int j=0; j<input.size[1]; ++j)
    		{
    			//convert the unknown input type to PixelType
				Image input_row;
				input2D.row(k*input.size[1]+j).convertTo(input_row, input_row.type());
				//copy to the even pixels
				PixelType * u = &upscaled(2*k, 2*j, 0);
				const PixelType * s = &input_row(0,0,0);
				for(int i=0; 2*i<upscaled.size[2]; ++i)
				{
					*u++ = *s++;
					u++;
				}
				//interpolate the odd pixels
				for(int i=0; 2*i+2<upscaled.size[2]; ++i)
				{
					*u = 0.5 * (*(u-1) + *(u+1));
					u++;
				}
    		}
    	//interpolate the odd lines of the even planes
    	for(int k=0; 2*k<upscaled.size[0]; ++k)
    		for(int j=0; 2*j+2<upscaled.size[1]; ++j)
    		{
    			OctaveFinder::PixelType * u = &upscaled(2*k, 2*j, 0),
    					* v = &upscaled(2*k, 2*j+1, 0),
    					* w = &upscaled(2*k, 2*j+2, 0);
    			for(int i=0; i<upscaled.size[2]; ++i)
    				*v++ = 0.5 * (*u++ + *w++);
    		}
    	//interpolate the odd planes
    	for(int k=0; 2*k+2<upscaled.size[0]; ++k)
    	{
    		OctaveFinder::PixelType * u = &upscaled(2*k, 0, 0),
						* v = &upscaled(2*k+1, 0, 0),
						* w = &upscaled(2*k+2, 0, 0);
			for(int i=0; i<upscaled.size[1]*upscaled.size[2]; ++i)
				*v++ = 0.5 * (*u++ + *w++);
    	}
    	return upscaled;
	}
    /**
     * \brief Correct radii by solving inter-particle coupling
     */
    void MultiscaleFinder3D::global_scale2radius(std::vector<Center3D > &centers) const
    {
    	//bonds and distances

    	//solve system 1

    	//solve system 2

    	//results

    }
}
