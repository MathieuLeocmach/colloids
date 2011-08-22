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
		this->small = cv::Mat_<double>(nrows, ncols, 0.0);
		this->upscaled = cv::Mat_<double>(2*nrows, 2*ncols, 0.0);
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
		this->small = cv::Mat_<double>(1, ncols, 0.0);
		this->upscaled = cv::Mat_<double>(1, 2*ncols, 0.0);
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
    	for(int j=0; 2*j<this->upscaled.cols; ++j)
    		for(int i=0; 2*i<this->upscaled.rows; ++i)
    			this->upscaled(2*i, 2*j) = small(i,j);
    	for(int j=0; 2*j<this->upscaled.cols; ++j)
			for(int i=0; 2*i+1<this->upscaled.rows && i+1<small.rows; ++i)
				this->upscaled(2*i+1, 2*j) = 0.5*(small(i,j)+small(i+1,j));
    	for(int j=0; 2*j+1<this->upscaled.cols && j+1<small.cols; ++j)
			for(int i=0; 2*i<this->upscaled.rows; ++i)
				this->upscaled(2*i, 2*j+1) = 0.5*(small(i,j)+small(i, j+1));
    	for(int j=0; 2*j+1<this->upscaled.cols && j+1<small.cols; ++j)
			for(int i=0; 2*i+1<this->upscaled.rows && i+1<small.rows; ++i)
				this->upscaled(2*i+1, 2*j+1) = 0.25*(small(i,j)+small(i+1,j)+small(i, j+1)+small(i+1, j+1));

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

	/**
	 * \brief Locate centers with subpixel and subscale resolution
	 */
	void MultiscaleFinder::subpix(std::vector<Center2D> &centers) const
	{
		centers.clear();
		//reserve memory for the center container
    	size_t n_centers = 0;
    	for(size_t o=0; o<this->octaves.size(); ++o)
    		n_centers += this->octaves[o]->get_nb_centers();
    	centers.reserve(n_centers);
    	//subpixel resolution
    	for(size_t o=0; o<this->octaves.size(); ++o)
    	{
    		std::vector<Center2D> v;
    		this->octaves[o]->subpix(v);
    		//correct the seam between octaves in sizing precision
    		if(o>0)
				for(size_t p=0; p<v.size(); ++p)
					this->seam(v[p], o-1);
    		//transform scale coordinate in size coordinate
    		for(size_t c=0; c< v.size(); ++c)
    			this->octaves[o]->scale(v[c]);
			//stack up
			for(size_t c=0; c<v.size(); ++c)
			{
				v[c][0] *= pow(2.0, (int)(o)-1);
				v[c][1] *= pow(2.0, (int)(o)-1);
				v[c].r *= pow(2.0, (int)(o)-1);
				centers.push_back(v[c]);
			}
    	}
	}
	/**
	 * Efficient processing pipeline from the input image to the output centers
	 */
	void MultiscaleFinder::get_centers(const cv::Mat & input, std::vector<Center2D>& centers)
	{
		this->fill(input);
		this->initialize_binary();
		this->subpix(centers);
	}

    const cv::Vec3i MultiscaleFinder::previous_octave_coords(const Center2D &v) const
    {
    	const double n = this->get_n_layers();
    	cv::Vec3i vi;
		for(int u=0; u<2; ++u)
			vi[u] = (int)(v[u]*2+0.5);
		vi[2] = int(v.r + n + 0.5);
		return vi;
    }

    const cv::Mat_<double> MultiscaleFinder2D::downscale(const size_t &o)
	{
    	//second to last Gaussian layer of octave o-1 has a blurring radius two time larger than the original
		cv::Mat_<double> roi2 = small(
				cv::Range(0, this->octaves[o]->get_width()),
				cv::Range(0, this->octaves[o]->get_height())
		);
		const cv::Mat_<double> & a = this->octaves[o-1]->get_layersG(this->octaves[o-1]->get_n_layers());
		for(int j=0; j<roi2.cols && 2*j+1<a.cols; ++j)
			for(int i=0; i<roi2.rows && 2*i+1<a.rows; ++i)
				roi2(i,j) = (a(2*i, 2*j) + a(2*i+1, 2*j) + a(2*i, 2*j+1) + a(2*i+1, 2*j+1))/4.0;
		return roi2;
	}
    const cv::Mat_<double> MultiscaleFinder1D::downscale(const size_t &o)
	{
		//second to last Gaussian layer of octave o-1 has a blurring radius two time larger than the original
		cv::Mat_<double> roi2 = small(
				cv::Range(0, this->octaves[o]->get_width()),
				cv::Range(0, this->octaves[o]->get_height())
		);
		const cv::Mat_<double> & a = this->octaves[o-1]->get_layersG(this->octaves[o-1]->get_n_layers());
		for(int i=0; i<roi2.cols; ++i)
			roi2(0, i) = (a(0, 2*i) + a(0, 2*i+1))/2.0;
		return roi2;
	}
    void MultiscaleFinder2D::seam(Center2D &v, const size_t &o) const
    {
		if(v.r<1)
			v.r = this->octaves[o]->scale_subpix(this->previous_octave_coords(v)) - this->get_n_layers();
	}
    void MultiscaleFinder1D::seam(Center2D &v, const size_t &o) const
	{
    	//nothing
	}
}
