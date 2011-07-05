#include "octavefinder.hpp"
#include <algorithm>

using namespace std;
using namespace Colloids;

OctaveFinder::OctaveFinder(const int nrows, const int ncols, const int nbLayers, const double &k):
		iterative_radii(nbLayers+2), iterative_gaussian_filters(nbLayers+2), sizes(nbLayers+3)
{
    this->layersG.reserve(nbLayers+3);
    for (int i = 0; i<nbLayers+3; ++i)
    	this->layersG.push_back(cv::Mat_<double>(nrows, ncols));
    this->layers.reserve(nbLayers+2);
	for (int i = 0; i<nbLayers+2; ++i)
		this->layers.push_back(cv::Mat_<double>(nrows, ncols));
	this->binary.reserve(nbLayers+2);
	for (int i = 0; i<nbLayers; ++i)
		this->binary.push_back(cv::Mat_<bool>::zeros(nrows, ncols));
	this->set_radius_preblur(k);
}

OctaveFinder::~OctaveFinder()
{
    //dtor
}

void Colloids::OctaveFinder::set_radius_preblur(const double &k)
{
	this->k = k;
	this->fill_iterative_radii(k);
}

void Colloids::OctaveFinder::fill(const cv::Mat &input)
{
	input.convertTo(this->layersG[0], this->layersG[0].type());
	//iterative Gaussian blur
	for(size_t i=0; i<this->layersG.size()-1; ++i)
		this->iterative_gaussian_filters[i].apply(
				this->layersG[i],
				this->layersG[i+1]
				);
	//difference of Gaussians
	for(size_t i=0; i<this->layers.size(); ++i)
		this->layers[i] = this->layersG[i+1] - this->layersG[i];
}

void Colloids::OctaveFinder::preblur_and_fill(const cv::Mat &input)
{
	cv::GaussianBlur(input, this->layersG[0], cv::Size(0,0), this->k);
	this->fill(this->layersG[0]);
}

void Colloids::OctaveFinder::initialize_binary()
{
	for(size_t i=0; i<this->binary.size(); ++i)
		this->binary[i].setTo(0);

	for(int l=0; l < this->get_n_layers(); ++l)
	{
		//Iterate through blocks of 3^3 pixels
		std::vector<cv::Mat_<double> > arrays;
		arrays.reserve(27);
		for(int k = l; k<l+3; ++k)
			for(int j=0; j<3; ++j)
				for(int i=0; i<3; ++i)
					arrays.push_back(this->layers[k](
							cv::Range(j, this->get_width()-2+j),
							cv::Range(i, this->get_height()-2+i)
							));
		//move the central pixel at the end
		std::swap(arrays[arrays.size()/2], arrays.back());
		std::vector<cv::MatConstIterator_<double> > its;
		std::transform(
				arrays.begin(), arrays.end(), std::back_inserter(its),
				std::mem_fun_ref<cv::MatConstIterator_<double> >(&cv::Mat_<double>::begin)
		);
		//ROI of the binary layer
		cv::Mat_<bool> broi = this->binary[l](
				cv::Range(1, this->get_width()-1),
				cv::Range(1, this->get_height()-1));
		for(cv::Mat_<bool>::iterator b_it= broi.begin(); b_it!=broi.end(); ++b_it)
		{
			//only negative pixels
			*b_it = *(its.back())<0;
			//non minimum suppression
			size_t i=0;
			while(*b_it && i<its.size()-1)
				*b_it = *(its.back()) < *(its[i++]);
			//increment the pointers
			for(std::vector<cv::MatConstIterator_<double> >::iterator it=its.begin(); it!=its.end(); ++it)
				(*it)++;
		}
	}
}

void Colloids::OctaveFinder::fill_iterative_radii(const double & k)
{
		//target blurring radii
        vector<double> sigmas(this->get_n_layers()+3);
        for(size_t i=0; i<sigmas.size(); ++i)
        	sigmas[i] = k * pow(2, i/double(this->get_n_layers()));
        //corresponding blob sizes
        vector<int>::iterator si = this->sizes.begin();
        for(vector<double>::const_iterator sig=sigmas.begin(); sig!=sigmas.end(); sig++)
        	*si++ = static_cast<int>((*sig * sqrt(2)) + 0.5);
        //iterative blurring radii
        transform(sigmas.begin(), sigmas.end(), sigmas.begin(), sigmas.begin(), multiplies<double>());
        for(size_t i=0; i<this->iterative_radii.size(); ++i)
        	this->iterative_radii[i] = sqrt(sigmas[i+1] - sigmas[i]);
        //iterative gaussian filters
        for(size_t i=0; i<this->iterative_radii.size(); ++i)
        	this->iterative_gaussian_filters[i] = *cv::createGaussianFilter
        	(
        			this->layersG[0].type(),
        			cv::Size(0,0),
        			this->iterative_radii[i]
        	);
}





