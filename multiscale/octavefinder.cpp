#include "octavefinder.hpp"
#include <algorithm>

using namespace cv;
using namespace std;
using namespace Colloids;

OctaveFinder::OctaveFinder(const int nrows, const int ncols, const int nbLayers, const double &k):
		iterative_radii(nbLayers+2), sizes(nbLayers+3)
{
    this->layersG.reserve(nbLayers+3);
    for (int i = 0; i<nbLayers+3; ++i)
    	this->layersG.push_back(Mat_<double>(nrows, ncols));
    this->layers.reserve(nbLayers+2);
	for (int i = 0; i<nbLayers+2; ++i)
		this->layers.push_back(Mat_<double>(nrows, ncols));
	this->binary.reserve(nbLayers+2);
	for (int i = 0; i<nbLayers+2; ++i)
		this->binary.push_back(Mat_<bool>(nrows, ncols));
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
	//iterative blur
	for(size_t i=0; i<this->layersG.size()-1; ++i)
		cv::GaussianBlur(
				this->layersG[i],
				this->layersG[i+1],
				cv::Size(0,0),
				this->iterative_radii[i]
				);
}

void Colloids::OctaveFinder::preblur_and_fill(const cv::Mat &input)
{
	cv::GaussianBlur(input, this->layersG[0], cv::Size(0,0), this->k);
	this->fill(this->layersG[0]);
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

}



