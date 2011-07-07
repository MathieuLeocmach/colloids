#include "octavefinder.hpp"
#include <algorithm>
#include <boost/array.hpp>

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

/**
 * \brief Detect local minima of the scale space
 *
 * Uses the dynamic block algorythm by Neubeck and Van Gool
 * a. Neubeck and L. Van Gool, 18th International Conference On Pattern Recognition (ICPRʼ06) 850-855 (2006).
 */
void Colloids::OctaveFinder::local_minima()
{
    //initialize
    for(size_t i = 0;i < this->binary.size();++i)
        this->binary[i].setTo(0);

    for(size_t k = 1;k < this->layers.size() - 1;k += 2)
        for(size_t j = 1;j < this->get_width() - 1;j += 2)
            for(size_t i = 1;i < this->get_height() - 1;i += 2)
            {
                size_t mi = i;
                size_t mj = j;
                size_t mk = k;
                //accept degenerated minima in the block, but select only one
                for(size_t k2 = k;k2 < k + 2;++k2)
                    for(size_t j2 = j;j2 < j + 2;++j2)
                        for(size_t i2 = i;i2 < i + 2;++i2)
                            if(this->layers[k2](j2, i2) <= this->layers[mk](mj, mi))
                            {
                                mi = i2;
                                mj = j2;
                                mk = k2;
                            }



                //maxima cannot be on the last layer
                if(mk > this->binary.size())
                    continue;

                bool *b = &this->binary[mk - 1](mj, mi);
                //consider only negative minima
                //with a value that is actually different from zero
                *b = (this->layers[mk](mj, mi) < 0) && (1 + pow(this->layers[mk](mj, mi), 2) > 1);
                //remove the minima if one of its neighbours outside the block has lower value
                for(size_t k2 = mk - 1;k2 < mk + 2 && *b;++k2)
                    for(size_t j2 = mj - 1;j2 < mj + 2 && *b;++j2)
                        for(size_t i2 = mi - 1;i2 < mi + 2 && *b;++i2)
                            if(k2 < k || j2 < j || i2 < i || k2 > k + 1 || j2 > j + 1 || i2 > i + 1)
                                *b = this->layers[mk](mj, mi) <= this->layers[k2](j2, i2);
            }


    //remove margins
    for(size_t i = 0;i < this->binary.size();++i)
    {
        this->binary[i].row(0).setTo(0);
        this->binary[i].row(this->get_width() - 1).setTo(0);
        this->binary[i].col(0).setTo(0);
        this->binary[i].col(this->get_height() - 1).setTo(0);
    }
}

void Colloids::OctaveFinder::remove_edges(const double & max_ratio)
{
	 for(size_t k = 1; k < this->layers.size() - 1; ++k)
		 for(size_t j = 1;j < this->get_width() - 1; ++j)
			 for(size_t i = 1;i < this->get_height() - 1; ++i)
			 {
				 bool *b = &this->binary[k - 1](j, i);
				 if(*b)
				 {
					 const double hess[3] = {
							 this->layers[k](j-1, i) - 2 * this->layers[k](j, i) + this->layers[k](j+1, i),
							 this->layers[k](j, i-1) - 2 * this->layers[k](j, i) + this->layers[k](j, i+1),
							 this->layers[k](j-1, i-1) + this->layers[k](j+1, i+1) - this->layers[k](j+1, i-1) - this->layers[k](j-1, i+1)
					 };
					 //determinant of the Hessian, for the coefficient see
					 //H Bay, a Ess, T Tuytelaars, and L Vangool,
					 //Computer Vision and Image Understanding 110, 346-359 (2008)
					 const double detH = hess[0]*hess[1] - pow(hess[2], 2),
							 ratio = pow(hess[0]+hess[1], 2) / (4.0 * hess[0]*hess[1]);
					 *b = !(detH<0 || ratio > max_ratio);
				 }
			 }
}

void Colloids::OctaveFinder::initialize_binary()
{
    local_minima();
    remove_edges();
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









