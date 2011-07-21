#include "octavefinder.hpp"
#include <boost/array.hpp>
#include <algorithm>
#include <stdexcept>


using namespace std;

namespace Colloids {

OctaveFinder::OctaveFinder(const int nrows, const int ncols, const int nbLayers, const double &preblur_radius):
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
	this->set_radius_preblur(preblur_radius);
}

OctaveFinder::~OctaveFinder()
{
    //dtor
}

void Colloids::OctaveFinder::set_radius_preblur(const double &k)
{
	this->preblur_radius = k;
	this->fill_iterative_radii(k);
}

void Colloids::OctaveFinder::fill(const cv::Mat &input)
{
	if(input.rows != this->get_width())
	{
		std::ostringstream os;
		os << "OctaveFinder::fill : the input's rows ("<<input.rows<<") must match the width of the finder ("<<this->get_width()<<")";
		throw std::invalid_argument(os.str().c_str());
	}
	if(input.cols != this->get_height())
	{
		std::ostringstream os;
		os << "OctaveFinder::fill : the input's cols ("<<input.cols<<") must match the height of the finder ("<<this->get_height()<<")";
		throw std::invalid_argument(os.str().c_str());
	}

	input.convertTo(this->layersG[0], this->layersG[0].type());
	//iterative Gaussian blur
	for(size_t i=0; i<this->layersG.size()-1; ++i)
		this->iterative_gaussian_filters[i].apply(
				this->layersG[i],
				this->layersG[i+1]
				);
	//difference of Gaussians
	for(size_t i=0; i<this->layers.size(); ++i)
		cv::subtract(this->layersG[i+1], this->layersG[i], this->layers[i]);
}

void Colloids::OctaveFinder::preblur_and_fill(const cv::Mat &input)
{
	if(input.rows != this->get_width())
	{
		std::ostringstream os;
		os << "OctaveFinder::preblur_and_fill : the input's rows ("<<input.rows<<") must match the width of the finder ("<<this->get_width()<<")";
		throw std::invalid_argument(os.str().c_str());
	}
	if(input.cols != this->get_height())
	{
		std::ostringstream os;
		os << "OctaveFinder::preblur_and_fill : the input's cols ("<<input.cols<<") must match the height of the finder ("<<this->get_height()<<")";
		throw std::invalid_argument(os.str().c_str());
	}
	cv::Mat_<double> blurred(input.size());
	cv::GaussianBlur(input, blurred, cv::Size(0,0), this->preblur_radius);
	this->fill(blurred);
}

/**
 * \brief Detect local minima of the scale space
 *
 * Uses the dynamic block algorythm by Neubeck and Van Gool
 * a. Neubeck and L. Van Gool, 18th International Conference On Pattern Recognition (ICPRʼ06) 850-855 (2006).
 */
void Colloids::OctaveFinder::initialize_binary(const double & max_ratio)
{
    //initialize
	this->centers_no_subpix.clear();
    for(size_t i = 0;i < this->binary.size();++i)
        this->binary[i].setTo(0);

    for(size_t k = 1;k < (size_t)(((this->layers.size() - 1)));k += 2)
            for(size_t j = 1;j < (size_t)(((this->get_width() - 1)));j += 2)
                for(size_t i = 1;i < (size_t)(((this->get_height() - 1)));i += 2){
                    size_t mi = i;
                    size_t mj = j;
                    size_t mk = k;
                    //accept degenerated minima in the block, but select only one
                    for(size_t k2 = k;k2 < k + 2;++k2)
                        for(size_t j2 = j;j2 < j + 2;++j2)
                            for(size_t i2 = i;i2 < i + 2;++i2)
                                if(this->layers[k2](j2, i2) <= this->layers[mk](mj, mi)){
                                    mi = i2;
                                    mj = j2;
                                    mk = k2;
                                }



                    //maxima cannot be on the last layer or on image edges
                    if(mk > this->binary.size() || !((this->sizes[mk] <= mj) && (mj < (size_t)(((this->get_width() - this->sizes[mk])))) && (this->sizes[mk] <= mi) && (mi < (size_t)(((this->get_height() - this->sizes[mk]))))))
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




                    //remove the local minima that are edges (elongated objects)
                    if(*b){
                        //hessian matrix
                        const double hess[3] = {this->layers[mk](mj - 1, mi) - 2 * this->layers[mk](mj, mi) + this->layers[mk](mj + 1, mi), this->layers[mk](mj, mi - 1) - 2 * this->layers[mk](mj, mi) + this->layers[mk](mj, mi + 1), this->layers[mk](mj - 1, mi - 1) + this->layers[mk](mj + 1, mi + 1) - this->layers[mk](mj + 1, mi - 1) - this->layers[mk](mj - 1, mi + 1)};
                        //determinant of the Hessian, for the coefficient see
                        //H Bay, a Ess, T Tuytelaars, and L Vangool,
                        //Computer Vision and Image Understanding 110, 346-359 (2008)
                        const double detH = hess[0] * hess[1] - pow(hess[2], 2), ratio = pow(hess[0] + hess[1], 2) / (4.0 * hess[0] * hess[1]);
                        *b = !(detH < 0 || ratio > max_ratio);
                        if(*b){
                            cv::Vec3i c;
                            c[0] = mi;
                            c[1] = mj;
                            c[2] = mk;
                            this->centers_no_subpix.push_back(c);
                        }
                    }

                }



    }

    double Colloids::OctaveFinder::gaussianResponse(const size_t & i, const size_t & j, const double & scale) const
    {
        if(scale < 0)
            throw std::invalid_argument("Colloids::OctaveFinder::gaussianResponse: the scale must be positive.");

        size_t k = (size_t)(scale);
        if (k>=this->layersG.size())
        	k = this->layersG.size()-1;
		if((scale - k) * (scale - k) + 1 == 1)
			return this->layersG[k](i, j);
		const double sigma = this->get_iterative_radius(scale, (double)k);
		//opencv is NOT dealing right with ROI (even if boasting about it), so we do it by hand
		const int m = ((int)(sigma*4+0.5)*2 + 1)|1;
		vector<double> gx(m, 0.0);
		cv::Mat_<double> kernel = cv::getGaussianKernel(m, sigma, this->layersG[0].type());
		for(int x=0; x<m; ++x)
			for(int y=0; y<m; ++y)
				gx[x] += layersG[k](i-x+m/2,j-y+m/2) * kernel(y,0);
		double resp = 0.0;
		for(int x=0; x<m; ++x)
			resp += gx[x] * kernel(x,0);

        return resp;
    }

    cv::Vec4d Colloids::OctaveFinder::single_subpix(const cv::Vec3i & ci) const
    {
        const size_t i = ci[0], j = ci[1], k = ci[2];
        //look for the pedestal
        double pedestal = this->layers[k](j, i);
        for(size_t k2 = k - 1;k2 < k + 2;++k2)
            for(size_t j2 = j - this->sizes[k];j2 < j + this->sizes[k] + 1;++j2)
                for(size_t i2 = i - this->sizes[k];i2 < i + this->sizes[k] + 1;++i2){
                    const double val = this->layers[k2](j2, i2);
                    if(val < 0 && val > pedestal)
                        pedestal = val;

                }


        cv::Vec4d c(0.0);
        if(pedestal == this->layers[k](j, i))
            pedestal = 0.0;

        //define the ROI: 3 pixels in the scale axis, according to the scale in space
        for(size_t k2 = k - 1;k2 < k + 2;++k2)
            for(size_t j2 = j - this->sizes[k];j2 < j + this->sizes[k] + 1;++j2)
                for(size_t i2 = i - this->sizes[k];i2 < i + this->sizes[k] + 1;++i2){
                    const double val = this->layers[k2](j2, i2);
                    if(val < 0){
                        const double v = val - pedestal;
                        c[0] += v * i2;
                        c[1] += v * j2;
                        //c[2] += v * k2;
                        c[3] += v;
                    }
                }



        for(int u = 0;u < 2;++u)
            c[u] /= c[3];

        //scale is better defined if we consider only the central pixel at different scales
        //Compute intermediate variables to do a quadratic estimate of the derivative
        boost::array<double,8> sublayerG;
		for(int u = 0;u < 3;++u)
			sublayerG[u] = this->gaussianResponse(i, j, k - 0.5 + u);

        boost::array<double,5> a = {{
        		this->layers[k - 1](j, i),
        		sublayerG[1] - sublayerG[0],
        		this->layers[k](j, i),
        		sublayerG[2] - sublayerG[1],
        		this->layers[k + 1](j, i)}};
        //Apply newton's method using quadratic estimate of the derivative
        c[2] = k - (-a[4] + 8*a[3] - 8*a[1] + a[0])/6.0 /(a[4]-2*a[2]+a[0]);
        //saturate the results
        if(c[2]>k+0.5)
        	c[2]= k + 0.5;
        if(c[2]<k-0.5)
        	c[2]= k- 0.5;
        //second round of newton's method
        if(c[2]>=2)
        {
        	if(c[2]<k)
        		c[2]= k- 0.5;
			for(size_t u = 0;u < sublayerG.size(); ++u)
				sublayerG[u] = this->gaussianResponse(i, j, c[2] - 1 + 0.5*u);
			for(size_t u =0; u<a.size();++u)
				a[u] = sublayerG[u+2] - sublayerG[u];
			c[2] -= (-a[4] + 8*a[3] - 8*a[1] + a[0])/6.0 /(a[4]-2*a[2]+a[0]);
        }
        if(c[2]<k-0.5)
			c[2]= k- 0.5;
        if(c[2]>k+0.5)
			c[2]= k + 0.5;
        return c;
    }

    std::vector<cv::Vec4d> Colloids::OctaveFinder::subpix() const
    {
        std::vector<cv::Vec4d> centers;
        centers.reserve(this->centers_no_subpix.size());
        //subpixel resolution in pixel units
        for(std::list<cv::Vec3i>::const_iterator ci = this->centers_no_subpix.begin();ci != this->centers_no_subpix.end();++ci)
            centers.push_back(this->single_subpix(*ci));

        return centers;
    }
    /**
	 * \brief Convert scale to real size
	 *
	 * Identify the maximum response of a difference of Gaussians for a perfect circle.
	 */
    void Colloids::OctaveFinder::scale(std::vector<cv::Vec4d> & centers) const
    {
        //convert scale to real size
        const double n = this->get_n_layers(), prefactor = 2.0 * this->preblur_radius * sqrt(log(2.0) / n / (pow(2.0, 2.0 / n) - 1));
        for(size_t c = 0;c < centers.size();++c)
            centers[c][2] = prefactor * pow(2.0, (centers[c][2] + 1) / n);

    }
    std::vector<cv::Vec4d> Colloids::OctaveFinder::operator ()(const cv::Mat & input, const bool preblur)
    {
        if(preblur)
            this->preblur_and_fill(input);

        else
            this->fill(input);

        this->initialize_binary();
        std::vector<cv::Vec4d> centers = this->subpix();
        this->scale(centers);
        return centers;
    }
    const double OctaveFinder::get_iterative_radius(const double & larger, const double & smaller) const
    {
    	return this->preblur_radius * sqrt(pow(2.0, 2.0*larger/this->get_n_layers()) - pow(2.0, 2.0*smaller/this->get_n_layers()));
    }

void Colloids::OctaveFinder::fill_iterative_radii(const double & k)
{
		//target blurring radii
        vector<double> sigmas(this->get_n_layers()+3);
        for(size_t i=0; i<sigmas.size(); ++i)
        	sigmas[i] = k * pow(2, i/double(this->get_n_layers()));
        //corresponding blob sizes
        vector<size_t>::iterator si = this->sizes.begin();
        for(vector<double>::const_iterator sig=sigmas.begin(); sig!=sigmas.end(); sig++)
        	*si++ = static_cast<size_t>((*sig * sqrt(2)) + 0.5);
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



}; //namespace





