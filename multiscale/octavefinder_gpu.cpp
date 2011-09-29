#include "octavefinder_gpu.hpp"
#include <boost/array.hpp>
#include <algorithm>
#include <stdexcept>


using namespace std;

namespace Colloids {

OctaveFinder_GPU::OctaveFinder_GPU(const int nrows, const int ncols, const int nbLayers, const double &preblur_radius):
		OctaveFinder(nrows, ncols, nbLayers, preblur_radius)
{

}

OctaveFinder_GPU::~OctaveFinder_GPU()
{
    //dtor
}

void Colloids::OctaveFinder_GPU::set_radius_preblur(const double &k)
{
	this->preblur_radius = k;
	this->fill_iterative_radii();
	this->initialize_gpu();
}

void Colloids::OctaveFinder_GPU::fill(const cv::Mat &input)
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

	input.convertTo(this->layersG.front(), this->layersG.front().type());
	if(this->use_gpu())
	{
		str.enqueueUpload(this->layersG.front(), this->gpu_layersG.front());
		//iterative Gaussian blur
		for(size_t i=0; i<gpu_layersG.size()-1; ++i)
		{
			this->gpu_iterative_gaussian_filters[i]->apply(this->gpu_layersG[i], this->gpu_layersG[i+1], cv::Rect(0,0,-1,-1), str);
			for(int j=1; j<this->n_iter[i]; ++j)
			{
				str.enqueueCopy(this->gpu_layersG[i+1], this->tmp);
				this->gpu_iterative_gaussian_filters[i]->apply(this->tmp, this->gpu_layersG[i+1], cv::Rect(0,0,-1,-1), str);
			}
		}

		for(size_t i=0; i<gpu_layersG.size(); ++i)
			str.enqueueDownload(this->gpu_layersG[i], this->layersG[i]);

		//difference of Gaussians
		for(size_t i=0; i<gpu_layers.size(); ++i)
			cv::gpu::subtract(this->gpu_layersG[i+1], this->gpu_layersG[i], this->gpu_layers[i]);//, str);
		for(size_t i=0; i<gpu_layers.size(); ++i)
			str.enqueueDownload(this->gpu_layers[i], this->layers[i]);
	}
	else
		OctaveFinder::fill(input);
}

/**
 * \brief Detect local minima of the scale space
 *
 * Uses the dynamic block algorythm by Neubeck and Van Gool
 * a. Neubeck and L. Van Gool, 18th International Conference On Pattern Recognition (ICPRʼ06) 850-855 (2006).
 */
void Colloids::OctaveFinder_GPU::initialize_binary(const double & max_ratio)
{
    //initialize
	this->centers_no_subpix.clear();
	if(this->use_gpu())
		this->str.waitForCompletion();
	/*if(this->use_gpu())
	{
		cv::gpu::Stream str;
		cv::gpu::GpuMat a(this->binary.front().size(), this->binary.front().type()),
				b(this->binary.front().size(), this->binary.front().type()),
				gpu_bin(this->binary.front().size(), this->binary.front().type());
		for(size_t k = 0;k < this->binary.size();++k)
		{
			if(2 * this->sizes[k+1] >= this->get_height())
			{
				this->binary[k].setTo(0);
				continue;
			}
			const int from = this->sizes[k+1], tow = this->get_width() - this->sizes[k+1], toh = this->get_height() - this->sizes[k+1];
			cv::gpu::GpuMat roi = this->gpu_layers[k+1](cv::Rect(from, from, tow, toh));
			cv::gpu::GpuMat roibin = b(cv::Rect(from, from, tow, toh));
			str.enqueueMemSet(gpu_bin, 1);
			for(size_t mk=k; mk<k+3;++mk)
				for(size_t mj=-1; mj<2;++mj)
					for(size_t mi=-1;mi<2;++mi)
					{
						cv::gpu::compare(
								roi,
								this->gpu_layers[mk](cv::Rect(from+mj, from+mi, tow+mj, toh+mi)),
								roibin, cv::CMP_LE, str);
						cv::gpu::bitwise_and(gpu_bin, b, this->tmp, cv::gpu::GpuMat(), str);
						str.enqueueCopy(this->tmp, gpu_bin);
					}
			str.enqueueDownload(gpu_bin, this->binary[k]);
			str.waitForCompletion();
		}
		for(size_t k = 1;k < (size_t)(((this->layers.size() - 1)));k++)
			for(size_t j = this->sizes[k];j < (size_t)(((this->get_width() - this->sizes[k])));j++)
			{
				bool *b = &this->binary[k - 1](j, this->sizes[k]);
				float *v = &this->layers[k](j, this->sizes[k]);
				for(size_t i = this->sizes[k];i < (size_t)(((this->get_height() -this->sizes[k])));i++)
				{
					if(*b && (*v<0) && (*v * *v +1 > 1))
					{
						const double hess[3] = {this->layers[k](j - 1, i) - 2 * this->layers[k](j, i) + this->layers[k](j + 1, i), this->layers[k](j, i - 1) - 2 * this->layers[k](j, i) + this->layers[k](j, i + 1), this->layers[k](j - 1, i - 1) + this->layers[k](j + 1, i + 1) - this->layers[k](j + 1, i - 1) - this->layers[k](j - 1, i + 1)};
						//determinant of the Hessian, for the coefficient see
						//H Bay, a Ess, T Tuytelaars, and L Vangool,
						//Computer Vision and Image Understanding 110, 346-359 (2008)
						const double detH = hess[0] * hess[1] - pow(hess[2], 2), ratio = pow(hess[0] + hess[1], 2) / (4.0 * hess[0] * hess[1]);
						*b = !((detH < 0 && 1+detH*detH > 1) || ratio > max_ratio);
						if(*b){
							cv::Vec3i c;
							c[0] = i;
							c[1] = j;
							c[2] = k;
							this->centers_no_subpix.push_back(c);
						}
					}
					b++;
					v++;
				}
			}
	}
	else*/
		OctaveFinder::initialize_binary(max_ratio);
}
void Colloids::OctaveFinder_GPU::initialize_gpu()
{
	this->gpu_layers.clear();
	this->gpu_layersG.clear();
	this->gpu_iterative_gaussian_filters.clear();
	this->n_iter.clear();
	//const int mMax = ((int)(this->iterative_radii.back()*4+0.5)*2 + 1)|1;
	//won't use GPU for small images
	if(this->layers.front().rows>32 && this->layers.front().cols>32)
	{
		this->gpu_layers.reserve(this->layers.size());
		for(size_t i=0; i<this->layers.size();++i)
			this->gpu_layers.push_back(cv::gpu::GpuMat(this->layers.front().size(), this->layers.front().type()));
		this->gpu_layersG.reserve(this->layersG.size());
		for(size_t i=0; i<this->layersG.size();++i)
			this->gpu_layersG.push_back(cv::gpu::GpuMat(this->layersG.front().size(), this->layersG.front().type()));
		this->gpu_iterative_gaussian_filters.resize(this->iterative_radii.size());
		this->n_iter.resize(this->iterative_radii.size());
		//iterative gaussian filters
		//GPU filters cannot have a kernel larger than 15, i.e. sigma<=1.625
		for(size_t i=0; i<this->iterative_radii.size(); ++i)
		{
			this->n_iter[i] = pow(this->iterative_radii[i]/1.625, 2)+1;
			const int ksize = ((int)(this->iterative_radii[i]/sqrt(this->n_iter[i])*4+0.5)*2 + 1)|1;
			this->gpu_iterative_gaussian_filters[i] = cv::gpu::createGaussianFilter_GPU
			(
					this->layersG[0].type(),
					cv::Size(ksize,ksize),
					this->iterative_radii[i]/sqrt(this->n_iter[i])
			);
		}
	}
}


}; //namespace





