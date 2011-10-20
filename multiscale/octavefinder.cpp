#include "octavefinder.hpp"
#include <boost/array.hpp>
#include <algorithm>
#include <numeric>
#include <stdexcept>


using namespace std;

namespace Colloids {

std::map<size_t, cv::Mat_<double> > OctaveFinder::kernels;

const cv::Mat_<double>& OctaveFinder::get_kernel(const double &sigma)
{
	const int m = ((int)(sigma*4+0.5)*2 + 1)|1;
	//grab the kernel if it already exists within 1% precision
	std::map<size_t, cv::Mat_<double> >::iterator ker = kernels.find(100*sigma);
	if(ker==kernels.end())
		ker = kernels.insert(std::make_pair(
				100*sigma,
				cv::getGaussianKernel(m, sigma, CV_32F)
		)).first;
	return ker->second;
}

OctaveFinder::OctaveFinder(const int nrows, const int ncols, const int nbLayers, const double &preblur_radius):
		iterative_radii(nbLayers+2), iterative_gaussian_filters(nbLayers+2), sizes(nbLayers+3)
{
    this->layersG.reserve(nbLayers+3);
    for (int i = 0; i<nbLayers+3; ++i)
    	this->layersG.push_back(Image(nrows, ncols, (PixelType)0));
    this->layers.reserve(nbLayers+2);
	for (int i = 0; i<nbLayers+2; ++i)
		this->layers.push_back(Image(nrows, ncols, (PixelType)0));
	this->binary.reserve(nbLayers);
	for (int i = 0; i<nbLayers; ++i)
		this->binary.push_back(cv::Mat_<bool>(nrows, ncols, (PixelType)0));
	this->set_radius_preblur(preblur_radius);
}

OctaveFinder3D::OctaveFinder3D(const int nplanes, const int nrows, const int ncols, const int nbLayers, const double &preblur_radius) :
	OctaveFinder(0, 0, nbLayers, preblur_radius), iterative_Zgaussian_filters(nbLayers+2)
{
	//random file name in the working directory
	this->path.reserve(30);
	this->path.push_back('_');
	this->path.push_back('_');
	for(int i=0; i<28;++i)
		this->path.push_back('a'+rand()%('Z'-'a'));
	//create a memory mapped file to contain the images data
	boost::iostreams::mapped_file_params params(this->path);
	const size_t nbpixels =  nplanes * nrows * ncols;
	params.new_file_size = nbpixels * sizeof(PixelType) * (2 * nbLayers + 5) + nbpixels * nbLayers;
	params.flags = boost::iostreams::mapped_file::readwrite;
	this->file.open(params);
	//create the images inside the memory mapped file
	int dims[3] = {nplanes, nrows, ncols};
	for (int i = 0; i<nbLayers+3; ++i)
		this->layersG[i] = cv::Mat(
				3, dims, layersG[i].type(),
				(void*)(this->file.data() + i * nbpixels * sizeof(PixelType))
				);
	for (int i = 0; i<nbLayers+2; ++i)
		this->layers[i] = cv::Mat(
				3, dims, layers[i].type(),
				(void*)(this->file.data() + (i + layersG.size()) * nbpixels * sizeof(PixelType))
				);
	for (int i = 0; i<nbLayers; ++i)
		this->binary[i] = cv::Mat(
				3, dims, binary[i].type(),
				(void*)(this->file.data() + (2 * nbLayers + 5) * nbpixels * sizeof(PixelType) + i * nbpixels)
				);
	//create supplementary image headers for layersG data of shape (nplanes, nrows*ncols)
	this->layersG2D.reserve(nbLayers+3);
	for (int i = 0; i<nbLayers+3; ++i)
		this->layersG2D.push_back(cv::Mat(
				nplanes, nrows*ncols, layersG[0].type(),
				(void*)this->layersG[i].data
				));
	this->set_radius_preblur(preblur_radius);
}

OctaveFinder::~OctaveFinder()
{
    //dtor
}

OctaveFinder3D::~OctaveFinder3D()
{
	remove(this->path.c_str());
}

void Colloids::OctaveFinder::set_radius_preblur(const double &k)
{
	this->preblur_radius = k;
	this->fill_iterative_radii();
	this->preblur_filter = cv::createGaussianFilter(
			this->layersG.front().type(),
			cv::Size(0,0),
			this->preblur_radius
	);
}

/**
 * Convert an image of any kind to PixelType, but already pre-blurred, and process it to fill internal buffers
 */
void Colloids::OctaveFinder::fill(const cv::Mat &input)
{
	Image temp;
	input.convertTo(temp, temp.type());
	this->fill(temp);
}
/**
 * Process the allready pre-blurred input in place to fill internal buffers
 */
void Colloids::OctaveFinder::fill(Image &input)
{
	if(input.dims != this->layersG.front().dims)
	{
		std::ostringstream os;
		os << "OctaveFinder::fill : the input's dimension ("<<input.dims<<") must match the dimension of the finder ("<<this->layersG.front().dims<<")";
		throw std::invalid_argument(os.str().c_str());
	}
	for(int d=0; d<input.dims; ++d)
		if(input.size[d] != this->layersG.front().size[d])
		{
			std::ostringstream os;
			os << "OctaveFinder::fill : the input's "<<d<< "th dimension ("<<input.size[d]<<") must match the one of the finder ("<<this->layersG.front().size[d]<<")";
			throw std::invalid_argument(os.str().c_str());
		}
	input.copyTo(this->layersG.front());
	this->_fill_internal(input);
}
/**
 * \brief Fill the layers (G and DoG) from the data in the first Gaussian layer
 */
void Colloids::OctaveFinder::_fill_internal(Image &temp)
{
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

void Colloids::OctaveFinder3D::_fill_internal(Image &temp)
{
	Image temp2D(
		temp.size[0], temp.size[1]*temp.size[2],
		(PixelType*)temp.data
	);
	//iterative Gaussian blur
	for(size_t i=0; i<this->layersG.size()-1; ++i)
	{
		//Z out of place
		this->iterative_Zgaussian_filters[i].apply(temp2D, temp2D);
		//X and Y inplace
		for(int k=0; k<this->layersG[i+1].size[0]; ++k)
		{
			cv::Mat slice(
					this->layersG[i+1].size[1],
					this->layersG[i+1].size[2],
					this->layersG[i+1].type(),
					(void*)&temp2D(k)
					);
			this->iterative_gaussian_filters[i].apply(slice, slice);
		}
		//difference of Gaussians (write directly to disk)
		cv::subtract(temp, this->layersG[i], this->layers[i]);
		//write gaussian layer to disk
		temp2D.copyTo(this->layersG2D[i+1]);
	}
}

void Colloids::OctaveFinder::preblur(Image &input)
{
	this->preblur_filter->apply(input, this->layersG.front());
}

void Colloids::OctaveFinder3D::preblur(Image &input)
{
	//Z inplace
	Image input2D(
			input.size[0], input.size[1]*input.size[2],
			(PixelType*)input.data
			);
	this->preblur_Zfilter->apply(input2D, input2D);
	//X and Y inplace
	for(int k=0; k<this->layersG.front().size[0]; ++k)
	{
		cv::Mat slice(
				this->layersG.front().size[1],
				this->layersG.front().size[2],
				this->layersG.front().type(),
				(void*)&input2D(k)
				);
		this->preblur_filter->apply(slice, slice);
	}
	//write to disk
	input.copyTo(this->layersG.front());

}
/**
 * Convert an image of any kind to PixelType and process it to fill internal buffers
 */
void Colloids::OctaveFinder::preblur_and_fill(const cv::Mat &input)
{
	Image temp;
	input.convertTo(temp, temp.type());
	this->preblur_and_fill(temp);
}
/**
 * Process the input to fill internal buffers
 */
void Colloids::OctaveFinder::preblur_and_fill(Image &input)
{
	if(input.dims != this->layersG.front().dims)
	{
		std::ostringstream os;
		os << "OctaveFinder::preblur_and_fill : the input's dimension ("<<input.dims<<") must match the dimension of the finder ("<<this->layersG.front().dims<<")";
		throw std::invalid_argument(os.str().c_str());
	}
	for(int d=0; d<input.dims; ++d)
		if(input.size[d] != this->layersG.front().size[d])
		{
			std::ostringstream os;
			os << "OctaveFinder::preblur_and_fill : the input's "<<d<< "th dimension ("<<input.size[d]<<") must match the one of the finder ("<<this->layersG.front().size[d]<<")";
			throw std::invalid_argument(os.str().c_str());
		}
	this->preblur(input);
	this->_fill_internal(input);
}

/**
 * \brief Detect local minima of the scale space
 *
 * Uses the dynamic block algorythm by Neubeck and Van Gool
 * a. Neubeck and L. Van Gool, 18th International Conference On Pattern Recognition (ICPRʼ06) 850-855 (2006).
 */
void Colloids::OctaveFinder::initialize_binary(const double & max_ratio)
{
	const int nblayers = this->binary.size();
    //initialize
	this->centers_no_subpix.clear();
    for(int i = 0;i < nblayers;++i)
        this->binary[i].setTo(0);

	for(int k = 1;k < nblayers+1;k += 2)
	{
		const Image & layer0 = this->layers[k], layer1 = this->layers[k+1];
		const int si = this->sizes[k];
		for(int j = this->sizes[k]+1;j < this->get_width() - si- 1;j += 2)
		{
			boost::array<const float*, 4> ngb_ptr = {{
					&layer0(j, si+1),
					&layer0(j+1, si+1),
					&layer1(j, si+1),
					&layer1(j+1, si+1)
			}};
			for(int i = si+1;i < this->get_height() -si - 1;i += 2){
				//copy the whole neighbourhood together for locality
				boost::array<float, 8> ngb = {{
						*ngb_ptr[0]++, *ngb_ptr[0]++,
						*ngb_ptr[1]++, *ngb_ptr[1]++,
						*ngb_ptr[2]++, *ngb_ptr[2]++,
						*ngb_ptr[3]++, *ngb_ptr[3]++
				}};
				const boost::array<float, 8>::const_iterator mpos = std::min_element(ngb.begin(), ngb.end());
				if(*mpos>=0.0)
					continue;
				const int ml = mpos-ngb.begin(),
					mi = i + !!(ml&1),
					mj = j + !!(ml&2),
					mk = k + !!(ml&4);


				//maxima cannot be on the last layer or on image edges
				if(mk > nblayers || !((this->sizes[mk] <= mj) && (mj < this->get_width() - this->sizes[mk]) && (this->sizes[mk] <= mi) && (mi < this->get_height() - this->sizes[mk])))
					continue;

				bool *b = &this->binary[mk - 1](mj, mi);
				//consider only negative minima
				//with a value that is actually different from zero
				*b = (*mpos < 0) && (1 + pow(*mpos, 2) > 1);
				//remove the minima if one of its neighbours outside the block has lower value
				for(int k2 = mk - 1;k2 < mk + 2 && *b;++k2)
					for(int j2 = mj - 1;j2 < mj + 2 && *b;++j2)
						for(int i2 = mi - 1;i2 < mi + 2 && *b;++i2)
							if(k2 < k || j2 < j || i2 < i || k2 > k + 1 || j2 > j + 1 || i2 > i + 1)
								*b = *mpos <= this->layers[k2](j2, i2);




				//remove the local minima that are edges (elongated objects)
				if(*b){
					//hessian matrix
					const double hess[3] = {
							this->layers[mk](mj - 1, mi) - 2 * this->layers[mk](mj, mi) + this->layers[mk](mj + 1, mi),
							this->layers[mk](mj, mi - 1) - 2 * this->layers[mk](mj, mi) + this->layers[mk](mj, mi + 1),
							this->layers[mk](mj - 1, mi - 1) + this->layers[mk](mj + 1, mi + 1) - this->layers[mk](mj + 1, mi - 1) - this->layers[mk](mj - 1, mi + 1)};
					//determinant of the Hessian, for the coefficient see
					//H Bay, a Ess, T Tuytelaars, and L Vangool,
					//Computer Vision and Image Understanding 110, 346-359 (2008)
					const double detH = hess[0] * hess[1] - pow(hess[2], 2),
							ratio = pow(hess[0] + hess[1], 2) / (4.0 * hess[0] * hess[1]);
					*b = !((detH < 0 && 1+detH*detH > 1) || ratio > max_ratio);
					if(*b){
						std::vector<int> c(3);
						c[0] = mi;
						c[1] = mj;
						c[2] = mk;
						this->centers_no_subpix.push_back(c);
					}
				}
			}
		}
	}
}
/**
 * \brief Detect local minima of the scale space
 *
 * Uses the dynamic block algorythm by Neubeck and Van Gool
 * a. Neubeck and L. Van Gool, 18th International Conference On Pattern Recognition (ICPRʼ06) 850-855 (2006).
 */
void Colloids::OctaveFinder1D::initialize_binary(const double & max_ratio)
{
    //initialize
	this->centers_no_subpix.clear();
    for(size_t i = 0;i < this->binary.size();++i)
        this->binary[i].setTo(0);

	for(int k = 1;k < (int)this->layers.size() - 1;k += 2)
	{
		boost::array<const float*, 2> ngb_ptr = {{
							&this->layers[k](0, this->sizes[k]+1),
							&this->layers[k+1](0, this->sizes[k]+1)
		}};
		for(int i = this->sizes[k]+1;i < this->get_height() -this->sizes[k]- 1;i += 2)
		{
			//copy the whole neighbourhood together for locality
			boost::array<float, 4> ngb = {{
					*ngb_ptr[0]++, *ngb_ptr[0]++,
					*ngb_ptr[1]++, *ngb_ptr[1]++
			}};
			const boost::array<float, 4>::const_iterator mpos = std::min_element(ngb.begin(), ngb.end());
			if(*mpos>=0.0)
					continue;
			const int ml = mpos-ngb.begin(),
				mi = i + !!(ml&1),
				mk = k + !!(ml&2);
			//maxima cannot be on the last layer or on image edges
			if(mk > (int)this->binary.size() || !((this->sizes[mk] <= mi) && (mi < (int)this->get_height() - this->sizes[mk])))
				continue;
			bool *b = &this->binary[mk - 1](0, mi);
			*b = (*mpos < 0) && (1 + pow(*mpos, 2) > 1);

			//remove the minima if one of its neighbours outside the block has lower value
			for(int k2 = mk - 1;k2 < mk + 2 && *b;++k2)
				for(int i2 = mi - 1;i2 < mi + 2 && *b;++i2)
					if(k2 < k || i2 < i || k2 > k + 1 || i2 > i + 1)
						*b = this->layers[mk](0, mi) <= this->layers[k2](0, i2);
			//Eliminating edge response in 1D is more tricky
			//We look for pixels where the ratio between the Laplacian and gradient values is large
			if(*b)
			{
				const float * v = &this->layersG[mk](0, mi);
				*b = std::abs((*(v+1) + *(v-1) - 2 * *v) / (*(v+1) - *(v-1))) > 0.5;
			}
			if(*b){
				std::vector<int> c(2);
				c[0] = mi;
				c[1] = mk;
				this->centers_no_subpix.push_back(c);
			}
		}
	}
}

void Colloids::OctaveFinder3D::initialize_binary(const double & max_ratio)
{
	const int nblayers = this->binary.size();
    //initialize
	this->centers_no_subpix.clear();
    for(int i = 0;i < nblayers;++i)
        this->binary[i].setTo(0);
    //prepare the cache
    Image cache(2*2*2, this->get_height(), (PixelType)0);

	for(int l = 1;l < nblayers+1;l += 2)
	{
		const Image & layer0 = this->layers[l], layer1 = this->layers[l+1];
		const int si = this->sizes[l];
		for(int k = si+1;k < layer0.size[0] - si- 1;k += 2)
		for(int j = si+1;j < layer0.size[1] - si- 1;j += 2)
		{
			//fill the cache for the whole line of blocks
			for(int cl=0; cl<2; ++cl)
				for(int ck=0; ck<2; ++ck)
				{
					const PixelType * in = &this->layers[cl+l](ck+k, j, 0);
					PixelType * out = &cache(2*(2*cl+ck), 0);
					for(int ci=0; ci<2*cache.cols; ++ci)
						*out++ = *in++;
				}
			for(int i = si+1;i < this->get_height() -si - 1;i += 2){
				//copy the whole block together for locality
				boost::array<PixelType, 16> ngb;
				for(int u=0;u<8;++u)
				{
					ngb[2*u] = cache(u,i);
					ngb[2*u+1] = cache(u,i+1);
				}
				const boost::array<PixelType, 16>::const_iterator mpos = std::min_element(ngb.begin(), ngb.end());
				if(*mpos>=0.0)
					continue;
				const int mm = mpos-ngb.begin(),
					mi = i + !!(mm&1),
					mj = j + !!(mm&2),
					mk = k + !!(mm&4),
					ml = l + !!(mm&8);


				//maxima cannot be on the last layer or on image edges
				if(ml > nblayers || !(
						(this->sizes[ml] <= mk) && (mk < layer0.size[0] - this->sizes[ml]) &&
						(this->sizes[ml] <= mj) && (mj < layer0.size[1] - this->sizes[ml]) &&
						(this->sizes[ml] <= mi) && (mi < layer0.size[2] - this->sizes[ml])
						))
					continue;

				//bool *b = &this->binary[ml - 1](mk, mj, mi);
				//consider only negative minima
				//with a value that is actually different from zero
				bool b = (*mpos < 0) && (1 + pow(*mpos, 2) > 1);
				//remove the minima if one of its neighbours outside the block has lower value
				for(int l2 = ml - 1;l2 < ml + 2 && b;++l2)
					for(int k2 = mk - 1;k2 < mk + 2 && b;++k2)
						for(int j2 = mj - 1;j2 < mj + 2 && b;++j2)
							for(int i2 = mi - 1;i2 < mi + 2 && b;++i2)
								if(l2 < l || k2 < k || j2 < j || i2 < i || l2 > l + 1 || k2 > k +1 || j2 > j + 1 || i2 > i + 1)
									b = *mpos <= this->layers[l2](k2, j2, i2);




				//remove the local minima that are edges (elongated objects) in XY
				if(b){
					//hessian matrix
					const double hess[3] = {
							this->layers[ml](mk, mj - 1, mi) - 2 * this->layers[ml](mk, mj, mi) + this->layers[ml](mk, mj + 1, mi),
							this->layers[ml](mk, mj, mi - 1) - 2 * this->layers[ml](mk, mj, mi) + this->layers[ml](mk, mj, mi + 1),
							this->layers[ml](mk, mj - 1, mi - 1) + this->layers[ml](mk, mj + 1, mi + 1) - this->layers[ml](mk, mj + 1, mi - 1) - this->layers[ml](mk, mj - 1, mi + 1)
					};
					//determinant of the Hessian, for the coefficient see
					//H Bay, a Ess, T Tuytelaars, and L Vangool,
					//Computer Vision and Image Understanding 110, 346-359 (2008)
					const double detH = hess[0] * hess[1] - pow(hess[2], 2),
							ratio = pow(hess[0] + hess[1], 2) / (4.0 * hess[0] * hess[1]);
					b = !((detH < 0 && 1+detH*detH > 1) || ratio > max_ratio);
					//*b &= abs(
					//		(this->layers[ml](mk+1, mj, mi) - this->layers[ml](mk-1, mj, mi)) /
					//		(this->layers[ml](mk+1, mj, mi)-2*this->layers[ml](mk, mj, mi)+this->layers[ml](mk-1, mj, mi))
					//		) < 1.0;
					if(b){
						std::vector<int> c(4);
						c[0] = mi;
						c[1] = mj;
						c[2] = mk;
						c[3] = ml;
						this->centers_no_subpix.push_back(c);
					}
				}
			}
		}
	}
	//make binary image coherent
	for(std::list<std::vector<int> >::const_iterator ci = this->centers_no_subpix.begin(); ci!= this->centers_no_subpix.end(); ++ci)
		this->binary[(*ci)[3]-1]((*ci)[2], (*ci)[1], (*ci)[0]) = true;
}

void Colloids::OctaveFinder::spatial_subpix(const std::vector<int> &ci, Center_base& c) const
{
        const int i = ci[0], j = ci[1], k = ci[2];
        //When particles environment is strongly asymmetric (very close particles),
        //it is better to find the maximum of Gausian rather than the minimum of DoG.
        //If possible, we use the Gaussian layer below the detected scale
        //to have better spatial resolution
        const Image & l = (k>0 ? this->layersG[k-1] : this->layersG[k]);
        const double a[4] = {
        		(l(j, i+1) - l(j, i-1))/2.0,
        		(l(j+1, i) - l(j-1, i))/2.0,
        		l(j, i+1) -2*l(j, i) + l(j, i-1),
        		l(j+1, i) -2*l(j, i) + l(j-1, i)
        };
        c[0] = i + 0.5 - (a[2]==0 ? 0 : a[0]/a[2]);
		c[1] = j + 0.5 - (a[3]==0 ? 0 : a[1]/a[3]);
		c.intensity = this->layers[k](j, i);
}
void Colloids::OctaveFinder1D::spatial_subpix(const std::vector<int> &ci, Center_base& c) const
{
        const int i = ci[0], k = ci.back();
        //When particles environment is strongly asymmetric (very close particles),
		//it is better to find the maximum of Gausian rather than the minimum of DoG.
		//If possible, we use the Gaussian layer below the detected scale
		//to have better spatial resolution
        const Image & l = (k>0 ? this->layersG[k-1] : this->layersG[k]);
        c[0] = ci[0] + 0.5 - (l(0, i+1) - l(0, i-1)) / 2.0 / (l(0, i+1) -2*l(0, i) + l(0, i-1));
        c.intensity = this->layers[k](0, i) - 0.25 * (c[0]-i) * (l(0, i+1) - l(0, i-1));
}
void Colloids::OctaveFinder3D::spatial_subpix(const std::vector<int> &ci, Center_base& c) const
{
        const int i = ci[0], j = ci[1], k = ci[2], l = ci[3];
        const double ds = abs(c.r - l);
        //When particles environment is strongly asymmetric (very close particles),
        //it is better to find the maximum of Gausian rather than the minimum of DoG.
        //If possible, we use the Gaussian layer below the detected scale
        //to have better spatial resolution
        const Image & lay = (l>0 ? this->layersG[l-1] : this->layersG[l]);
        PixelType const * vg = &lay(k, j, i);
        for(size_t d=0; d<3; ++d)
        {
        	const size_t step = this->layers[l].step[2-d]/sizeof(PixelType);
        	const double a[3] = {*(vg-step), *vg, *(vg+step)};
        	double shift = (a[2] - a[0]) /2.0	/ (a[0] - 2 * a[1] + a[2]);
        	//in some rare cases (z edges), the shift may be very large if computed from gaussian layers
        	if(abs(shift)>0.5)
        	{
        		PixelType const * v = &this->layers[l](k, j, i);
        		double b[3] = {*(v-step), *v, *(v+step)};
        		shift = (b[2] - b[0]) /2.0	/ (b[0] - 2 * b[1] + b[2]);
        	}
			if(ds>0.25)
				shift *= (c.r<1.6)?0.985:0.99;
			else if(ds<0.1)
				shift *= 1.005;
        	c[d] = ci[d] + 0.4375 - shift;
        }
		c.intensity = this->layers[l](k, j, i);
}
double Colloids::OctaveFinder::gaussianResponse(const std::vector<int> &ci, const double & scale) const
{
        if(scale < 0)
            throw std::invalid_argument("Colloids::OctaveFinder::gaussianResponse: the scale must be positive.");
        if(ci.size()<2)
        	throw std::invalid_argument("Colloids::OctaveFinder::gaussianResponse: coordinates must be at least 2D.");

        size_t k = (size_t)(scale);
        if (k>=this->layersG.size())
        	k = this->layersG.size()-1;
		if((scale - k) * (scale - k) + 1 == 1)
			return this->layersG[k](ci[1], ci[0]);
		const double sigma = this->get_iterative_radius(scale, (double)k);
		//opencv is NOT dealing right with ROI (even if boasting about it), so we do it by hand
		const cv::Mat_<double>& kernel = get_kernel(sigma);
		const int m = kernel.rows;
		vector<double> gx(m, 0.0);
		const int xmin = max(0, (int)ci[0]+m/2+1-this->get_width()),
				xmax = min(m, (int)ci[0]+m/2+1),
				ymin = max(0, (int)ci[1]+m/2+1-this->get_height()),
				ymax = min(m, (int)ci[1]+m/2+1);
		const double *ker = &kernel(ymin,0);
		for(int y=ymin; y<ymax; ++y)
		{
			const OctaveFinder::PixelType * v = &layersG[k](ci[1]-y+m/2, ci[0]-xmin+m/2);
			for(int x=xmin; x<xmax; ++x)
				gx[x] += *v-- * *ker;
			ker++;
		}

		double resp = 0.0;
		ker = &kernel(0,0);
		for(int x=0; x<m; ++x)
			resp += gx[x] * *ker++;

        return resp;
}

double Colloids::OctaveFinder1D::gaussianResponse(const std::vector<int> &ci, const double & scale) const
{
	if(scale < 0)
		throw std::invalid_argument("Colloids::OctaveFinder::gaussianResponse: the scale must be positive.");
	if(ci.size()<1)
		throw std::invalid_argument("Colloids::OctaveFinder::gaussianResponse: coordinates must be at least 1D.");

	size_t k = (size_t)(scale);
	if (k>=this->layersG.size())
		k = this->layersG.size()-1;
	if((scale - k) * (scale - k) + 1 == 1)
		return this->layersG[k](0, ci[0]);
	const double sigma = this->get_iterative_radius(scale, (double)k);
	//opencv is NOT dealing right with ROI (even if boasting about it), so we do it by hand
	const cv::Mat_<double>& kernel = get_kernel(sigma);
	const int m = kernel.rows;
	double resp = 0.0;
	for(int x=0; x<m; ++x)
		resp += layersG[k](0, cv::borderInterpolate(ci[0]-x+m/2, this->get_height(), cv::BORDER_DEFAULT)) * kernel(x,0);
	return resp;
}

double Colloids::OctaveFinder3D::gaussianResponse(const std::vector<int> &ci, const double & scale) const
{
        if(scale < 0)
            throw std::invalid_argument("Colloids::OctaveFinder::gaussianResponse: the scale must be positive.");
        if(ci.size()<3)
        	throw std::invalid_argument("Colloids::OctaveFinder::gaussianResponse: coordinates must be at least 3D.");

        size_t k = (size_t)(scale);
        if (k>=this->layersG.size())
        	k = this->layersG.size()-1;
		if((scale - k) * (scale - k) + 1 == 1)
			return this->layersG[k](ci[2], ci[1], ci[0]);
		const double sigma = this->get_iterative_radius(scale, (double)k);
		//opencv is NOT dealing right with ROI (even if boasting about it), so we do it by hand
		const cv::Mat_<double>& kernel = get_kernel(sigma);
		const int m = kernel.rows;
		Image im(m, m, 0.0f);
		std::vector<double> gx(m, 0.0);
		const int xmin = max(0, (int)ci[0]+m/2+1-this->get_width()),
				xmax = min(m, (int)ci[0]+m/2+1),
				ymin = max(0, (int)ci[1]+m/2+1-this->get_height()),
				ymax = min(m, (int)ci[1]+m/2+1),
				zmin = max(0, (int)ci[2]+m/2+1-this->get_depth()),
				zmax = min(m, (int)ci[2]+m/2+1);

		const double *ker = &kernel(zmin,0);
		for(int z=zmin; z<zmax; ++z)
		{
			for(int y=ymin; y<ymax; ++y)
			{
				const OctaveFinder::PixelType * v = &layersG[k](ci[2]-z+m/2, ci[1]-y+m/2, ci[0]-xmin+m/2);
				OctaveFinder::PixelType	* re = &im(y, xmin);
				for(int x=xmin; x<xmax; ++x)
					*re++ += *v-- * *ker;
			}
			ker++;
		}
		ker = &kernel(ymin,0);
		for(int y=ymin; y<ymax; ++y)
		{
			const OctaveFinder::PixelType * v = &im(y, xmin);
			for(int x=xmin; x<xmax; ++x)
				gx[x] += *v++ * *ker;
			ker++;
		}

		double resp = 0.0;
		ker = &kernel(0,0);
		for(int x=0; x<m; ++x)
			resp += gx[x] * *ker++;

        return resp;
}

double Colloids::OctaveFinder::scale_subpix(const std::vector<int> &ci) const
{
    	const int &l = ci.back();
		//scale is better defined if we consider only the central pixel at different scales
		//Compute intermediate variables to do a quadratic estimate of the derivative
    	boost::array<double,8> sublayerG;
		for(size_t u = 0;u < sublayerG.size(); ++u)
			sublayerG[u] = this->gaussianResponse(ci, l - 1 + 0.5*u);
    	/*std::vector<double> scales(8);
    	for(size_t u = 0;u < scales.size(); ++u)
    		scales[u] = l - 1 + 0.5*u;
    	std::vector<double> sublayerG = this->gaussianResponse(ci, scales);*/
		boost::array<double,5> a;
		for(size_t u =0; u<a.size();++u)
			a[u] = sublayerG[u+2] - sublayerG[u];
		//Apply newton's method using quadratic estimate of the derivative
		double s = l - (-a[4] + 8*a[3] - 8*a[1] + a[0])/6.0 /(a[4]-2*a[2]+a[0]);
		//saturate the results
		if(s>l+0.5)
			s= l + 0.5;
		if(s<l-0.5)
			s= l- 0.5;
		//second round of newton's method
		if(s>=1)
		{
			if(s+0.1<l)
			{
				s= l- 0.5;
				for(size_t u = 0;u < sublayerG.size(); ++u)
					sublayerG[u] = this->gaussianResponse(ci, s - 1 + 0.5*u);
				//for(size_t u = 0;u < sublayerG.size(); ++u)
				//	scales[u] = s - 1 + 0.5*u;
				//sublayerG = this->gaussianResponse(ci, scales);
				for(size_t u =0; u<a.size();++u)
					a[u] = sublayerG[u+2] - sublayerG[u];
				s -= (-a[4] + 8*a[3] - 8*a[1] + a[0])/6.0 /(a[4]-2*a[2]+a[0]);
			}
		}
		else
		{
			//for sizes significantly below the sampled scale
			//the linear estimate of the derivative is (marginally) better
			if(s+0.25<l)
				s = l - (a[3] - a[1])/(a[4]-2*a[2]+a[0]);
		}
		if(s<l-0.5)
			s= l- 0.5;
		if(s>l+0.5)
			s= l + 0.5;
		return s;
}
double Colloids::OctaveFinder1D::scale_subpix(const std::vector<int> &ci) const
{
	//Empirical correction
	//return ci[2] + 1.1*(OctaveFinder::scale_subpix(ci)-ci[2]) - 0.1/ci[2]- 0.1;//-0.025*this->layers.size();
	const size_t l = ci.back();
	double h = 1.0/3.0;
	boost::array<double,7> a;
	for(int u=0; u<(int)a.size();++u)
		a[u] = this->gaussianResponse(ci, l - 3*h + u*h);
	double s = 2*h * (a[5] -2*a[3] + a[1])/(a[6] -3*a[4] +3*a[2] -a[0]);
	return l - 1.05*s + 0.08*pow(s,2) - pow(2,-2/(double)this->get_n_layers())+0.025*l -0.025;
}
double Colloids::OctaveFinder3D::scale_subpix(const std::vector<int> &ci) const
{
    	const int &l = ci.back();
		//scale is better defined if we consider only the central pixel at different scales
    	double shift = 0.0;
		//use only a linear estimate of the derivative
		boost::array<double,3> a;
		for(int u = l-1;u < l+2; ++u)
			a[u-l+1] = this->layers[u](ci[2], ci[1], ci[0]);
		shift = - (a[2] - a[0]) /2.0 /(a[2] -2*a[1] +a[0]);
		/*if(shift>-0.4 && shift<0.3)
		{
			shift -= 0.01*l;
		}*/
		if(shift>-0.4)
		{
			if(shift>-0.3 && shift<0.2)
				shift -= 0.015*l;
			else if(shift<0.3)
				shift -= 0.005*l;
		}
		//avoid overflow
		if(shift < -0.5)
			shift = - 0.5;
		if(shift > 0.5)
			shift = 0.5;
		return l+shift;
}

void Colloids::OctaveFinder::single_subpix(const std::vector<int> &ci, Center_base &c) const
{
	c.r = this->scale_subpix(ci);
	this->spatial_subpix(ci, c);
}

const double OctaveFinder::get_iterative_radius(const double & larger, const double & smaller) const
{
	return this->preblur_radius * sqrt(pow(2.0, 2.0*larger/this->get_n_layers()) - pow(2.0, 2.0*smaller/this->get_n_layers()));
}

const std::vector<int> OctaveFinder::get_center_pixel(const size_t n) const
{
	std::list<std::vector<int> >::const_iterator it = this->centers_no_subpix.begin();
	std::advance(it, n);
	return *(it);
}
/**
 * \brief Eliminate pixel centers duplicated at the seam between two Octaves
 */
void OctaveFinder::seam_binary(OctaveFinder & other)
{
	const double sizefactor = this->get_height()/other.get_height();
	OctaveFinder & highRes = sizefactor>1?(*this):other, &lowRes = sizefactor>1?other:(*this);
	const double sf = sizefactor>1?sizefactor:(1.0/sizefactor);
	//centers in highRes that corresponds to a lower intensity in lowRes
	std::list<std::vector<int> >::iterator c = highRes.centers_no_subpix.begin();
	while(c!=highRes.centers_no_subpix.end())
	{
		if(
			(c->back() == (int)highRes.get_n_layers()) &&
			(lowRes.binary.front()((*c)[1]/sf+0.5, (*c)[0]/sf+0.5)) &&
			(highRes.layers[c->back()]((*c)[1], (*c)[0]) > lowRes.layers[1]((*c)[1]/sf+0.5, (*c)[0]/sf+0.5))
			)
		{
			highRes.binary[c->back()-1]((*c)[1], (*c)[0]) = false;
			c = highRes.centers_no_subpix.erase(c);
		}
		else c++;
	}
	//centers in lowRes that corresponds to a lower intensity in highRes
	c = lowRes.centers_no_subpix.begin();
	while(c!=lowRes.centers_no_subpix.end())
	{
		if(c->back() == 1)
		{
			//we must look at all pixels of highRes that overlap with the pixel of lowRes
			bool *b = &lowRes.binary.front()((*c)[1], (*c)[0]);
			const PixelType vb = lowRes.layers[1]((*c)[1], (*c)[0]);
			for(size_t j=max(0, (int)(((*c)[1]-1)*sf)); j<(size_t)(((*c)[1]+1)*sf) && j<(size_t)highRes.get_width(); ++j)
				for(size_t i=max(0, (int)(((*c)[0]-1)*sf)); i<(size_t)(((*c)[0]+1)*sf) && i<(size_t)highRes.get_height(); ++i)
					*b &= !(highRes.binary.back()(j, i) && (vb > highRes.layers[highRes.get_n_layers()](j, i)));

			if(*b)
				c++;
			else
				c = lowRes.centers_no_subpix.erase(c);
		}
		else c++;
	}

}

void OctaveFinder1D::seam_binary(OctaveFinder & other)
{
	const double sizefactor = this->get_height()/other.get_height();
	OctaveFinder1D & highRes = sizefactor>1?(*this):dynamic_cast<OctaveFinder1D&>(other),
			&lowRes = sizefactor>1?dynamic_cast<OctaveFinder1D&>(other):(*this);
	const double sf = sizefactor>1?sizefactor:(1.0/sizefactor);
	//centers in highRes that corresponds to a lower intensity in lowRes
	std::list<std::vector<int> >::iterator c = highRes.centers_no_subpix.begin();
	while(c!=highRes.centers_no_subpix.end())
	{
		if(
			(c->back() == (int)highRes.get_n_layers()) &&
			(lowRes.binary.front()(0, (*c)[0]/sf+0.5)) &&
			(highRes.layers[c->back()](0, (*c)[0]) > lowRes.layers[1](0, (*c)[0]/sf+0.5))
			)
		{
			highRes.binary[c->back()-1](0, (*c)[0]) = false;
			c = highRes.centers_no_subpix.erase(c);
		}
		else c++;
	}
	//centers in lowRes that corresponds to a lower intensity in highRes
	c = lowRes.centers_no_subpix.begin();
	while(c!=lowRes.centers_no_subpix.end())
	{
		if(c->back() == 1)
		{
			//we must look at all pixels of highRes that overlap with the pixel of lowRes
			bool *b = &lowRes.binary.front()(0, (*c)[0]);
			const PixelType vb = lowRes.layers[1](0, (*c)[0]);
			for(size_t i=max(0, (int)(((*c)[0]-1)*sf)); i<(size_t)(((*c)[0]+1)*sf) && i<(size_t)highRes.get_height(); ++i)
				*b &= !(highRes.binary.back()(0, i) && (vb > highRes.layers[highRes.get_n_layers()](0, i)));

			if(*b)
				c++;
			else
				c = lowRes.centers_no_subpix.erase(c);
		}
		else c++;
	}

}
void OctaveFinder3D::seam_binary(OctaveFinder & other)
{
	const double sizefactor = this->get_height()/other.get_height();
	OctaveFinder3D & highRes = sizefactor>1?(*this):dynamic_cast<OctaveFinder3D&>(other),
				&lowRes = sizefactor>1?dynamic_cast<OctaveFinder3D&>(other):(*this);
	const double sf = sizefactor>1?sizefactor:(1.0/sizefactor);
	//centers in highRes that corresponds to a lower intensity in lowRes
	std::list<std::vector<int> >::iterator c = highRes.centers_no_subpix.begin();
	while(c!=highRes.centers_no_subpix.end())
	{
		if(
			(c->back() == (int)highRes.get_n_layers()) &&
			(lowRes.binary.front()((*c)[2]/sf+0.5, (*c)[1]/sf+0.5, (*c)[0]/sf+0.5)) &&
			(highRes.layers[c->back()]((*c)[2], (*c)[1], (*c)[0]) > lowRes.layers[1]((*c)[2]/sf+0.5, (*c)[1]/sf+0.5, (*c)[0]/sf+0.5))
			)
		{
			highRes.binary[c->back()-1]((*c)[2], (*c)[1], (*c)[0]) = false;
			c = highRes.centers_no_subpix.erase(c);
		}
		else c++;
	}
	//centers in lowRes that corresponds to a lower intensity in highRes
	c = lowRes.centers_no_subpix.begin();
	while(c!=lowRes.centers_no_subpix.end())
	{
		if(c->back() == 1)
		{
			//we must look at all pixels of highRes that overlap with the pixel of lowRes
			bool *b = &lowRes.binary.front()((*c)[2], (*c)[1], (*c)[0]);
			const PixelType vb = lowRes.layers[1]((*c)[2], (*c)[1], (*c)[0]);
			for(size_t k=max(0, (int)(((*c)[2]-1)*sf)); k<(size_t)(((*c)[2]+1)*sf) && k<(size_t)highRes.get_depth(); ++k)
				for(size_t j=max(0, (int)(((*c)[1]-1)*sf)); j<(size_t)(((*c)[1]+1)*sf) && j<(size_t)highRes.get_width(); ++j)
					for(size_t i=max(0, (int)(((*c)[0]-1)*sf)); i<(size_t)(((*c)[0]+1)*sf) && i<(size_t)highRes.get_height(); ++i)
						*b &= !(highRes.binary.back()(k, j, i) && (vb > highRes.layers[highRes.get_n_layers()](k, j, i)));

			if(*b)
				c++;
			else
				c = lowRes.centers_no_subpix.erase(c);
		}
		else c++;
	}

}

void Colloids::OctaveFinder::fill_iterative_radii()
{
		//target blurring radii
        vector<double> sigmas(this->get_n_layers()+3);
        for(size_t i=0; i<sigmas.size(); ++i)
        	sigmas[i] = this->preblur_radius * pow(2, i/double(this->get_n_layers()));
        //corresponding blob sizes
        const double n = this->get_n_layers();
        if(this->get_height()==1 || this->get_width()==1)
        	this->prefactor = sqrt(2.0 * log(2.0) / n / (pow(2.0, 2.0 / n) - 1));
        else
        	this->prefactor = 2.0 * sqrt(log(2.0) / n / (pow(2.0, 2.0 / n) - 1));
        vector<int>::iterator si = this->sizes.begin();
        for(vector<double>::const_iterator sig=sigmas.begin(); sig!=sigmas.end(); sig++)
        	*si++ = static_cast<int>((*sig * prefactor) + 0.5);
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

void Colloids::OctaveFinder3D::fill_iterative_radii()
{
	const double n = this->get_n_layers();
	//target blurring radii
	vector<double> sigmas(this->get_n_layers()+3);
	for(size_t i=0; i<sigmas.size(); ++i)
		sigmas[i] = this->preblur_radius * pow(2, i/n);
	//corresponding blob sizes
	this->prefactor = sqrt(6.0 * log(2.0) / n / (pow(2.0, 2.0 / n) - 1));
	vector<int>::iterator si = this->sizes.begin();
	for(vector<double>::const_iterator sig=sigmas.begin(); sig!=sigmas.end(); sig++)
		*si++ = static_cast<int>((*sig * prefactor) + 0.5);
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
	cv::Mat_<double> kx(1,1, 1.0);
	for(size_t i=0; i<this->iterative_radii.size(); ++i)
		this->iterative_Zgaussian_filters[i] = *createSeparableLinearFilter
		(
				this->layersG[0].type(), this->layersG[0].type(),
				kx, get_kernel(this->iterative_radii[i]),
				cv::Point(-1,-1), 0, cv::BORDER_DEFAULT
		);
	this->preblur_Zfilter = createSeparableLinearFilter
	(
			this->layersG[0].type(), this->layersG[0].type(),
			kx, get_kernel(this->preblur_radius),
			cv::Point(-1,-1), 0, cv::BORDER_DEFAULT
	);
}



}; //namespace





