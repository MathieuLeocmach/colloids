#include "octavefinder.hpp"
#include <boost/array.hpp>
#include <algorithm>
#include <numeric>
#include <stdexcept>
# ifdef _OPENMP
#include <omp.h>
#endif


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
		this->binary.push_back(cv::Mat_<bool>(nrows, ncols, (bool)0));
	this->set_radius_preblur(preblur_radius);
}

OctaveFinder3D::OctaveFinder3D(const int nplanes, const int nrows, const int ncols, const int nbLayers, const double &preblur_radius, bool incore) :
	OctaveFinder(0, 0, nbLayers, preblur_radius), iterative_Zgaussian_filters(nbLayers+2), ZXratio(1.0), halfZpreblur(false), deconv(false)
{
	const size_t nbpixels =  nplanes * nrows * ncols;
	if(incore)
	{
		data = new char[nbpixels * sizeof(PixelType) * (nbLayers + 3)];
	}
	else
	{
		//random file name in the working directory
		do
		{
			this->path.clear();
			this->path.reserve(30);
			this->path.push_back('_');
			this->path.push_back('_');
			for(int i=0; i<28;++i)
				this->path.push_back('a'+rand()%('Z'-'a'));
		} while(std::ifstream(path.c_str()).good());
		//create a memory mapped file to contain the images data
		boost::iostreams::mapped_file_params params(this->path);
		params.new_file_size = nbpixels * sizeof(PixelType) * (nbLayers + 3);
		params.flags = boost::iostreams::mapped_file::readwrite;
		this->file.open(params);
		this->data = this->file.data();
	}
	//create the images inside the memory mapped file
	int dims[3] = {nplanes, nrows, ncols};
	for (int i = 0; i<nbLayers+3; ++i)
		this->layersG[i] = cv::Mat(
				3, dims, layersG[i].type(),
				(void*)(this->data + i * nbpixels * sizeof(PixelType))
				);
	for (int i = 0; i<nbLayers+2; ++i)
		this->layers[i] = cv::Mat(
				3, dims, layers[i].type(),
				(void*)(this->data + (i + layersG.size()) * nbpixels * sizeof(PixelType))
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
	if(!this->path.empty())
		remove(this->path.c_str());
	else
		delete[] this->data;
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

/**
 * \brief Gaussian blur of each XY plane of a 3D image
 */
void inplace_blurXY(cv::Mat &im, const double &radius)
{
	cv::Mat temp2D(
		im.size[0], im.size[1]*im.size[2],
		im.type(), (void*)im.data
	);
	# ifdef _OPENMP
	#pragma omp parallel
	{
		cv::Ptr<cv::FilterEngine> filter = createSeparableLinearFilter
		(
			im.type(), im.type(),
			OctaveFinder::get_kernel(radius), OctaveFinder::get_kernel(radius),
			cv::Point(-1,-1), 0, cv::BORDER_DEFAULT
		);
		#pragma omp for
		for(int k=0; k<im.size[0]; ++k)
		{
			cv::Mat slice(im.size[1], im.size[2], im.type(), (void*)temp2D.ptr(k));
			filter->apply(slice, slice);
		}
	}
	#else
	cv::Ptr<cv::FilterEngine> filterXY = createSeparableLinearFilter
	(
		im.type(), im.type(),
		OctaveFinder::get_kernel(radius), OctaveFinder::get_kernel(radius),
		cv::Point(-1,-1), 0, cv::BORDER_DEFAULT
	);
	for(int k=0; k<im.size[0]; ++k)
	{
		cv::Mat slice(im.size[1], im.size[2], im.type(), (void*)temp2D.ptr(k));
		filterXY->apply(slice, slice);
	}
	#endif
}

void inplace_blur3D(cv::Mat &im, const double &radius, const double &ZXratio)
{
	cv::Mat temp2D(
		im.size[0], im.size[1]*im.size[2],
		im.type(), (void*)im.data
	);
	//Z
	# ifdef _OPENMP
	#pragma omp parallel
	{
		cv::Mat_<double> kx(1,1, 1.0);
		cv::Ptr<cv::FilterEngine> filter = createSeparableLinearFilter
		(
			im.type(), im.type(),
			kx, OctaveFinder::get_kernel(radius/ZXratio),
			cv::Point(-1,-1), 0, cv::BORDER_DEFAULT
		);
		int sectionsize = temp2D.cols/omp_get_num_threads();
		cv::Rect roi = cv::Rect(
				omp_get_thread_num()*sectionsize, 0,
				sectionsize, temp2D.rows
				)&cv::Rect(0,0,temp2D.cols, temp2D.rows);
		cv::Mat dst(temp2D, roi);
		filter->apply(dst, dst);
	}
	#else
	//Z
	cv::Mat_<double> kx(1,1, 1.0);
	cv::Ptr<cv::FilterEngine> filterZ = createSeparableLinearFilter
	(
		im.type(), im.type(),
		kx, OctaveFinder::get_kernel(radius/ZXratio),
		cv::Point(-1,-1), 0, cv::BORDER_DEFAULT
	);
	filterZ->apply(temp2D, temp2D);
	#endif
	//XY
	inplace_blurXY(im, radius);
}

void Colloids::OctaveFinder3D::set_deconv(bool value)
{
	if(value && ((int)this->deconvKernel.size() != this->get_height()/2+1))
		throw std::invalid_argument("Load the kernel before");
	this->deconv=value;
}

void Colloids::OctaveFinder3D::load_deconv_kernel(const std::vector<PixelType> &kernel)
{
	if((int)kernel.size() != this->get_height()/2+1)
		throw std::invalid_argument("The kernel size must match the last dimension of the finder");
	this->deconvKernel = kernel;
}

void Colloids::OctaveFinder3D::_fill_internal(Image &temp)
{
	// TODO	Add here the z-deconvolution code from known kernel

	//iterative Gaussian blur
	for(size_t i=0; i<this->layersG.size()-1; ++i)
	{
		inplace_blur3D(temp, this->iterative_radii[i], this->ZXratio);
		//write gaussian layer to disk
		temp.copyTo(this->layersG[i+1]);
	}
}

void Colloids::OctaveFinder::preblur(Image &input)
{
	this->preblur_filter->apply(input, this->layersG.front());
}

void Colloids::OctaveFinder3D::preblur(Image &input)
{
	if(this->halfZpreblur)
		inplace_blur3D(input, this->preblur_radius, 0.5);
	else
		inplace_blur3D(input, this->preblur_radius, this->ZXratio);
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
			boost::array<const PixelType*, 4> ngb_ptr = {{
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
		boost::array<const PixelType*, 2> ngb_ptr = {{
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
				const PixelType * v = &this->layersG[mk](0, mi);
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
	const int nblayers = this->layersG.size()-3;
    //initialize
	this->centers_no_subpix.clear();
	this->centers.clear();

	//In 3D, the Gaussian layers are stored in (memory mapped) files,
	//so we have to access the data in large chunks to avoid disk latency.
	//We cannot load all layers in memory, but loading all the scales is better to compute DoG.
	//We load 8 consecutive planes at all scales:
	//	3 planes to compute the third order estimate of z derivative of Gaussian (subpixel resolution in z)
	//	2 planes to have dynamic blocks of depth 2 in DoG
	//	3 planes to compute the third order estimate of z derivative of Gaussian
	CircularZ4D circ(this->layersG.size(), this->layersG.front().size[1], this->layersG.front().size[2]);
	//initial fill of the 6 first planes at every scale
	for(int l=0; l<(int)this->layersG.size(); ++l)
		circ.loadplanes(&this->layersG[l](this->sizes[1]-3, 0, 0), l, -3, 6);

	//dynamic block algorithm in 4D
	for(int k=this->sizes[1]; k<this->layersG.front().size[0] - this->sizes[1]-1; k += 2)
	{
		//load the 2 next planes at every scale
		for(int l=0; l<(int)this->layersG.size(); ++l)
			circ.loadplanes(&this->layersG[l](k+3, 0, 0), l);

		//look for local minima in DoG
		for(int l=1; l<(int)this->layersG.size()-2; l+=2)
			for(int j = this->sizes[l]; j < this->layersG.front().size[1] - this->sizes[l]- 1; j += 2)
				for(int i = this->sizes[l]; i < this->layersG.front().size[2] - this->sizes[l]- 1; i += 2)
				{
					//DoG block
					int ml, mk, mj, mi;
					OctaveFinder::PixelType value;
					circ.blockmin(l, j, i, ml, mk, mj, mi, value);
					//consider only negative minima with a value that is actually different from zero
					if(value>0 || (1+value*value==1) ||
							//minima cannot be on the last layer or on image edges
							(ml > nblayers) || !(
							(this->sizes[ml] <= mk+k) && (mk+k < this->layersG.front().size[0] - this->sizes[ml]) &&
							(this->sizes[ml] <= mj) && (mj < this->layersG.front().size[1] - this->sizes[ml]) &&
							(this->sizes[ml] <= mi) && (mi < this->layersG.front().size[2] - this->sizes[ml])
							))
						continue;
					//compare with the DoG pixels outside the block
					if(!circ.is_localmin(l, j, i, ml, mk, mj, mi, value))
						continue;
					//remove the local minima that are edges (elongated objects) in XY
					if(circ.is_edge(ml, mk, mj, mi, max_ratio))
						continue;

					//we finally have a valid local minimum that we register as a center
					std::vector<int> ci(4);
					ci[0] = mi;
					ci[1] = mj;
					ci[2] = mk+k;
					ci[3] = ml;
					this->centers_no_subpix.push_back(ci);
					//subpixel resolution
					Center3D c;
					c.intensity = value;
					c[0] = ci[0] + circ.shift(ml, mk, mj, mi, 0);
					c[1] = ci[1] + circ.shift(ml, mk, mj, mi, 1);
					c[2] = ci[2] + circ.shift(ml, mk, mj, mi, 2);
					c.r = ml + circ.shift(ml, mk, mj, mi, 3);
					this->centers.push_back(c);
				} //end of finding local minima

		//prepare next step
		++circ;
	}
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
        //When particles environment is strongly asymmetric (very close particles),
        //it is better to find the maximum of Gausian rather than the minimum of DoG.
        //If possible, we use the Gaussian layer below the detected scale
        //to have better spatial resolution
        const Image & lay = this->layersG[l-1];
        PixelType const * vg = &lay(k, j, i);
        for(size_t d=0; d<3; ++d)
        {
        	const size_t step = this->layersG[l].step[2-d]/sizeof(PixelType);
        	//const double a[3] = {*(vg-step), *vg, *(vg+step)};
        	//double shift = (a[2] - a[0]) /2.0	/ (a[0] - 2 * a[1] + a[2]);
        	//const double a[5] = {*(vg-2*step), *(vg-step), *vg, *(vg+step), *(vg+2*step)};
        	//double shift = (-a[4] + 8*a[3] - 8*a[1] +a[0]) /12.0 / (a[3] - 2 * a[2] + a[1]);
        	const double a[7] = {*(vg-3*step), *(vg-2*step), *(vg-step), *vg, *(vg+step), *(vg+2*step), *(vg+3*step)};
        	double shift = (a[6] - 9*a[5] + 45*a[4] - 45*a[2] + 9*a[1] -a[0]) /60.0 / (a[6]/90 -3*a[5]/20 + 1.5*a[4] - 49*a[3]/18 + 1.5*a[2] -3*a[1]/20 + a[0]/90);
        	//prevent overflow
        	if(shift>1)
        		shift=1;
        	if(shift<-1)
        		shift=-1;
        	c[d] = ci[d] + 0.4375 - shift;
        }
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
		boost::array<double,3> b;
		for(int u = l-1;u < l+2; ++u)
			b[u-l+1] = this->layersG[u+1](ci[2], ci[1], ci[0])-this->layersG[u](ci[2], ci[1], ci[0]);
		shift = - (b[2] - b[0]) /2.0 /(b[2] -2*b[1] +b[0]);
		//empirical correction
		shift -= 0.045*l/(double)this->get_n_layers();
		shift = 1.07*shift + 0.235*shift*shift -0.29*pow(shift,3) -0.30*pow(shift,4);
		/*if(shift>-0.4 && shift<0.3)
		{
			shift -= 0.01*l;
		}*/
		/*if(shift>-0.4)
		{
			if(shift>-0.3 && shift<0.2)
				shift -= 0.015*l;
			else if(shift<0.3)
				shift -= 0.005*l;
		}*/
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
				kx, get_kernel(this->iterative_radii[i]/this->ZXratio),
				cv::Point(-1,-1), 0, cv::BORDER_DEFAULT
		);
	this->preblur_Zfilter = createSeparableLinearFilter
	(
			this->layersG[0].type(), this->layersG[0].type(),
			kx, get_kernel(this->preblur_radius/this->ZXratio),
			cv::Point(-1,-1), 0, cv::BORDER_DEFAULT
	);
}

CircularZ4D::CircularZ4D(int nLayers, int nrows, int ncols)
{
	z0=3;
	int dims[4] = {nLayers, 8, nrows, ncols};
	this->gaussians = Image(4, dims, PixelType(0));
}

void CircularZ4D::loadplanes(const CircularZ4D::PixelType* input, const int &l, const int & k, const int & nplanes)
{
	if((k+z0+8)%8 +nplanes > 8)
		throw std::invalid_argument("That much planes crosses the periodic boundary conditions");
	std::copy(
			input,
			input+nplanes*this->gaussians.size[3]*this->gaussians.size[2],
			(PixelType*)(
						this->gaussians.data
						+ l*this->gaussians.step[0]
						+ ((k+z0+8)%8) * this->gaussians.step[1])
			);
}
void CircularZ4D::blockmin(const int &l, const int &j, const int &i, int &ml, int &mk, int &mj, int &mi, CircularZ4D::PixelType& value) const
{
	//DoG block starting in (l,0, j, i)
	boost::array<PixelType, 16> block;
	boost::array<PixelType, 16>::iterator bl_it = block.begin();
	for(int l2=0; l2<2; ++l2)
		for(int k2=0; k2<2; ++k2)
			for(int j2=0; j2<2; ++j2)
				for(int i2=0; i2<2; ++i2)
					*bl_it++ = this->getDoG(l+l2, k2, j+j2, i+i2);
	//position of the minimum inside the block
	const boost::array<PixelType, 16>::const_iterator mpos = std::min_element(block.begin(), block.end());
	const int mm = mpos-block.begin();
	mi = i + !!(mm&1),
	mj = j + !!(mm&2),
	mk = !!(mm&4),
	ml = l + !!(mm&8);
	value = *mpos;
}
bool CircularZ4D::is_localmin(
		const int &l, const int &j, const int &i,
		const int &ml, const int &mk, const int &mj, const int &mi,
		const CircularZ4D::PixelType &value) const
{
	//remove the minima if one of its neighbours outside the block has lower value
	for(int l2 = ml - 1;l2 < ml + 2; ++l2)
		for(int k2 = mk - 1;k2 < mk + 2; ++k2)
			for(int j2 = mj - 1;j2 < mj + 2; ++j2)
				for(int i2 = mi - 1;i2 < mi + 2; ++i2)
					if(l2 < l || k2 < 0 || j2 < j || i2 < i || l2 > l + 1 || k2 > 1 || j2 > j + 1 || i2 > i + 1)
						if(value > this->getDoG(l2, k2, j2, i2))
							return false;
	return true;
}
/**
 * Test if a given the local minimum is an edge (elongated objects) in XY
 */
bool CircularZ4D::is_edge(const int &ml, const int &mk, const int &mj, const int &mi, const double & max_ratio) const
{
	//hessian matrix in XY
	const double hess[3] = {
			this->getDoG(ml, mk, mj-1, mi) - 2* this->getDoG(ml, mk, mj, mi) + this->getDoG(ml, mk, mj+1, mi),
			this->getDoG(ml, mk, mj, mi-1) - 2* this->getDoG(ml, mk, mj, mi) + this->getDoG(ml, mk, mj, mi+1),
			this->getDoG(ml, mk, mj-1, mi-1) + this->getDoG(ml, mk, mj+1, mi+1)
			- this->getDoG(ml, mk, mj+1, mi-1) - this->getDoG(ml, mk, mj-1, mi+1)
	};
	//determinant of the Hessian, for the coefficient see
	//H Bay, a Ess, T Tuytelaars, and L Vangool,
	//Computer Vision and Image Understanding 110, 346-359 (2008)
	const double detH = hess[0] * hess[1] - pow(hess[2], 2);
	if(detH<0 && 1+detH*detH>1)
		return true;
	const double ratio = pow(hess[0] + hess[1], 2) / (4.0 * hess[0] * hess[1]);
	if(ratio > max_ratio)
		return true;
	return false;
}
double CircularZ4D::shift(const int &ml, const int &mk, const int &mj, const int &mi, const int&d) const
{
	double shift = 0.0;
	if(d>2)
	{
		//subscale resolution uses DoG
		boost::array<double,3> a;
		for(int u = 0;u < 3; ++u)
			a[u] = this->getDoG(ml-1+u, mk, mj, mi);
		//use only a linear estimate of the derivative
		shift = - (a[2] - a[0]) /2.0 /(a[2] -2*a[1] +a[0]);
		//empirical correction
		shift -= 0.045*ml/(double)(this->gaussians.size[0]-2);
		shift = 1.07*shift + 0.235*shift*shift -0.29*pow(shift,3) -0.30*pow(shift,4);
		//avoid overflow
		if(shift < -0.5)
			shift = - 0.5;
		if(shift > 0.5)
			shift = 0.5;
	}
	else
	{
		//use the Gaussian layer below the detected scale
        //to have better spatial resolution when particles are close
		boost::array<double, 7> a;
		switch(d)
		{
		case 0:
			for(size_t u=0; u<a.size(); ++u)
				a[u] = this->getG(ml-1, mk, mj, mi-3+u);
			break;
		case 1:
			for(size_t u=0; u<a.size(); ++u)
				a[u] = this->getG(ml-1, mk, mj-3+u, mi);
			break;
		case 2:
			for(size_t u=0; u<a.size(); ++u)
				a[u] = this->getG(ml-1, mk-3+u, mj, mi);
			break;
		}
		//third order estimate of the first and second derivatives
		shift = (a[6] - 9*a[5] + 45*a[4] - 45*a[2] + 9*a[1] -a[0]) /60.0 / (a[6]/90 -3*a[5]/20 + 1.5*a[4] - 49*a[3]/18 + 1.5*a[2] -3*a[1]/20 + a[0]/90);
		//prevent overflow
		if(shift>1)
			shift=1;
		if(shift<-1)
			shift=-1;
		shift = 0.4375 - shift;
	}
	return shift;
}

}; //namespace





