#ifndef OCTAVEFINDER_H
#define OCTAVEFINDER_H

#include "center.h"
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <list>
#include <numeric>

namespace Colloids
{

    class OctaveFinder : boost::noncopyable
    {
        public:
			typedef float PixelType;
			typedef cv::Mat_<PixelType> Image;

            OctaveFinder(const int nrows=256, const int ncols=256, const int nbLayers=3, const double &preblur_radius=1.6);
            virtual ~OctaveFinder();

            //accessors
            inline const int & get_width() const {return this->layers[0].size[this->layers[0].dims-2];};
            inline const int & get_height() const {return this->layers[0].size[this->layers[0].dims-1]; };
            inline const size_t get_n_layers() const {return this->layers.size()-2;};
            const double & get_radius_preblur() const {return this->preblur_radius;}
            virtual void set_radius_preblur(const double &k=1.6);
            const double & get_prefactor() const {return this->prefactor;}
            const double & get_iterative_radius(const size_t l) const {return this->iterative_radii[l];};
            const double get_iterative_radius(const double &larger, const double &smaller) const;
            const size_t & get_size(const size_t l) const {return this->sizes[l];};
            inline const cv::Mat_<bool> & get_binary(const size_t l) const {return binary[l-1];};
			inline const Image & get_layersG(const size_t l) const {return layersG[l];}
			inline const Image & get_layers(const size_t l) const {return layers[l];}
			inline const size_t get_nb_centers() const {return centers_no_subpix.size();}
			const std::vector<int> get_center_pixel(const size_t n) const;
			static const cv::Mat_<double>& get_kernel(const double &sigma);

            //processing
            void fill(const cv::Mat &input);
            void preblur_and_fill(const cv::Mat &input);
            virtual void initialize_binary(const double &max_ratio = 1.2);
            template<int D>
            inline void subpix(std::vector<Center<D> > &centers) const;
            void single_subpix(const std::vector<int> &ci, Center_base &c) const;
            virtual void spatial_subpix(const std::vector<int> &ci, Center_base& c) const;
            virtual double scale_subpix(const std::vector<int> &ci) const;
            void scale(Center_base &c) const{
            	c.r = this->prefactor * this->preblur_radius * pow(2.0, (c.r + 1) / this->get_n_layers());
            };
            template<int D>
            std::vector<Center<D> > get_centers(const cv::Mat &input, const bool preblur=false);
            virtual double gaussianResponse(const std::vector<int> &ci, const double & scale) const;
            virtual std::vector<double> gaussianResponse(const std::vector<int> &ci, const std::vector<double> & scales) const;
            virtual void seam_binary(OctaveFinder & other);


    protected:
            std::vector<Image > layersG, layers;
            std::vector<cv::Mat_<bool> > binary;
            std::vector<double> iterative_radii;
            std::vector<cv::FilterEngine> iterative_gaussian_filters;
            std::vector<size_t> sizes;
            std::list<std::vector<int> > centers_no_subpix;
            double preblur_radius, prefactor;
            static std::map<size_t, cv::Mat_<double> > kernels;
            cv::Ptr<cv::FilterEngine> preblur_filter;

            virtual void _fill_internal(Image &temp);
            virtual void preblur(Image &input);
            virtual void fill_iterative_radii();
    };

    class OctaveFinder1D : public OctaveFinder
    {
		public:
			OctaveFinder1D(const int ncols=256, const int nbLayers=3, const double &preblur_radius=1.6) :
				OctaveFinder(1, ncols, nbLayers, preblur_radius){};

			virtual void initialize_binary(const double &max_ratio = 1.1);
			virtual void spatial_subpix(const std::vector<int> &ci, Center_base& c) const;
			virtual double scale_subpix(const std::vector<int> &ci) const;
			virtual double gaussianResponse(const std::vector<int> &ci, const double & scale) const;
			virtual std::vector<double> gaussianResponse(const std::vector<int> &ci, const std::vector<double> & scales) const;
			virtual void seam_binary(OctaveFinder & other);
    };

    class OctaveFinder3D : public OctaveFinder
    {
		public:
			OctaveFinder3D(const int nplanes=256, const int nrows=256, const int ncols=256, const int nbLayers=3, const double &preblur_radius=1.6);
			virtual ~OctaveFinder3D();

			inline const int & get_depth() const {return this->layers[0].size[this->layers[0].dims-3];};

			virtual void initialize_binary(const double &max_ratio = 1.1);
			virtual void spatial_subpix(const std::vector<int> &ci, Center_base& c) const;
			//virtual double scale_subpix(const std::vector<int> &ci) const;
			virtual double gaussianResponse(const std::vector<int> &ci, const double & scale) const;
			virtual std::vector<double> gaussianResponse(const std::vector<int> &ci, const std::vector<double> & scales) const;
			virtual void seam_binary(OctaveFinder & other);

		protected:
			boost::iostreams::mapped_file file;
			std::string path;
			std::vector<Image > layersG2D;
			std::vector<cv::FilterEngine> iterative_Zgaussian_filters;
			cv::Ptr<cv::FilterEngine> preblur_Zfilter;

			virtual void fill_iterative_radii();
			virtual void preblur(Image &input);
			virtual void _fill_internal(Image &temp);
    };

    template<int D>
    inline void OctaveFinder::subpix(std::vector<Center<D> > &centers) const
	{
		centers.clear();
		centers.resize(this->centers_no_subpix.size());
		//subpixel resolution in pixel units
		std::list<std::vector<int> >::const_iterator ci = this->centers_no_subpix.begin();
		for(size_t c=0; c<centers.size(); ++c)
			this->single_subpix(*ci++, centers[c]);
	}
    template<int D>
    std::vector<Center<D> > OctaveFinder::get_centers(const cv::Mat & input, const bool preblur)
    {
    	if(preblur)
    		this->preblur_and_fill(input);

    	else
    		this->fill(input);

    	this->initialize_binary();
    	std::vector<Center<D> > centers;
    	this->subpix(centers);
    	for(size_t c=0; c<centers.size(); ++c)
    		this->scale(centers[c]);
    	return centers;
    }
    template<>
	inline void OctaveFinder::subpix(std::vector<Center3D > &centers) const
	{
    	typedef std::list<std::vector<int> >::const_iterator list_it;
		centers.clear();
		centers.resize(this->centers_no_subpix.size());
		//spatial subpix
		list_it ci = this->centers_no_subpix.begin();
		for(size_t c=0; c<centers.size(); ++c)
			this->spatial_subpix(*ci++, centers[c]);

		//subscale responses for each particle
		Image sublayerG(this->centers_no_subpix.size(), this->layersG.size()*2-1, (PixelType)0);
		//easy coefficients : the ones from existing Gaussian layers
		ci = this->centers_no_subpix.begin();
		for(size_t c=0; c<this->centers_no_subpix.size(); ++c)
		{
			for(size_t l=0; l<this->layersG.size(); ++l)
				sublayerG(c, 2*l) = this->layersG[l]((*ci)[2], (*ci)[1], (*ci)[0]);
			ci++;
		}
		//additional coefficients : intermediate Gaussian layers must be partially computed
		int dims[3] = {this->layersG.front().size[0], this->layersG.front().size[1], this->layersG.front().size[2]};
		Image layG2D = Image::zeros(dims[0], dims[1]*dims[2]);
		cv::Mat_<double> kx(1,1, 1.0);
		for(size_t l=0; l+1<this->layersG.size(); ++l)
		{
			const cv::Mat_<double> &kernel = get_kernel(this->get_iterative_radius(l+0.5, (double)l));
			//fully compute the Z filter
			cv::Ptr<cv::FilterEngine> filter = createSeparableLinearFilter
			(
				this->layersG[0].type(), this->layersG[0].type(),
				kx, kernel,
				cv::Point(-1,-1), 0, cv::BORDER_DEFAULT
			);
			Image lG2D(dims[0], dims[1]*dims[2], (PixelType*)this->layersG[l].data);
			filter->apply(lG2D, layG2D);
			//now compute the XY only where needed (at the center's position)
			Image layG = cv::Mat(3, dims, lG2D.type(), layG2D.data);
			std::vector<double> gx(kernel.rows);
			ci = this->centers_no_subpix.begin();
			for(size_t c=0; c<this->centers_no_subpix.size(); ++c)
			{
				std::fill(gx.begin(), gx.end(), 0);
				const int xmin = std::max(0, (*ci)[0]-kernel.rows/2),
						xmax = std::min(layG.size[2], (*ci)[0]+kernel.rows/2+1);
				for(int y=std::max(0, (*ci)[1]-kernel.rows/2); y<std::min(layG.size[1], (*ci)[1]+kernel.rows/2+1); ++y)
					gx[y-(*ci)[1]+kernel.rows/2] = std::inner_product(
							&layG((*ci)[2], y, xmin),
							&layG((*ci)[2], y, xmin) + xmax - xmin,
							&kernel(0, xmin - (*ci)[0]+kernel.rows/2),
							0.0);
				sublayerG(c, 2*l+1) = std::inner_product(gx.begin(), gx.end(), &kernel(0,0), 0.0);
				ci++;
			}
		}
		//Estimate scale by Newton's method
		std::vector<double> dog(this->layers.size()*2-1);
		ci = this->centers_no_subpix.begin();
		for(size_t c=0; c<this->centers_no_subpix.size(); ++c)
		{
			//DoG layers values
			for(size_t u =0; u<dog.size();++u)
				dog[u] = sublayerG(c, u+2) - sublayerG(c, u);
			const int l = ci->back();
			//Apply newton's method using quadratic estimate of the derivative
			double s = l - (-dog[2*l+2] + 8*dog[2*l+1] - 8*dog[2*l-1] + dog[2*l-2])/6.0 /(dog[2*l+2]-2*dog[2*l]+dog[2*l-2]);
			//saturate the results
			if(s>l+0.5)
				s= l + 0.5;
			if(s<l-0.5)
				s= l- 0.5;
			//second round of newton's method
			if(s>=1)
			{
				if(s+0.1<l)
					//estimate near l-0.5
					s = l - 0.5 - (-dog[2*l+1] + 8*dog[2*l] - 8*dog[2*l-2] + dog[2*l-3])/6.0 /(dog[2*l+1]-2*dog[2*l-1]+dog[2*l-3]);
			}
			else
			{
				//for sizes significantly below the sampled scale
				//the linear estimate of the derivative is (marginally) better
				if(s+0.25<l)
					s = l - (dog[2*l+1] - dog[2*l-1])/(dog[2*l+2]-2*dog[2*l]+dog[2*l-2]);
			}
			if(s<l-0.5)
				s= l- 0.5;
			if(s>l+0.5)
				s= l + 0.5;
			centers[c].r = s;
			ci++;
		}
	}
};
#endif // OCTAVEFINDER_H
