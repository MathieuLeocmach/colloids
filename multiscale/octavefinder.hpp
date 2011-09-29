#ifndef OCTAVEFINDER_H
#define OCTAVEFINDER_H

#include "center.h"
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <list>

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
			static const cv::Mat_<double>& get_kernel(const double &sigma);

            //processing
            void fill(const cv::Mat &input);
            void preblur_and_fill(const cv::Mat &input);
            virtual void initialize_binary(const double &max_ratio = 1.2);
            void subpix(std::vector<Center2D>& centers) const;
            void single_subpix(const cv::Vec3i & ci, Center2D &c) const;
            virtual Center2D spatial_subpix(const cv::Vec3i & ci) const;
            virtual double scale_subpix(const cv::Vec3i & ci) const;
            void scale(Center2D &c) const{
            	c.r = this->prefactor * this->preblur_radius * pow(2.0, (c.r + 1) / this->get_n_layers());
            };
            std::vector<Center2D> operator()(const cv::Mat &input, const bool preblur=false);
            virtual double gaussianResponse(const size_t &j, const size_t &i, const double & scale) const;
            void seam_binary(OctaveFinder & other);


    protected:
            std::vector<Image > layersG, layers;
            std::vector<cv::Mat_<bool> > binary;
            std::vector<double> iterative_radii;
            std::vector<cv::FilterEngine> iterative_gaussian_filters;
            std::vector<size_t> sizes;
            std::list<cv::Vec3i> centers_no_subpix;
            double preblur_radius, prefactor;
            static std::map<size_t, cv::Mat_<double> > kernels;
            cv::Ptr<cv::FilterEngine> preblur_filter;

            virtual void _fill_internal();
            virtual void preblur(const cv::Mat &input);
            virtual void fill_iterative_radii();
    };

    class OctaveFinder1D : public OctaveFinder
    {
		public:
			OctaveFinder1D(const int ncols=256, const int nbLayers=3, const double &preblur_radius=1.6) :
				OctaveFinder(1, ncols, nbLayers, preblur_radius){};

			virtual void initialize_binary(const double &max_ratio = 1.1);
			virtual Center2D spatial_subpix(const cv::Vec3i & ci) const;
			virtual double scale_subpix(const cv::Vec3i & ci) const;
			virtual double gaussianResponse(const size_t &j, const size_t &i, const double & scale) const;
    };

    class OctaveFinder3D : public OctaveFinder
    {
		public:
			OctaveFinder3D(const int nplanes=256, const int nrows=256, const int ncols=256, const int nbLayers=3, const double &preblur_radius=1.6);
			virtual ~OctaveFinder3D();

			inline const int & get_depth() const {return this->layers[0].size[this->layers[0].dims-3];};

			/*virtual void initialize_binary(const double &max_ratio = 1.1);
			virtual Center2D spatial_subpix(const cv::Vec3i & ci) const;
			virtual double scale_subpix(const cv::Vec3i & ci) const;
			virtual double gaussianResponse(const size_t &j, const size_t &i, const double & scale) const;*/

		protected:
			boost::iostreams::mapped_file file;
			std::string path;
			std::vector<Image > layersG2D;
			std::vector<cv::FilterEngine> iterative_Zgaussian_filters;
			cv::Ptr<cv::FilterEngine> preblur_Zfilter;

			virtual void fill_iterative_radii();
			virtual void preblur(const cv::Mat &input);
			virtual void _fill_internal();
    };
};
#endif // OCTAVEFINDER_H
