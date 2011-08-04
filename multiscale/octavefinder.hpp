#ifndef OCTAVEFINDER_H
#define OCTAVEFINDER_H

#include "center.h"
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <boost/ptr_container/ptr_vector.hpp>
#include <list>
#include <memory>

namespace Colloids
{

    class OctaveFinder : boost::noncopyable
    {
        public:
            OctaveFinder(const int nrows=256, const int ncols=256, const int nbLayers=3, const double &preblur_radius=1.6);
            virtual ~OctaveFinder();

            //accessors
            inline const int & get_width() const {return this->layers[0].rows; };
            inline const int & get_height() const {return this->layers[0].cols; };
            inline const size_t get_n_layers() const {return this->layers.size()-2;};
            const double & get_radius_preblur() const {return this->preblur_radius;}
            void set_radius_preblur(const double &k=1.6);
            const double & get_prefactor() const {return this->prefactor;}
            const double & get_iterative_radius(const size_t l) const {return this->iterative_radii[l];};
            const double get_iterative_radius(const double &larger, const double &smaller) const;
            const size_t & get_size(const size_t l) const {return this->sizes[l];};
            inline const cv::Mat_<bool> & get_binary(const size_t l) const {return binary[l-1];};
			inline const cv::Mat_<double> & get_layersG(const size_t l) const {return layersG[l];}
			inline const cv::Mat_<double> & get_layers(const size_t l) const {return layers[l];}
			inline const size_t get_nb_centers() const {return centers_no_subpix.size();}

            //processing
            void fill(const cv::Mat &input);
            void preblur_and_fill(const cv::Mat &input);
            virtual void initialize_binary(const double &max_ratio = 1.1);
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
            std::vector<cv::Mat_<double> > layersG, layers;
            std::vector<cv::Mat_<bool> > binary;
            std::vector<double> iterative_radii;
            std::vector<cv::FilterEngine> iterative_gaussian_filters;
            std::vector<size_t> sizes;
            std::list<cv::Vec3i> centers_no_subpix;
            double preblur_radius, prefactor;

            void fill_iterative_radii(const double &k=1.6);
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
};
#endif // OCTAVEFINDER_H
