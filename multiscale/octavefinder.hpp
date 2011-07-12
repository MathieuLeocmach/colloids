#ifndef OCTAVEFINDER_H
#define OCTAVEFINDER_H

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <boost/ptr_container/ptr_vector.hpp>
#include <list>

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
            const double & get_iterative_radius(const size_t l) const {return this->iterative_radii[l];};
            const size_t & get_size(const size_t l) const {return this->sizes[l];};
            inline const cv::Mat_<bool> & get_binary(const size_t l) const {return binary[l-1];};
			inline const cv::Mat_<double> & get_layersG(const size_t l) const {return layersG[l];}
			inline const cv::Mat_<double> & get_layers(const size_t l) const {return layers[l];}

            //processing
            void fill(const cv::Mat &input);
            void preblur_and_fill(const cv::Mat &input);
            void initialize_binary(const double &max_ratio = 1.1);
            std::vector<cv::Vec4d> subpix() const;
            cv::Vec4d single_subpix(const cv::Vec3i & ci) const;
            void scale(std::vector<cv::Vec4d> &centers) const;
            std::vector<cv::Vec4d> operator()(const cv::Mat &input, const bool preblur=false);


    protected:
        private:
            std::vector<cv::Mat_<double> > layersG, layers;
            std::vector<cv::Mat_<bool> > binary;
            std::vector<double> iterative_radii;
            std::vector<cv::FilterEngine> iterative_gaussian_filters;
            std::vector<size_t> sizes;
            std::list<cv::Vec3i> centers_no_subpix;
            double preblur_radius;

            void fill_iterative_radii(const double &k=1.6);
    };
};
#endif // OCTAVEFINDER_H
