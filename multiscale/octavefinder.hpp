#ifndef OCTAVEFINDER_H
#define OCTAVEFINDER_H

#include <opencv/cv.h>
#include <opencv/highgui.h>

namespace Colloids
{

    class OctaveFinder
    {
        public:
            OctaveFinder(const int nrows=256, const int ncols=256, const int nbLayers=3, const double &k=1.6);
            virtual ~OctaveFinder();

            //accessors
            inline const int & get_width() const {return this->layers[0].rows; };
            inline const int & get_height() const {return this->layers[0].cols; };
            inline const int get_n_layers() const {return this->layers.size()-2;};
            const double & get_radius_preblur() const {return this->k;}
            void set_radius_preblur(const double &k=1.6);
            const double & get_iterative_radius(const size_t l) const {return this->iterative_radii[l];};
            const int & get_size(const size_t l) const {return this->sizes[l];};
            inline const cv::Mat_<bool> & get_binary(const size_t l) const {return binary[l];};
			inline const cv::Mat_<double> & get_layersG(const size_t l) const {return layersG[l];}
			inline const cv::Mat_<double> & get_layers(const size_t l) const {return layers[l];}

            //processing
            void fill(const cv::Mat &input);
            void preblur_and_fill(const cv::Mat &input);
            void initialize_binary();


    protected:
        private:
            std::vector<cv::Mat_<double> > layersG, layers;
            std::vector<cv::Mat_<bool> > binary;
            std::vector<double> iterative_radii;
            std::vector<cv::FilterEngine> iterative_gaussian_filters;
            std::vector<int> sizes;
            double k;
            std::vector<cv::Mat_<double> > block_min_val;
            std::vector< cv::Mat_<cv::Vec3i> > block_min_pos;

            void fill_iterative_radii(const double &k=1.6);
    };
};
#endif // OCTAVEFINDER_H
