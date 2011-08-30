#ifndef OCTAVEFINDER_GPU_H
#define OCTAVEFINDER_GPU_H

#include "octavefinder.hpp"
#include <opencv2/gpu/gpu.hpp>

namespace Colloids
{

    class OctaveFinder_GPU : OctaveFinder
    {
        public:
			typedef float PixelType;
			typedef cv::Mat_<float> Image;

            OctaveFinder_GPU(const int nrows=256, const int ncols=256, const int nbLayers=3, const double &preblur_radius=1.6);
            virtual ~OctaveFinder_GPU();

            //accessors
			inline const bool use_gpu() const {return !gpu_layersG.empty();}
			inline const std::vector<int>& get_n_iter() const {return n_iter;}
			virtual void set_radius_preblur(const double &k=1.6);

            //processing
			virtual void fill(const cv::Mat &input);
			virtual void initialize_binary(const double &max_ratio = 1.1);

    protected:
            std::vector<cv::gpu::GpuMat> gpu_layersG, gpu_layers;
            std::vector<cv::Ptr<cv::gpu::FilterEngine_GPU> > gpu_iterative_gaussian_filters;
            std::vector<int> n_iter;
            cv::gpu::GpuMat tmp;
            cv::gpu::Stream str;

            void initialize_gpu();
    };
};
#endif // OCTAVEFINDER_GPU_H
