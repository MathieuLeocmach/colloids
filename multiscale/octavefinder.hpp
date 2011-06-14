#ifndef OCTAVEFINDER_H
#define OCTAVEFINDER_H

#include <opencv/cv.h>
#include <opencv/highgui.h>

namespace Colloids
{

    class OctaveFinder
    {
        public:
            const int width, height, n_layers;
            OctaveFinder(int nrows=256, int ncols=256, int nbLayers=3);
            virtual ~OctaveFinder();
        protected:
        private:
            cv::Mat layersG, layers;
    };
};
#endif // OCTAVEFINDER_H
