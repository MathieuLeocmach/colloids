#include "octavefinder.hpp"

using namespace cv;
using namespace std;
using namespace Colloids;

OctaveFinder::OctaveFinder(int nrows, int ncols, int nbLayers) :
    width(nrows), height(ncols), n_layers(nbLayers)
{
    int sizesG[] = {nbLayers+3, nrows, ncols};
    int sizes[] = {nbLayers+2, nrows, ncols};
    this->layersG = Mat(3, *sizesG, CV_64FC1);
    this->layersG = Mat(3, *sizes, CV_64FC1);
}

OctaveFinder::~OctaveFinder()
{
    //dtor
}
