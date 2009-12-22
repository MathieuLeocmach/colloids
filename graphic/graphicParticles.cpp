#include "graphicParticles.hpp"

using namespace std;

/** \brief constructor out of datafile */

GraphicParticles::GraphicParticles(const string &filename, const double &rad) : Particles(filename)
{
    radius = rad;
}

/**
    \brief constructor out of the coordinates of the non-zero pixels of a binary image
    \param img The image with some non-zero pixels indication the centers' positions
    \param rad Real radius of the particles and also the extend of the margin not to be searched in
*/
GraphicParticles::GraphicParticles(const cimg_library::CImg<bool> &img, const double &rad) : Particles(0,0.0)
{
    radius = rad;
    bb = imageBounds(img);

    int margin = (int)rad, marginZ;
    if(img.dimz()<2+2*margin) marginZ = 0;
    else marginZ = margin;

    for(int i=margin;i+margin<img.dimx();++i)
      for(int j=margin;j+margin<img.dimy();++j)
        for(int k=marginZ;k+marginZ<img.dimz();++k)
          if(img(i,j,k,0)!=0)
          {
            double a[] = {(double)i,(double)j,(double)k};
            push_back(valarray<double>(a,3));
          }
    return;
}


/**
    \brief remove particles out of the bounding box
*/
void GraphicParticles::clean(void)
{
    if(bb.edges[2].second-bb.edges[2].first>2*radius+1)
    {
         deque< valarray<double> > newCenters;
         for(iterator p=begin();p!=end();++p)
           if(bb.overlaps(bounds(*p,radius)))
             newCenters.push_back(*p);

         assign(newCenters.begin(),newCenters.end());
    }
}
