/**
 * \file graphicParticles.hpp
 * \brief Defines classes for particles interacting with pictures
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 26 Novemeber 2008
 *
 * Define functions and classes relative to the particles for the particle tracking code
 *
 */

#ifndef graphic_particles_H
#define graphic_particles_H

#include "../particles.hpp"
#include <CImg.h>
#include <list>

/**
    \brief defines a set of particles in relation with CImg
    Drawing and display functions. Results from tracking.
*/
class GraphicParticles : public Particles
{
    public:
        GraphicParticles(const cimg_library::CImg<bool> &img, const double &rad =0);
        GraphicParticles(const std::string &filename, const double &rad);

        template<typename T>
        void drawSpheres(cimg_library::CImg<T> &img,const T &val) const;

        template<typename T>
        void displaySpheres(const cimg_library::CImg<T> &img,const T &val) const;

        template<typename T>
        cimg_library::CImg<T> getRepresentation(const T &val,const T &bg) const;

        template<typename T>
        cimg_library::CImg<T> get_neigbourhood(const size_t p,const cimg_library::CImg<T> &img,const double &rad) const;

        template<typename T>
        void sortDecreasingIntensity(const cimg_library::CImg<T> &img);

        template<typename T>
        void centroid(const cimg_library::CImg<T> &img);

        void clean(void);
};

/**
    \brief Makes a bounding box out of the dimensions of an image
*/
template <typename T>
BoundingBox imageBounds(const cimg_library::CImg<T> &img,const double &Zratio=1.0)
{
	BoundingBox bb;

	for(size_t i=0;i<3;++i)
        bb.edges[i].first  = 0;

	bb.edges[0].second  = img.dimx()-1;
	bb.edges[1].second  = img.dimy()-1;
	bb.edges[2].second  = (img.dimz()-1)*Zratio;

	return bb;
}

/**
    \brief compare the intensities of two centers
*/
inline bool compIntensities(const std::pair< float,std::valarray<double> > &c,const std::pair< float,std::valarray<double> > &d)
{
        return c.first<d.first;
}

/**
    \brief draw a sphere.
*/
template <typename T>
void drawSphere(cimg_library::CImg<T> &img, const int x, const int y, const double z, const double r,const T val)
{
    const T value[] = {val};
    const float radTwo = pow(r,2.0);

    if(img.dimz()>1)
        for(int k=(int)(std::max(z-r,0.0));k<z+r && k<img.dimz();++k)
            img.get_shared_plane(k).draw_circle(x,y,(int)sqrt(radTwo-pow((double)k-z,2.0)),value);
    else
        img.draw_circle(x,y,(int)r,value);
}
/**
    \brief draw an ellipsoid.
*/
template <typename T>
void drawEllipsoid(cimg_library::CImg<T> &img, const int x, const int y, const double z, const double r, const double rz,const T val)
{
    const T value[] = {val};
    const float radzTwo = pow(rz,2.0);

    if(img.dimz()>1)
        for(int k=(int)(std::max(z-rz,0.0));k<z+rz && k<img.dimz();++k)
            img.get_shared_plane(k).draw_circle(x,y,(int)(r*sqrt(1.0-(k-z)*(k-z)/radzTwo)),value);
    else
        img.draw_circle(x,y,(int)r,value);
}

/**
    \brief draw a sphere representing a center. Can be used to erase a center's neigbourhood.
*/
template <typename T>
void drawSphere(const std::valarray<double> &Center,double rad,cimg_library::CImg<T> &img,const T &val)
{
    drawSphere(img,(int)Center[0],(int)Center[1],Center[2],rad,val);
}

/**
    \brief draw the particles as spheres on an image.
    \param img Image to draw to
    \param val Color of the circle
*/
template <typename T>
void GraphicParticles::drawSpheres(cimg_library::CImg<T> &img,const T &val) const
{
    for(const_iterator p=begin();p!=end();++p)
        drawSphere((*p),radius,img,val);
}

/**
    \brief Display found centers superimposed to the experimental image.
    \param img Experimental image
    \param val value of the center color
*/
template <typename T>
void GraphicParticles::displaySpheres(const cimg_library::CImg<T> &img,const T &val) const
{
    cimg_library::CImg<T> final(img.width,img.height,img.depth,2,0);

    drawSpheres(final,val);
    final.get_shared_channel(1)=img;
    final.display("image+spheres");
}

template<typename T>
cimg_library::CImg<T> GraphicParticles::getRepresentation(const T &val,const T &bg) const
{
    cimg_library::CImg<T> img((size_t)bb.edges[0].second,(size_t)bb.edges[1].second,(size_t)bb.edges[2].second,1,bg);
    drawSpheres(img,val);
    return img;
}

/**
    \brief get the pixels near to a particle
    \param p Index of the particle
    \param img image to crop from
    \param rad radius of the neigbourhood
*/
template <typename T>
cimg_library::CImg<T> GraphicParticles::get_neigbourhood(const size_t p,const cimg_library::CImg<T> &img,const double &rad) const
{
    if(img.dimz()>2*rad+1)
        return img.get_crop((int)(at(p)[0]-rad),(int)(at(p)[1]-rad),(int)(at(p)[2]-rad),(int)(at(p)[0]+rad),(int)(at(p)[1]+rad),(int)(at(p)[2]+rad),false);
    return img.get_crop((int)(at(p)[0]-rad),(int)(at(p)[1]-rad),(int)(at(p)[2]),(int)(at(p)[0]+rad),(int)(at(p)[1]+rad),(int)(at(p)[2]),false);
}

/**
    \brief sort the particles by decreasing neigbourhood intensity
*/
template<typename T>
void GraphicParticles::sortDecreasingIntensity(const cimg_library::CImg<T> &img)
{
    std::list< std::pair< double,std::valarray<double> > > li(0);
    for(size_t p=0;p<size();++p)
        li.push_back(std::make_pair(-get_neigbourhood(p,img,1.0).sum(), at(p)));

    li.sort(compIntensities);
    iterator p=begin();
    for(std::list< std::pair< double, std::valarray<double> > >::iterator l = li.begin(); l!=li.end(); ++l)
        (*p++)=(*l).second;

}

/** \brief sub pixel resolution */
template<typename T>
void GraphicParticles::centroid(const cimg_library::CImg<T> &img)
{
    cimg_library::CImg<double> Ix;
    cimg_library::CImgList<double> slide;
    std::vector< std::vector<double> > disp(3, std::vector<double>(size(),0.0));
    std::vector<double> meanDisp(3,0.0);
    cimg_library::CImg<T> sample;
    if(img.dimz()>2)
    {
        Ix.assign(3,3,3).fill(-1.0,0.0,1.0);
        slide = Ix << Ix.get_permute_axes("yxzv") << Ix.get_permute_axes("zyxv");
        sample.assign(3,3,3);
    }
    else
    {
        Ix.assign(3,3).fill(-1.0,0.0,1.0);
        slide = Ix << Ix.get_permute_axes("yxzv");
        sample.assign(3,3);
    }

    for(size_t c=0;c<size();++c)
    {
        sample.assign(get_neigbourhood(c,img,1.0));
        //const float M = sample.max(), m=sample.min();
        const double norm = sample(1,1); //3*(M-m)/(4*M+5*m)*sample.sum();
        for(size_t i=0;i<slide.size;++i)
        {
            disp[i][c] = sample.dot(slide[i])/norm;
            meanDisp[i] += std::abs(disp[i][c]);
        }
    }
    const double nor = 0.25*size();
    for(size_t i=0;i<slide.size;++i)
        meanDisp[i] /= nor;

    for(size_t c=0;c<size();++c)
        for(size_t i=0;i<slide.size;++i)
            at(c)[i] += disp[i][c]/meanDisp[i];
}
#endif

