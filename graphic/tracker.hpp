/**
 * \file tracker.hpp
 * \brief Define classes to track particles in 3D image series
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 7 april 2009
 *
 * Object oriented implementation of the code elaborated in 2008.
 * Based on Crocker & Grier algorithm
 *
 */

#include "lifFile.hpp"
#include "../files_series.hpp"
#include "graphicParticles.hpp"
#include <CImg.h>

/** \brief basic tracker class containing the tracking algorythm*/
class Tracker
{
    void FFTapplyMask();

    public:
        cimg_library::CImg<float> image;
        cimg_library::CImg<bool> FFTmask;
        fftwf_complex *data;
        fftwf_plan forward_plan,backward_plan;

        size_t channel;
        bool view,quiet;
        double displayRadius;
        double Zratio;

        void load3D(const size_t &time);
        void makeBandPassMask(const double &radiusMin, const double &radiusMax, const double &zradiusMin, const double &zradiusMax);
        GraphicParticles trackXYZ(const double &threshold=0);
        std::vector<GraphicParticles*> granulometry(const double &radiusMin, const double &radiusMax);
		static void saveWisdom(const std::string & filename);
		static bool loadWisdom(const std::string & filename);
};

/** \brief Tracker able to import data from Leica LIF file*/
class LifTracker : public Tracker
{
    public:
        std::auto_ptr<LifFile> lif;
        size_t serie;

        explicit LifTracker(std::auto_ptr_ref<LifFile> l, const size_t &s,const size_t &ch=0, unsigned flags=FFTW_ESTIMATE);
        void load3D(const size_t &time);
        ~LifTracker();
};

/** \brief Tracker able to import data from a serie of 2D image files */
class FileSerieTracker : public Tracker
{
    public:
        TokenTree *fileSerie;

        explicit FileSerieTracker(TokenTree &tt,const double &Zra,const size_t &xsize, const size_t &ysize,const size_t &zsize,const size_t &ch, unsigned flags=FFTW_ESTIMATE);

        void load3D(const size_t &time);
};
