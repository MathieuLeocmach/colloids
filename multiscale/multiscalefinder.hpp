/*
 * multiscalefinder.hpp
 *
 *  Created on: 11 juil. 2011
 *      Author: mathieu
 */

#include "octavefinder.hpp"
#include <boost/ptr_container/ptr_vector.hpp>

#ifndef MULTISCALEFINDER_HPP_
#define MULTISCALEFINDER_HPP_

namespace Colloids {

class MultiscaleFinder : boost::noncopyable
{
public:
	typedef OctaveFinder::Image Image;
	typedef OctaveFinder::PixelType PixelType;

	virtual ~MultiscaleFinder();

	//accessors
	inline const size_t get_n_octaves() const {return this->octaves.size();};
	inline const OctaveFinder & get_octave(const size_t l) const {return *this->octaves[l];};
	virtual const size_t get_width() const =0;
	const size_t get_height() const {return this->octaves[0]->get_height()/2; };
	inline const size_t get_n_layers() const {return this->octaves[0]->get_n_layers();};
	const double & get_radius_preblur() const {return this->octaves[0]->get_radius_preblur();}
	const double & get_prefactor() const {return this->octaves[0]->get_prefactor();}
	void set_radius_preblur(const double &k=1.6);
	//processing
	void fill(const cv::Mat &input);
	void initialize_binary();
	template<int D>	inline void subpix(std::vector<Center<D> >& centers) const;
	template<int D>
	inline void get_centers(const cv::Mat &input, std::vector<Center<D> >& centers);
	template<int D> const std::vector<int> previous_octave_coords(const Center<D> &v) const;
	virtual Image downscale(const size_t &o) const = 0;
	virtual Image upscale(const cv::Mat &input) const;
	template<int D> inline void seam(Center<D> &v, const size_t &o) const{};


protected:
	std::vector<OctaveFinder*> octaves;
	//Image small, upscaled;
	MultiscaleFinder(){};
};

class MultiscaleFinder2D : public MultiscaleFinder
{
public:
	MultiscaleFinder2D(const int nrows=256, const int ncols=256, const int nbLayers=3, const double &preblur_radius=1.6);
	virtual const size_t get_width() const {return this->octaves[0]->get_width()/2; };
	virtual Image downscale(const size_t &o) const;
};

class MultiscaleFinder1D : public MultiscaleFinder
{
public:
	MultiscaleFinder1D(const int ncols=256, const int nbLayers=3, const double &preblur_radius=1.6);
	virtual const size_t get_width() const {return 1; };
	virtual Image downscale(const size_t &o) const;
};

class MultiscaleFinder3D : public MultiscaleFinder
{
public:
	MultiscaleFinder3D(const int nplanes=256, const int nrows=256, const int ncols=256, const int nbLayers=3, const double &preblur_radius=1.6);
	virtual const size_t get_width() const {return this->octaves[0]->get_width()/2; };
	const size_t get_depth() const {return dynamic_cast<OctaveFinder3D*>(this->octaves[0])->get_depth()/2; };
	virtual Image downscale(const size_t &o) const;
	virtual Image upscale(const cv::Mat &input) const;
};

/**
 * \brief Locate centers with subpixel and subscale resolution
 */
template<int D>
void MultiscaleFinder::subpix(std::vector<Center<D> > &centers) const
{
	centers.clear();
	//reserve memory for the center container
	size_t n_centers = 0;
	for(size_t o=0; o<this->octaves.size(); ++o)
		n_centers += this->octaves[o]->get_nb_centers();
	centers.reserve(n_centers);
	//subpixel resolution
	for(size_t o=0; o<this->octaves.size(); ++o)
	{
		std::vector<Center<D> > v;
		this->octaves[o]->subpix(v);
		//correct the seam between octaves in sizing precision
		if(o>0)
			for(size_t p=0; p<v.size(); ++p)
				this->seam(v[p], o-1);
		//transform scale coordinate in size coordinate
		for(size_t c=0; c< v.size(); ++c)
			this->octaves[o]->scale(v[c]);
		//stack up
		for(size_t c=0; c<v.size(); ++c)
		{
			for(int d=0; d<D; ++d)
			v[c][d] *= pow(2.0, (int)(o)-1);
			v[c].r *= pow(2.0, (int)(o)-1);
			centers.push_back(v[c]);
		}
	}
}

/**
 * Efficient processing pipeline from the input image to the output centers
 */
template<int D>
inline void MultiscaleFinder::get_centers(const cv::Mat & input, std::vector<Center<D> >& centers)
{
	this->fill(input);
	this->initialize_binary();
	this->subpix(centers);
}

template<int D>
const std::vector<int> MultiscaleFinder::previous_octave_coords(const Center<D> &v) const
{
	const double n = this->get_n_layers();
	std::vector<int> vi(D+1);
	for(int u=0; u<D; ++u)
		vi[u] = (int)(v[u]*2+0.5);
	vi[D] = int(v.r + n + 0.5);
	return vi;
}
template<>
inline void MultiscaleFinder::seam(Center2D &v, const size_t &o) const
{
	if(v.r<1)
	{
		std::vector<int> ci = this->previous_octave_coords(v);
		v.r = this->octaves[o]->scale_subpix(ci) - this->get_n_layers();
	}
}
}

#endif /* MULTISCALEFINDER_HPP_ */
