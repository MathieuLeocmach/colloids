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

	MultiscaleFinder(const int nrows=256, const int ncols=256, const int nbLayers=3, const double &preblur_radius=1.6);
	virtual ~MultiscaleFinder();

	//accessors
	inline const size_t get_n_octaves() const {return this->octaves.size();};
	inline const OctaveFinder & get_octave(const size_t l) const {return *this->octaves[l];};
	inline const size_t get_width() const {return this->octaves[0]->get_width()/2; };
	inline const size_t get_height() const {return this->octaves[0]->get_height()/2; };
	inline const size_t get_n_layers() const {return this->octaves[0]->get_n_layers();};
	const double & get_radius_preblur() const {return this->octaves[0]->get_radius_preblur();}
	void set_radius_preblur(const double &k=1.6);
	std::vector<cv::Vec4d> operator()(const cv::Mat &input);


private:
	std::vector<OctaveFinder*> octaves;
	cv::Mat_<double> small, upscaled;

};

}

#endif /* MULTISCALEFINDER_HPP_ */
