/*
 * liffinder.h
 *
 *  Created on: 11 août 2011
 *      Author: mathieu
 */

#ifndef LIFFINDER_H_
#define LIFFINDER_H_

#include "lifFile.hpp"
#include "reconstructor.h"
#include "multiscalefinder.hpp"

namespace Colloids {

class LocatorFromLif {
public:
	typedef Reconstructor::OutputType Centers;
	explicit LocatorFromLif(LifSerie * serie);
	virtual ~LocatorFromLif();

	//accessors
	const Reconstructor & get_reconstructor() const {return this->rec;};
	const cv::Mat_<unsigned char> & get_slice() const {return this->slice;};
	const size_t get_z() const {return get_reconstructor().size();}
	const size_t size() const {return total_t;}
	const size_t & get_t() const {return t;}

	//processing
	void clear();
	void fill_next_slice();
	void fill_time_step();
	void get_centers(Centers &centers);

private:
	LifSerie * serie;
	std::auto_ptr<MultiscaleFinder2D> finder;
	Reconstructor rec;
	cv::Mat_<unsigned char> slice;
	std::vector<size_t> dims;
	size_t t, total_t, total_z;
	std::istreambuf_iterator<char> input;
};

}

#endif /* LIFFINDER_H_ */
