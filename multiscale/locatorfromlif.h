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
	explicit LocatorFromLif(LifSerie * serie);
	virtual ~LocatorFromLif();

	//accessors
	const Reconstructor & get_reconstructor() const {return this->rec;};
	const cv::Mat_<unsigned char> & get_slice() const {return this->slice;};

	//processing
	void fill_next_slice();

private:
	LifSerie * serie;
	std::auto_ptr<MultiscaleFinder2D> finder;
	Reconstructor rec;
	cv::Mat_<unsigned char> slice;
	std::vector<size_t> dims;
	size_t t, total_t, z;
	std::istreambuf_iterator<char> input;
};

}

#endif /* LIFFINDER_H_ */
