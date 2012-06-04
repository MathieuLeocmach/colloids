/*
 * deconvolution.cpp
 *
 *  Created on: 4 juin 2012
 *      Author: mathieu
 */

#include "deconvolution.hpp"

namespace Colloids {

	Convolver::Convolver(int size, int step) : _size(size), _fourier_size(size/2+1)
	{
		this->fourier = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * this->_fourier_size);
	}

}
