/*
 * deconvolution.hpp
 *
 *  Created on: 13 avr. 2012
 *      Author: mathieu
 */

#ifndef DECONVOLUTION_HPP_
#define DECONVOLUTION_HPP_

#include<vector>
#include<complex>
#include<fftw3.h>

namespace Colloids {

	class Convolver
	{
	public:
		//constructor
		Convolver(int size, int step);

		//accessors
		int size(){return this->_size;}
		int fourier_size(){return this->_fourier_size;}

		//processing
		void spectrum(const float* input, float* output);
		void operator()(float* input, float* kernel);

	protected:
		int _size;
		int _fourier_size;
		fftwf_complex* fourier;
		fftw_plan forward, backward;
	};

}

#endif /* DECONVOLUTION_HPP_ */
