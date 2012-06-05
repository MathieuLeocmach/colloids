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
		//constructor, destructor
		Convolver(unsigned long int size);
		~Convolver();

		//accessors
		int size(){return this->_size;}
		int fourier_size(){return this->_fourier_size;}

		//processing
		void spectrum(const float* input, const int step, float* output);
		void operator()(float* input, const int step, float* kernel);

	protected:
		unsigned long int _size;
		unsigned long int _fourier_size;
		float* real;
		fftwf_complex* fourier;
		fftwf_plan forward, backward;
	};

}

#endif /* DECONVOLUTION_HPP_ */
