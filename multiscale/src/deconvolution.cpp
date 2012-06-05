/*
 * deconvolution.cpp
 *
 *  Created on: 4 juin 2012
 *      Author: mathieu
 */

#include "deconvolution.hpp"
#include <algorithm>

namespace Colloids {

	Convolver::Convolver(unsigned long int size) : _size(size), _fourier_size(size/2+1)
	{
		this->real = fftwf_alloc_real(this->_size);
		this->fourier = fftwf_alloc_complex(this->_fourier_size);
		this->forward = fftwf_plan_dft_r2c_1d(this->_size, this->real, this->fourier,
                FFTW_MEASURE);
		this->backward = fftwf_plan_dft_c2r_1d(this->_size, this->fourier, this->real,
                FFTW_MEASURE);
	}

	Convolver::~Convolver()
	{
		fftwf_destroy_plan(this->forward);
		fftwf_destroy_plan(this->backward);
		fftwf_free(real); fftwf_free(fourier);
	}

	void Convolver::spectrum(const float *input, int step, float *output)
	{
		this->fill(input, step);
		fftwf_execute(this->forward);
		for (unsigned long int i=0; i<this->_fourier_size; ++i)
			*output++ = std::norm(*reinterpret_cast<std::complex<float>*>(this->fourier+i));
	}

	void Convolver::operator()(float* input, const int step, const float* kernel)
	{
		this->fill(input, step);
		fftwf_execute(this->forward);
		for(unsigned long int i=0; i<this->_fourier_size; ++i)
		{
			const float factor = *kernel++ / this->_size;
			this->fourier[i][0] *= factor;
			this->fourier[i][1] *= factor;
		}
		fftwf_execute(this->backward);
		for (unsigned long int i=0; i<this->_size; ++i)
			*input++ = this->real[i];
	}
	void Convolver::fill(const float* input, const int step)
	{
		for (unsigned long int i=0; i<this->_size; ++i)
		{
			this->real[i] = *input;
			input += step;
		}
	}

}
