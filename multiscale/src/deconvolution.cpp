/*
 * deconvolution.cpp
 *
 *  Created on: 4 juin 2012
 *      Author: mathieu
 */

#include "deconvolution.hpp"
#include <algorithm>
#include <stdexcept>

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

	std::vector<float> get_spectrum_1d(const cv::Mat_<float> &im, int axis)
	{
		if(axis >= im.dims)
			throw std::invalid_argument("Matrix dimension is too small to compute the spectrum along this axis");
		assert(im.isContinuous());
		Convolver co(im.size[axis]);
		std::vector<float> spectrum(co.size());
		std::vector<double> tot(co.size(), 0.0);
		unsigned long int step = im.step1(axis);
		//whatever the real dimension, we fall back to a 3d situation where the dimension of interest is y
		//and either x or z can be of size 1
		int nbplanes = 1;
		for(int d=0; d<axis; ++d)
			nbplanes *= im.size[axis];
		int planestep = im.total()/nbplanes;
		//for each plane
		for(int i=0; i<nbplanes; ++i)
		{
			//for each line
			for(size_t j=0; j<step; ++j)
			{
				co.spectrum(reinterpret_cast<float* const>(im.data) + i*planestep + j, step, &spectrum[0]);
				for(size_t u=0; u<spectrum.size(); ++u)
					tot[u] += spectrum[u];
			}
		}
		const double icount = 1.0 / (nbplanes * step);
		for(size_t i=0; i<tot.size(); ++i)
			spectrum[i] = tot[i]*icount;
		return spectrum;
	}

}
