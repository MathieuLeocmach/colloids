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
		this->real = (float *) fftwf_malloc(sizeof(float) * this->_size);
		this->fourier = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * this->_fourier_size);
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

	std::vector<float> get_spectrum_1d(const cv::Mat_<float> &im, const int axis)
	{
		if(axis >= im.dims)
			throw std::invalid_argument("Matrix dimension is too small to compute the spectrum along this axis");
		assert(im.isContinuous());
		Convolver co(im.size[axis]);
		std::vector<float> spectrum(co.fourier_size());
		std::vector<double> tot(co.fourier_size(), 0.0);
		unsigned long int step = im.step1(axis);
		//whatever the real dimension, we fall back to a 3d situation where the axis of interest is y
		//and either x or z can be of size 1
		int nbplanes = 1;
		for(int d=0; d<axis; ++d)
			nbplanes *= im.size[d];
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

	std::vector<float> get_deconv_kernel(const cv::Mat_<float> &im, const int good_axis, const int bad_axis, const double size_ratio)
	{
		std::vector<float> bad_sp = get_spectrum_1d(im, bad_axis);
		std::vector<float> good_sp = get_spectrum_1d(im, good_axis);
		//linear interpolation of good_sp to take into account the voxel size ratio
		const double qratio = bad_sp.size() * size_ratio / good_sp.size();
		std::vector<float> scaled(std::min(bad_sp.size(), (size_t)((good_sp.size()-1)/qratio)), 0);
		for(size_t i=0; i<scaled.size(); ++i)
		{
			const double q = i * qratio;
			const size_t nq = q;
			const double dq = q-nq;
			scaled[i] = (1.0-dq) * good_sp[nq] + dq * good_sp[nq+1];
		}
		//deconvolution kernel
		std::vector<float> kernel(bad_sp.size(), 1.0f);
		for(size_t i=0; i<bad_sp.size() && i<scaled.size(); ++i)
			if(1.0f+bad_sp[i]*bad_sp[i]>1.0f)
				kernel[i] = sqrt(scaled[i]/bad_sp[i]);
		return kernel;
	}

	void convolve(cv::Mat_<float> &im, const int axis, const float* kernel)
	{
		if(axis >= im.dims)
			throw std::invalid_argument("Matrix dimension is too small to convolve along this axis");
		assert(im.isContinuous());
		Convolver co(im.size[axis]);
		unsigned long int step = im.step1(axis);
		//whatever the real dimension, we fall back to a 3d situation where the axis of interest is y
		//and either x or z can be of size 1
		int nbplanes = 1;
		for(int d=0; d<axis; ++d)
			nbplanes *= im.size[d];
		int planestep = im.total()/nbplanes;
		//for each plane
		for(int i=0; i<nbplanes; ++i)
			//for each line
			for(size_t j=0; j<step; ++j)
				co(reinterpret_cast<float*>(im.data) + i*planestep + j, step, kernel);
	}
}
