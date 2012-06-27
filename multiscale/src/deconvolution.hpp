/*
 * deconvolution.hpp
 *
 *  Created on: 13 avr. 2012
 *      Author: mathieu
 */

#ifndef DECONVOLUTION_HPP_
#define DECONVOLUTION_HPP_

#include <opencv2/core/core.hpp>
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
		const fftwf_complex* get_fourier(){return this->fourier;}

		//processing
		/**
		 * \brief Compute the spectrum of a (possibly discontinuous) input
		 */
		void spectrum(const float* input, const int step, float* output);
		/**
		 * \brief Convolve in place the (possibly discontinuous) input with the kernel (given in Fourier space)
		 */
		void operator()(float* input, const int step, const float* kernel);
		/**
		 * \brief set the window function to Hanning
		 */
		void set_hanning();
		void unset_window(){this->window.clear();}
		bool windowing(){return !this->window.empty();}

	protected:
		unsigned long int _size;
		unsigned long int _fourier_size;
		float* real;
		fftwf_complex* fourier;
		fftwf_plan forward, backward;
		std::vector<double> window;
		void fill(const float* input, const int step);
	};

	std::vector<float> get_spectrum_1d(const cv::Mat_<float> &im, const int axis=0, const bool windowing=true);
	std::vector<float> get_deconv_kernel(const cv::Mat_<float> &im, const int good_axis, const int bad_axis, const double size_ratio=1.0);
	void convolve(cv::Mat_<float> &im, const int axis, const float* kernel);

}

#endif /* DECONVOLUTION_HPP_ */
