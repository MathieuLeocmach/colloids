/*
 * deconvolution.cpp
 *
 *  Created on: 4 juin 2012
 *      Author: mathieu
 */

#define BOOST_TEST_DYN_LINK

#include "../src/deconvolution.hpp"
#include "../src/multiscalefinder.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>
#include <boost/array.hpp>


using namespace Colloids;
using namespace boost::posix_time;

BOOST_AUTO_TEST_SUITE( Deconvolution )
	BOOST_AUTO_TEST_CASE( convolver )
	{
		//tests both even and odd sizes
		boost::array<int,2> sizes ={{11,10}};
		boost::array<float,6> spectrum;
		for (boost::array<int,2>::const_iterator s=sizes.begin(); s!=sizes.end(); ++s)
		{
			Convolver co(*s);
			BOOST_REQUIRE_EQUAL(co.size(), *s);
			BOOST_REQUIRE_EQUAL(co.fourier_size(), 6);
			std::vector<float> input(*s);
			//spectrum from empty input
			std::fill(input.begin(), input.end(), 0.0f);
			co.spectrum(&input[0], 1, &spectrum[0]);
			BOOST_CHECK_EQUAL(*std::max_element(spectrum.begin(), spectrum.end()), 0);
			BOOST_CHECK_EQUAL(*std::min_element(spectrum.begin(), spectrum.end()), 0);
			//step input
			input[0]=1;
			input[1]=1;
			co.spectrum(&input[0], 1, &spectrum[0]);
			BOOST_CHECK_EQUAL(spectrum[0], 4.0);
			//striding input
			std::vector<float> input2(2*(*s));
			std::fill(input2.begin(), input2.end(), 0.0f);
			input2[0]=1;
			input2[1]=1;
			input2[2]=1;
			co.spectrum(&input2[0], 2, &spectrum[0]);
			BOOST_CHECK_EQUAL(spectrum[0], 4.0);
			//simplest convolution: remove of the DC coefficient = subtract the average
			boost::array<float,6> kernel ={{0,1,1,1,1,1}};
			co(&input[0], 1, &kernel[0]);
			BOOST_CHECK_CLOSE(input[0], 1-2.0 / *s, 1e-5);
			BOOST_CHECK_CLOSE(input[2], -2.0 / *s, 1e-5);
		}
		//for a size of 10 the result is exactly null
		BOOST_CHECK_EQUAL(spectrum[5], 0.0);
		fftwf_cleanup();
	}
	BOOST_AUTO_TEST_CASE( spectrum )
	{
		cv::Mat_<float>input(32, 32);
		input.setTo(0.0f);
		cv::circle(input, cv::Point(16, 16), 4, 255, -1);
		std::vector<float> spy = get_spectrum_1d(input, 0);
		std::vector<float> spx = get_spectrum_1d(input, 1);
		BOOST_CHECK_GT(*std::max_element(spx.begin(), spx.end()), 0);
		BOOST_CHECK_CLOSE(spx[0], 668538.28125, 1e-5);
		BOOST_REQUIRE_EQUAL(spy.size(), spx.size());
		for(size_t i=0; i<spy.size(); ++i)
			BOOST_CHECK_MESSAGE(spx[i] == spy[i], "at i="<<i<<"\t"<<spx[i]<<" != "<<spy[i]);
		//isotropic blur
		cv::GaussianBlur(input, input, cv::Size(0,0), 2*1.6);
		spy = get_spectrum_1d(input, 0);
		spx = get_spectrum_1d(input, 1);
		BOOST_CHECK_GT(*std::max_element(spx.begin(), spx.end()), 0);
		BOOST_CHECK_CLOSE(spx[0], 364005.77, 1e-2);
		BOOST_REQUIRE_EQUAL(spy.size(), spx.size());
		for(size_t i=0; i<spy.size(); ++i)
			BOOST_CHECK_CLOSE(spx[i], spy[i], 0.4);

	}
	BOOST_AUTO_TEST_CASE( Gaussian_y )
	{
		cv::Mat_<float>input(32, 32);
		input.setTo(0.0f);
		cv::circle(input, cv::Point(16, 16), 4, 255, -1);
		//anisotropic blur
		cv::GaussianBlur(input, input, cv::Size(0,0), 2*1.6, 1.6);
		//deconvolution kernel
		std::vector<float> kernel = get_deconv_kernel(input, 1, 0);
		BOOST_CHECK_LT(*std::max_element(kernel.begin(), kernel.end()), 100);
		BOOST_CHECK_GT(*std::min_element(kernel.begin(), kernel.end()), 0);
		//deconvolve y
		convolve(input, 0, &kernel[0]);
		BOOST_CHECK_GT(*std::max_element(input.begin(), input.end()), 0);
		BOOST_CHECK_GT(input(16,16), input(15,16));
		BOOST_CHECK_GT(input(16,16), input(16,15));
		BOOST_CHECK_CLOSE(input(16,15), input(15,16), 2);
	}
	BOOST_AUTO_TEST_CASE( rectangular_y )
	{
		cv::Mat_<float>input(32, 30);
		input.setTo(0.0f);
		cv::circle(input, cv::Point(15, 16), 4, 255, -1);
		//anisotropic blur
		cv::GaussianBlur(input, input, cv::Size(0,0), 2*1.6, 1.6);
		//deconvolution kernel
		std::vector<float> kernel = get_deconv_kernel(input, 1, 0);
		BOOST_CHECK_LT(*std::max_element(kernel.begin(), kernel.end()), 100);
		BOOST_CHECK_GT(*std::min_element(kernel.begin(), kernel.end()), 0);
		//deconvolve y
		convolve(input, 0, &kernel[0]);
		BOOST_CHECK_GT(*std::max_element(input.begin(), input.end()), 0);
		BOOST_CHECK_GT(input(16,15), input(15,15));
		BOOST_CHECK_GT(input(16,15), input(16,14));
		BOOST_CHECK_CLOSE(input(16,14), input(15,15), 2);
	}
	BOOST_AUTO_TEST_CASE( rectangular_x )
	{
		cv::Mat_<float>input(30, 32);
		input.setTo(0.0f);
		cv::circle(input, cv::Point(16, 15), 4, 255, -1);
		//anisotropic blur
		cv::GaussianBlur(input, input, cv::Size(0,0), 2*1.6, 1.6);
		//deconvolution kernel
		std::vector<float> kernel = get_deconv_kernel(input, 1, 0);
		BOOST_CHECK_LT(*std::max_element(kernel.begin(), kernel.end()), 100);
		BOOST_CHECK_GT(*std::min_element(kernel.begin(), kernel.end()), 0);
		//deconvolve y
		convolve(input, 0, &kernel[0]);
		BOOST_CHECK_GT(*std::max_element(input.begin(), input.end()), 0);
		BOOST_CHECK_GT(input(15,16), input(14,16));
		BOOST_CHECK_GT(input(15,16), input(15,15));
		BOOST_CHECK_CLOSE(input(15,15), input(14,16), 2);
	}
	BOOST_AUTO_TEST_CASE( sampling )
	{
		cv::Mat_<float>input(32, 32), sinput(16,32);
		input.setTo(0.0f);
		sinput.setTo(0.0f);
		cv::circle(input, cv::Point(16, 16), 4, 255, -1);
		//anisotropic blur
		cv::GaussianBlur(input, input, cv::Size(0,0), 2*1.6, 1.6);
		//remove half of the lines
		for(int j=0; j<16; ++j)
			std::copy(&input(2*j,0), (&input(2*j,0))+32, &sinput(j,0));
		BOOST_REQUIRE_CLOSE(std::accumulate(sinput.begin(), sinput.end(), 0.0), 6247.4463, 1e-2);
		//deconvolution kernel
		std::vector<float> kernel = get_deconv_kernel(sinput, 1, 0, 2.0);
		BOOST_REQUIRE_EQUAL(kernel.size(), 9);
		BOOST_CHECK_LT(*std::max_element(kernel.begin(), kernel.end()), 100);
		BOOST_CHECK_GT(*std::min_element(kernel.begin(), kernel.end()), 0);
		//deconvolve y
		convolve(sinput, 0, &kernel[0]);
		BOOST_CHECK_GT(*std::max_element(sinput.begin(), sinput.end()), 0);
		BOOST_CHECK_GT(sinput(8,16), sinput(7,16));
		BOOST_CHECK_GT(sinput(8,16), sinput(8,15));
		BOOST_CHECK_GT(sinput(8,15), sinput(7,16));
		BOOST_CHECK_LT(sinput(8,13), sinput(7,16));
		BOOST_CHECK_CLOSE(sinput(8,14), sinput(7,16), 10);
	}
	BOOST_AUTO_TEST_CASE( feed_finder )
	{
		MultiscaleFinder3D finder(66, 64, 62);
		int dims[3] = {66,64,62};
		OctaveFinder::Image input(3, dims);
		input.setTo(1);
		std::vector<float> kernel = get_deconv_kernel(input, 0, 2, 1.0);
		BOOST_REQUIRE_EQUAL(kernel.size(), 62/2+1);
		finder.load_deconv_kernel(kernel);
	}
BOOST_AUTO_TEST_SUITE_END()

