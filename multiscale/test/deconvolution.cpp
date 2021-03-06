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

void images_are_close(const cv::Mat &a, const cv::Mat &b, float precision=1e-5);

BOOST_AUTO_TEST_SUITE( Deconvolution )
	BOOST_AUTO_TEST_CASE( convolver )
	{
		//tests both even and odd sizes
		boost::array<int,2> sizes ={{11,10}};
		boost::array<float,6> spectrum;
		for (boost::array<int,2>::const_iterator s=sizes.begin(); s!=sizes.end(); ++s)
		{
			Convolver co(*s);
			BOOST_REQUIRE(!co.windowing());
			BOOST_REQUIRE_EQUAL(co.size(), *s);
			BOOST_REQUIRE_EQUAL(co.fourier_size(), 6);
			std::vector<float> input(*s);
			//spectrum from empty input
			std::fill(input.begin(), input.end(), 0.0f);
			co.spectrum(&input[0], 1, &spectrum[0]);
			BOOST_CHECK_EQUAL(*std::max_element(spectrum.begin(), spectrum.end()), 0);
			BOOST_CHECK_EQUAL(*std::min_element(spectrum.begin(), spectrum.end()), 0);
			BOOST_REQUIRE(!co.windowing());
			//step input
			input[0]=1;
			input[1]=1;
			co.spectrum(&input[0], 1, &spectrum[0]);
			BOOST_CHECK_EQUAL(spectrum[0], 4.0);
			BOOST_CHECK_GE(*std::min_element(spectrum.begin(), spectrum.end()), 0);
			//striding input
			std::vector<float> input2(2*(*s));
			std::fill(input2.begin(), input2.end(), 0.0f);
			input2[0]=1;
			input2[1]=1;
			input2[2]=1;
			co.spectrum(&input2[0], 2, &spectrum[0]);
			BOOST_CHECK_GE(*std::min_element(spectrum.begin(), spectrum.end()), 0);
			BOOST_CHECK_EQUAL(spectrum[0], 4.0);
			//simplest convolution: remove of the DC coefficient = subtract the average
			boost::array<float,6> kernel ={{0,1,1,1,1,1}};
			co(&input[0], 1, &kernel[0]);
			BOOST_CHECK_CLOSE(input[0], 1-2.0 / *s, 1e-5);
			BOOST_CHECK_CLOSE(input[2], -2.0 / *s, 1e-5);
			//convolution with a step>1
			co(&input2[0], 2, &kernel[0]);
			BOOST_CHECK_CLOSE(input2[0], 1-2.0 / *s, 1e-5);
			BOOST_CHECK_CLOSE(input2[1], 1, 1e-5);
			BOOST_CHECK_CLOSE(input2[2], 1-2.0 / *s, 1e-5);
			BOOST_CHECK_CLOSE(input2[3], 0, 1e-5);
			BOOST_CHECK_CLOSE(input2[4], -2.0 / *s, 1e-5);
			//windowing ON/OFF
			BOOST_REQUIRE(!co.windowing());
			co.set_hanning();
			BOOST_REQUIRE(co.windowing());
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
		std::vector<float> spy = get_spectrum_1d(input, 0, false);
		std::vector<float> spx = get_spectrum_1d(input, 1, false);
		BOOST_CHECK_GE(*std::min_element(spx.begin(), spx.end()), 0);
		BOOST_CHECK_GT(*std::max_element(spx.begin(), spx.end()), 0);
		BOOST_CHECK_CLOSE(spx[0], 516072.4365234375, 1e-5);
		BOOST_REQUIRE_EQUAL(spy.size(), spx.size());
		for(size_t i=0; i<spy.size(); ++i)
			BOOST_CHECK_MESSAGE(spx[i] == spy[i], "at i="<<i<<"\t"<<spx[i]<<" != "<<spy[i]);
		//isotropic blur
		cv::GaussianBlur(input, input, cv::Size(0,0), 2*1.6);
		spy = get_spectrum_1d(input, 0, false);
		spx = get_spectrum_1d(input, 1, false);
		BOOST_CHECK_GT(*std::max_element(spx.begin(), spx.end()), 0);
		BOOST_CHECK_CLOSE(spx[0], 211539.93515454759, 1e-2);
		BOOST_REQUIRE_EQUAL(spy.size(), spx.size());
		for(size_t i=0; i<spy.size(); ++i)
			BOOST_CHECK_CLOSE(spx[i], spy[i], 0.4);
		//including windowing
		input.setTo(0.0f);
		cv::circle(input, cv::Point(16, 16), 4, 255, -1);
		spy = get_spectrum_1d(input, 0);
		spx = get_spectrum_1d(input, 1);
		BOOST_CHECK_GT(*std::max_element(spx.begin(), spx.end()), 0);
		BOOST_CHECK_CLOSE(spx[0], 469144.84933865839, 1e-5);
		BOOST_REQUIRE_EQUAL(spy.size(), spx.size());
		for(size_t i=0; i<spy.size(); ++i)
			BOOST_CHECK_MESSAGE(spx[i] == spy[i], "at i="<<i<<"\t"<<spx[i]<<" != "<<spy[i]);
	}
	BOOST_AUTO_TEST_CASE( convolution )
	{
		cv::Mat_<float>input(32, 32);
		input.setTo(0.0f);
		cv::circle(input, cv::Point(16, 16), 4, 255, -1);
		cv::Mat_<float> original = input.clone();
		//Convolve by identity kernel
		std::vector<float> kernel(17, 1.0f);
		convolve(input, 0, &kernel[0]);
		images_are_close(input, original, 1e-2);
		//Convolve to double the mean
		kernel.assign(17,2);
		convolve(input, 0, &kernel[0]);
		images_are_close(0.5*input, original, 1e-2);
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
		BOOST_CHECK_CLOSE(input(16,15), input(15,16), 1.5);
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
		BOOST_CHECK_CLOSE(input(16,14), input(15,15), 1.5);
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
		BOOST_CHECK_CLOSE(input(15,15), input(14,16), 1.5);
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
		BOOST_CHECK_CLOSE(sinput(8,14), sinput(7,16), 8);
	}
	BOOST_AUTO_TEST_CASE( feed_finder )
	{
		MultiscaleFinder3D finder(66, 64, 62);
		int dims[3] = {66,64,62};
		OctaveFinder::Image input(3, dims);
		input.setTo(0);
		drawsphere(input, 32, 32, 32, 4.0, (OctaveFinder::PixelType)1.0);
		//anisotropic blur
		inplace_blur3D(input, 0.5, 8);
		std::vector<float> kernel = get_deconv_kernel(input, 2, 0, 1.0);
		BOOST_REQUIRE_EQUAL(kernel.size(), 66/2+1);
		finder.load_deconv_kernel(kernel);
		finder.set_deconv();
		std::vector<Center3D> v;
		finder.get_centers(input, v);
		BOOST_REQUIRE_EQUAL(v.size(), 1);
		BOOST_CHECK_CLOSE(v[0][0], 32, 2);
		BOOST_CHECK_CLOSE(v[0][1], v[0][0], 0.01);
		BOOST_CHECK_CLOSE(v[0][2], v[0][1], 0.01);
	}
	BOOST_AUTO_TEST_CASE( Z_elong )
	{
		int dims[3] = {55, 50, 50};
		cv::Mat_<uchar> image(3, dims);
		image.setTo(0);
		//read test data from disk
		std::ifstream imf("test_input/Z_elong.raw");
		BOOST_REQUIRE_MESSAGE(imf.good(), "could not find test_input/Z_elong.raw");
		imf.read((char*)image.data, 55*50*50);
		imf.close();
		//track in 3D
		MultiscaleFinder3D finder(55, 50, 50, 3, 1.6, true);
		finder.disable_Octave0();
		std::vector<Center3D> centers;
		finder.get_centers(image, centers);
		//check track output
		BOOST_REQUIRE(!centers.empty());
		BOOST_CHECK_EQUAL(centers.size(), 4);
		//two centers should be closer than the sum of their radii
		int tooclose=0;
		for(size_t p=0; p<centers.size()-1; ++p)
			for(size_t q=p+1; q<centers.size(); ++q)
			{
				const double disq = centers[q] - centers[p];
				if(disq < pow(centers[p].r + centers[q].r, 2))
				{
					++tooclose;
					BOOST_WARN_MESSAGE(false, "dist = "<<sqrt(disq)<<"\tR_"<<p<<" = "<<centers[p].r<<"\tR_"<<q<<" = "<<centers[q].r);
				}
			}
		{
			std::ofstream omfc("test_output/Z_elong.xyz");
			for(size_t p=0; p<centers.size(); ++p)
			{
				for(int d=0; d<3; ++d)
					omfc<<centers[p][d]<<"\t";
				omfc<<centers[p].r<<"\n";
			}
		}
		BOOST_CHECK_EQUAL(tooclose, 1);
		//save first octave as reference
		std::ofstream omf("test_output/Z_elong_octaves.raw");
		omf.write(
				(const char*)finder.get_octave(1).get_layersG(0).data,
				55*50*50*6*sizeof(OctaveFinder::PixelType));
		omf.close();

		//Influence of deconvolution on OctaveFinder
		OctaveFinder::Image input;
		image.convertTo(input, input.type());
		std::vector<float> kernel = get_deconv_kernel(input, 2, 0, 1.0);
		std::ofstream kernelf("test_output/Z_elong.kernel");
		std::copy(
			kernel.begin(), kernel.end(),
			std::ostream_iterator<float>(kernelf, "\t"));
		kernelf.close();
		inplace_blur3D(input, 1.6, 1.0);
		convolve(input, 0, &(kernel[0]));
		OctaveFinder3D fiXY(55, 50, 50, 3, 1.6, true);
		fiXY.fill(input);
		fiXY.initialize_binary();
		centers.clear();
		fiXY.subpix(centers);
		BOOST_REQUIRE(!centers.empty());
		BOOST_CHECK_GE(centers.size(), 4);
		for(size_t c=0; c< centers.size(); ++c)
			fiXY.scale(centers[c]);
		//no two centers should be closer than the sum of their radii
		tooclose=0;
		for(size_t p=0; p<centers.size()-1; ++p)
			for(size_t q=p+1; q<centers.size(); ++q)
			{
				const double disq = centers[q] - centers[p];
				if(disq < pow(centers[p].r + centers[q].r, 2))
				{
					++tooclose;
					BOOST_WARN_MESSAGE(false, "dist = "<<sqrt(disq)<<"\tR_"<<p<<" = "<<centers[p].r<<"\tR_"<<q<<" = "<<centers[q].r);
				}
			}
		BOOST_CHECK_EQUAL(tooclose, 0);

		//influence of deconvolution on Multiscale finder
		image.convertTo(input, input.type());
		finder.load_deconv_kernel(kernel);
		finder.set_deconv();
		centers.clear();
		finder.get_centers(input, centers);
		std::ofstream omf2("test_output/Z_elong_deconv_octaves.raw");
		omf2.write(
				(const char*)finder.get_octave(1).get_layersG(0).data,
				55*50*50*6*sizeof(OctaveFinder::PixelType));
		omf2.close();
		//no two centers should be closer than the sum of their radii
		tooclose=0;
		for(size_t p=0; p<centers.size()-1; ++p)
			for(size_t q=p+1; q<centers.size(); ++q)
			{
				const double disq = centers[q] - centers[p];
				if(disq < pow(centers[p].r + centers[q].r, 2))
				{
					++tooclose;
					BOOST_WARN_MESSAGE(false, "dist = "<<sqrt(disq)<<"\tR_"<<p<<" = "<<centers[p].r<<"\tR_"<<q<<" = "<<centers[q].r);
				}
			}
		{
			std::ofstream omfc("test_output/Z_elong_deconv.xyz");
			for(size_t p=0; p<centers.size(); ++p)
			{
				for(int d=0; d<3; ++d)
					omfc<<centers[p][d]<<"\t";
				omfc<<centers[p].r<<"\n";
			}
		}
		BOOST_CHECK_EQUAL(tooclose, 0);
	}
	BOOST_AUTO_TEST_CASE( gel )
	{
		int dims[3] = {31, 64, 64};
		cv::Mat_<uchar> image(3, dims);
		image.setTo(0);
		//read test data from disk
		std::ifstream imf("test_input/gel.raw");
		BOOST_REQUIRE_MESSAGE(imf.good(), "could not find test_input/gel.raw");
		imf.read((char*)image.data, 31*64*64);
		imf.close();
		//deconvolution kernel
		std::vector<float> kernel = get_deconv_kernel(image, 2, 0, 1.0);
		//track in 3D
		MultiscaleFinder3D finder(31, 64, 64);
		finder.load_deconv_kernel(kernel);
		finder.set_deconv();
		std::vector<Center3D> centers;
		finder.get_centers(image, centers);
		//output deconvolved image
		std::ofstream imout("test_output/gel_deconv.raw", std::ios_base::binary);
		imout.write((char*)finder.get_octave(1).get_layersG(0).data, 31*64*64*sizeof(OctaveFinder::PixelType));
		imout.close();
		//check track output
		BOOST_REQUIRE(!centers.empty());
		int contains = 0;
		for(size_t p=0; p<centers.size(); ++p)
		{
			if(centers[p][0]>30 && centers[p][0]<39 &&
					centers[p][1]>15 && centers[p][1]<25 &&
					centers[p][2]>11 && centers[p][2]<25)
				++contains;
		}
		BOOST_REQUIRE_MESSAGE(contains>0, "Bridge particle not detected");
		BOOST_CHECK_MESSAGE(contains<2, ""<<contains<<" bridge particles detected");
		removeOverlapping(centers);
		std::ofstream out("test_output/gel_deconv.csv");
		contains = 0;
		for(size_t p=0; p<centers.size(); ++p)
		{
			if(centers[p][0]>30 && centers[p][0]<39 &&
					centers[p][1]>15 && centers[p][1]<25 &&
					centers[p][2]>11 && centers[p][2]<25)
				++contains;
			for(int d=0;d<3;++d)
				out<<centers[p][d]<<";";
			out<<centers[p].r<<";"<<centers[p].intensity<<"\n";
		}
		BOOST_CHECK_MESSAGE(contains==1, ""<<contains<<" bridge particles detected after overlap removal");
	}
BOOST_AUTO_TEST_SUITE_END()

