#define BOOST_TEST_MODULE multiscale test
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "octavefinder.hpp"

using namespace Colloids;

void images_are_close(const cv::Mat &a, const cv::Mat &b, double precision=1e-9)
{
	cv::Mat_<double> M = cv::abs(a)+cv::abs(b);
	cv::Mat_<double> diff = cv::abs(a-b) / (*std::max_element(M.begin(), M.end()));
	double u;
	for(int i=0; i<a.rows; i++)
		for(int j=0; j<a.cols; j++)
		{
			u = diff(i,j);
			BOOST_CHECK_MESSAGE(u<precision, "at x=" << i <<" y=" << j);
			BOOST_CHECK_SMALL(u, precision);
		}
}


BOOST_AUTO_TEST_SUITE( octave )

BOOST_AUTO_TEST_SUITE( octave_constructors )

	BOOST_AUTO_TEST_CASE( octave_constructor_square )
    {
        OctaveFinder finder;
        BOOST_CHECK_EQUAL(finder.get_width(), 256);
        BOOST_CHECK_EQUAL(finder.get_height(), 256);
        BOOST_CHECK_EQUAL(finder.get_n_layers(), 3);
    }
	BOOST_AUTO_TEST_CASE( octave_constructor_rectangular_flat )
    {
        OctaveFinder finder(128, 512, 1, 1.0);
        BOOST_CHECK_EQUAL(finder.get_width(), 128);
        BOOST_CHECK_EQUAL(finder.get_height(), 512);
        BOOST_CHECK_EQUAL(finder.get_n_layers(), 1);
    }
	BOOST_AUTO_TEST_CASE( octave_constructor_rectangular_ten_layers )
    {
		OctaveFinder finder(128, 512, 10, 2.8);
		BOOST_CHECK_EQUAL(finder.get_width(), 128);
		BOOST_CHECK_EQUAL(finder.get_height(), 512);
		BOOST_CHECK_EQUAL(finder.get_n_layers(), 10);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( iterative_radius )

	BOOST_AUTO_TEST_CASE(iterative_radius_default)
	{
		OctaveFinder finder;
		BOOST_CHECK_EQUAL(finder.get_size(0), 2);
		BOOST_CHECK_EQUAL(finder.get_size(1), 3);
		BOOST_CHECK_EQUAL(finder.get_size(2), 4);
		BOOST_CHECK_EQUAL(finder.get_size(3), 5);
		BOOST_CHECK_EQUAL(finder.get_size(4), 6);
		BOOST_CHECK_EQUAL(finder.get_size(5), 7);
		BOOST_CHECK_CLOSE(finder.get_iterative_radius(4), 3.09001559, 1e-6);
		finder.set_radius_preblur(1.0);
		BOOST_CHECK_EQUAL(finder.get_size(0), 1);
		BOOST_CHECK_EQUAL(finder.get_size(1), 2);
		BOOST_CHECK_EQUAL(finder.get_size(2), 2);
		BOOST_CHECK_EQUAL(finder.get_size(3), 3);
		BOOST_CHECK_EQUAL(finder.get_size(4), 4);
		BOOST_CHECK_EQUAL(finder.get_size(5), 4);
	}
	BOOST_AUTO_TEST_CASE(iterative_radius_flat)
    {
        OctaveFinder finder(128, 512, 1, 1.0);
        BOOST_CHECK_EQUAL(finder.get_size(0), 1);
		BOOST_CHECK_EQUAL(finder.get_size(1), 3);
		BOOST_CHECK_EQUAL(finder.get_size(2), 6);
		BOOST_CHECK_EQUAL(finder.get_size(3), 11);
    }
	BOOST_AUTO_TEST_CASE(iterative_radius_ten_layers)
	{
		OctaveFinder finder(128, 512, 10, 2.8);
		BOOST_CHECK_CLOSE(finder.get_iterative_radius(0), 1.07971992, 1e-6);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( octave_fill )

	BOOST_AUTO_TEST_CASE( fill_test )
	{
		OctaveFinder finder(64, 64);
		cv::Mat_<double>input(64, 64), other(64, 64);
		//the finder should contain a copy of the input data
		input.setTo(1);
		finder.fill(input);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], 64*64, 1e-9);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(input)[0], 1e-9);
		images_are_close(finder.get_layersG(0), input);
		//the internal layersG[0] should not be the same object as the input, but a deep copy
		input.setTo(0);
		input(32,32)=1.0;
		BOOST_CHECK_NE(cv::sum(finder.get_layersG(0))[0], 1.0);
		BOOST_CHECK_NE(cv::sum(finder.get_layersG(0))[0], cv::sum(input)[0]);
		finder.fill(input);
		//Gaussian blur should be normalized
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(finder.get_layersG(1))[0], 1e-9);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(1))[0], cv::sum(finder.get_layersG(2))[0], 1e-9);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(2))[0], cv::sum(finder.get_layersG(3))[0], 1e-9);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(3))[0], cv::sum(finder.get_layersG(4))[0], 1e-9);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(finder.get_layersG(4))[0], 1e-9);
		//The second to last Gaussian layer should be a good approxiation
		//to the input blurred by a two time larger radius than the preblur
		cv::GaussianBlur(input, other, cv::Size(0,0), 2*1.6);
		finder.preblur_and_fill(input);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(3))[0], cv::sum(other)[0], 1e-9);
		BOOST_CHECK_CLOSE(finder.get_layersG(3)(32,32), other(32,32), 1e-2);
		images_are_close(finder.get_layersG(3), other, 1e-4);
	}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
