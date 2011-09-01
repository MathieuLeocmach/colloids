#define BOOST_TEST_MODULE multiscale test
#define BOOST_TEST_DYN_LINK

#include "multiscalefinder.hpp"
#include "traj.hpp"
#include "reconstructor.h"
#include "locatorfromlif.h"
#include <fstream>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
#include <boost/array.hpp>
#include <set>
#include <limits>

using namespace Colloids;

void images_are_close(const cv::Mat &a, const cv::Mat &b, float precision=1e-5)
{
	OctaveFinder::Image  M = cv::abs(a)+cv::abs(b);
	double peak = *std::max_element(M.begin(), M.end());
	BOOST_REQUIRE_CLOSE(cv::sum(a)[0], cv::sum(b)[0], precision/peak);
	OctaveFinder::Image  diff = cv::abs(a-b) / peak;
	cv::MatConstIterator_<float> u = diff.begin();
	for(int i=0; i<a.rows; i++)
		for(int j=0; j<a.cols; j++)
		{
			BOOST_REQUIRE_MESSAGE(*u<precision, "at x=" << i <<" y=" << j);
			BOOST_CHECK_SMALL(*u, precision);
			u++;
		}
}

template<class T>
struct by_coordinate : std::binary_function<const T&, const T&, bool>
{
	size_t N;
	by_coordinate(size_t n): N(n){};
	bool operator()(const T & a, const T &b)
	{
		return a[N] < b[N];
	}
};

BOOST_AUTO_TEST_SUITE( twoD )

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
		BOOST_CHECK_EQUAL(finder.get_size(2), 3);
		BOOST_CHECK_EQUAL(finder.get_size(3), 4);
		BOOST_CHECK_EQUAL(finder.get_size(4), 5);
		BOOST_CHECK_EQUAL(finder.get_size(5), 6);
		BOOST_CHECK_CLOSE(finder.get_iterative_radius(4), 3.09001559, 1e-6);
		BOOST_CHECK_CLOSE(finder.get_iterative_radius(4), finder.get_iterative_radius(5.0, 4.0), 1e-9);
		finder.set_radius_preblur(1.0);
		BOOST_CHECK_EQUAL(finder.get_size(0), 1);
		BOOST_CHECK_EQUAL(finder.get_size(1), 2);
		BOOST_CHECK_EQUAL(finder.get_size(2), 2);
		BOOST_CHECK_EQUAL(finder.get_size(3), 3);
		BOOST_CHECK_EQUAL(finder.get_size(4), 3);
		BOOST_CHECK_EQUAL(finder.get_size(5), 4);
		BOOST_CHECK_CLOSE(finder.get_iterative_radius(4), finder.get_iterative_radius(5.0, 4.0), 1e-9);
	}
	BOOST_AUTO_TEST_CASE(iterative_radius_flat)
    {
        OctaveFinder finder(128, 512, 1, 1.0);
        BOOST_CHECK_EQUAL(finder.get_size(0), 1);
		BOOST_CHECK_EQUAL(finder.get_size(1), 2);
		BOOST_CHECK_EQUAL(finder.get_size(2), 4);
		BOOST_CHECK_EQUAL(finder.get_size(3), 8);
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
		OctaveFinder::Image input(64, 64), other(64, 64);
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
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(finder.get_layersG(1))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(1))[0], cv::sum(finder.get_layersG(2))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(2))[0], cv::sum(finder.get_layersG(3))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(3))[0], cv::sum(finder.get_layersG(4))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(finder.get_layersG(4))[0], 1e-5);
		//The second to last Gaussian layer should be a good approxiation
		//to the input blurred by a two time larger radius than the preblur
		cv::GaussianBlur(input, other, cv::Size(0,0), 2*1.6);
		finder.preblur_and_fill(input);
		BOOST_REQUIRE_CLOSE(cv::sum(finder.get_layersG(3))[0], cv::sum(other)[0], 1e-5);
		BOOST_CHECK_CLOSE(finder.get_layersG(3)(32,32), other(32,32), 1e-2);
		images_are_close(finder.get_layersG(3), other, 1e-4);
		//Sum of layers should reconstruct blurred image
		other = finder.get_layersG(0) + finder.get_layers(0);
		BOOST_REQUIRE_CLOSE(cv::sum(finder.get_layersG(1))[0], cv::sum(other)[0], 1e-9);
		images_are_close(finder.get_layersG(1), other);
		other += finder.get_layers(1);
		images_are_close(finder.get_layersG(2), other);
		other += finder.get_layers(2);
		images_are_close(finder.get_layersG(3), other);
		other += finder.get_layers(3);
		images_are_close(finder.get_layersG(4), other);
		//the inner data should not change type when filled with something else than double
		cv::Mat_<unsigned char>input2(64, 64);
		input2.setTo(0);
		input2(32,32)=1;
		finder.fill(input2);
		BOOST_CHECK_EQUAL(finder.get_layersG(0).type(), input.type());
	}
	BOOST_AUTO_TEST_CASE( fill_rectangular )
	{
		OctaveFinder finder(51, 103);
		OctaveFinder::Image input(51, 103), other(51, 103);
		//the finder should contain a copy of the input data
		input.setTo(0);
		input(32,32)=1.0;
		finder.fill(input);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(finder.get_layersG(1))[0], 1e-4);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(1))[0], cv::sum(finder.get_layersG(2))[0], 1e-4);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(2))[0], cv::sum(finder.get_layersG(3))[0], 1e-4);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(3))[0], cv::sum(finder.get_layersG(4))[0], 1e-4);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(finder.get_layersG(4))[0], 1e-4);
		cv::GaussianBlur(input, other, cv::Size(0,0), 2*1.6);
		finder.preblur_and_fill(input);
		BOOST_REQUIRE_CLOSE(cv::sum(finder.get_layersG(3))[0], cv::sum(other)[0], 1e-4);
		BOOST_CHECK_CLOSE(finder.get_layersG(3)(32,32), other(32,32), 1e-2);
		images_are_close(finder.get_layersG(3), other, 1e-4);
	}
	BOOST_AUTO_TEST_CASE( fill_speed )
	{
		OctaveFinder finder;
		OctaveFinder::Image input(256, 256);
		input.setTo(0);
		cv::circle(input, cv::Point(128, 128), 4, 1.0, -1);
		boost::progress_timer ti;
		for (size_t i=0; i<100; ++i)
			finder.fill(input);
		std::cout<<"100 fill in ";
	}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( local_max )

	BOOST_AUTO_TEST_CASE( single_circle )
	{
		OctaveFinder finder;
		OctaveFinder::Image input(256, 256);
		input.setTo(0);
		cv::circle(input, cv::Point(128, 128), 4, 1.0, -1);
		finder.preblur_and_fill(input);
		//with a radius of 4, the maximum should be in layer 2
		BOOST_CHECK_GE(finder.get_layers(0)(128, 128), finder.get_layers(1)(128, 128));
		BOOST_CHECK_GE(finder.get_layers(1)(128, 128), finder.get_layers(2)(128, 128));
		BOOST_CHECK_LE(finder.get_layers(2)(128, 128), finder.get_layers(3)(128, 128));
		BOOST_CHECK_LE(finder.get_layers(3)(128, 128), finder.get_layers(4)(128, 128));
		finder.initialize_binary();

		//The minimum should be at the center of the circle
		double min_layers[5];
		for (int l=0; l<5; ++l)
			min_layers[l] =  *std::min_element(finder.get_layers(l).begin(), finder.get_layers(l).end());
		double global_min = *std::min_element(min_layers, min_layers+5);
		BOOST_CHECK_CLOSE(finder.get_layers(2)(128, 128), global_min, 1e-9);

		BOOST_CHECK_EQUAL(finder.get_binary(2)(128, 128), true);
		BOOST_CHECK_EQUAL(finder.get_binary(2)(129, 128), false);
		BOOST_CHECK_EQUAL(finder.get_binary(2)(128, 129), false);
		BOOST_CHECK_EQUAL(finder.get_binary(2)(127, 128), false);
		BOOST_CHECK_EQUAL(finder.get_binary(2)(128, 127), false);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(2))[0], 1);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(1))[0], 0);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(3))[0], 0);
		/*cv::namedWindow("truc");
		cv::imshow("truc", 255*finder.get_binary(2));
		cv::waitKey();*/
		//gaussian response
		BOOST_CHECK_LT(finder.gaussianResponse(128, 128, 2.5), finder.gaussianResponse(128, 128, 2.0));
		BOOST_CHECK_LT(finder.gaussianResponse(128, 128, 3.5), finder.gaussianResponse(128, 128, 2.5));
		BOOST_CHECK_LT(finder.gaussianResponse(128, 128, 3.5), finder.gaussianResponse(128, 128, 3.0));
		BOOST_CHECK_GT(finder.gaussianResponse(128, 128, 2.5)-finder.gaussianResponse(128, 128, 1.5), finder.get_layers(2)(128, 128));
		BOOST_CHECK_GT(finder.gaussianResponse(128, 128, 3.5)-finder.gaussianResponse(128, 128, 2.5), finder.get_layers(2)(128, 128));
		//lower bound
		BOOST_CHECK_LT(finder.gaussianResponse(128, 128, 1.5), finder.gaussianResponse(128, 128, 1.0));
		BOOST_CHECK_LT(finder.gaussianResponse(128, 128, 0.5), finder.gaussianResponse(128, 128, 0));
		//further than the top layer
		BOOST_CHECK_LT(finder.gaussianResponse(128, 128, 10.5), finder.gaussianResponse(128, 128, 0));
		BOOST_CHECK_LT(finder.gaussianResponse(128, 128, 15.5), finder.gaussianResponse(128, 128, 5));
	}

	BOOST_AUTO_TEST_CASE( empty_rectangular )
	{
		OctaveFinder finder(128, 512, 10, 2.8);
		OctaveFinder::Image input(128, 512);
		input.setTo(0);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		for(int i=1; i<10; ++i)
			BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(i))[0], 0);
		input.setTo(-1);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		for(int i=1; i<10; ++i)
			BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(i))[0], 0);
	}
	BOOST_AUTO_TEST_CASE( local_max_rectangular )
	{
		OctaveFinder finder(151, 103);
		OctaveFinder::Image input(151, 103), other(151, 103);
		//the finder should contain a copy of the input data
		input.setTo(0);
		cv::circle(input, cv::Point(32, 32), 4, 1.0, -1);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		BOOST_CHECK_EQUAL(finder.get_binary(2)(32, 32), true);
		BOOST_CHECK_EQUAL(finder.get_binary(2)(33, 32), false);
		BOOST_CHECK_EQUAL(finder.get_binary(2)(32, 33), false);
		BOOST_CHECK_EQUAL(finder.get_binary(2)(31, 32), false);
		BOOST_CHECK_EQUAL(finder.get_binary(2)(32, 31), false);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(2))[0], 1);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(1))[0], 0);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(3))[0], 0);
		/*cv::imshow("layer",-finder.get_layers(2));
		cv::imshow("binary", 255*finder.get_binary(2));
		cv::waitKey();*/
	}

	BOOST_AUTO_TEST_CASE( multiple_circles_on_edges )
	{
		OctaveFinder finder;
		OctaveFinder::Image input(256, 256);
		input.setTo(0);
		//circles on the edges
		for(int i=0; i<8; ++i)
		{
			cv::circle(input, cv::Point(i-2, (i+1)*32), 4, 1.0, -1);
			cv::circle(input, cv::Point((i+1)*32, 257-i), 4, 1.0, -1);
		}
		for(int i=0; i<7; ++i)
			BOOST_CHECK_CLOSE(input((i+1)*32, 0), 1.0, 1e-9);
		BOOST_CHECK_CLOSE(input(0, 7*32), 0, 1e-9);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		//only the circles further than 1 px from the border can be detected
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(2))[0], 2);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(1))[0], 0);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(3))[0], 0);
		BOOST_CHECK_EQUAL(finder.get_binary(2)(7*32, 7-3), true);
		BOOST_CHECK_EQUAL(finder.get_binary(2)(258-7, 7*32), true);
		std::vector<Center2D> v;
		finder.subpix(v);
	}

	BOOST_AUTO_TEST_CASE( multiple_circles_various_sizes )
	{
		OctaveFinder finder;
		OctaveFinder::Image input(256, 256);
		input.setTo(0);
		for(int i=0; i<7; ++i)
			cv::circle(input, cv::Point(128, (i+1)*32), i+1, 1.0, -1);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(1))[0], 1);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(2))[0], 1);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(3))[0], 1);
		BOOST_CHECK_LE(finder.get_layers(3)(5*32, 128), 0);
		BOOST_CHECK(1+pow(finder.get_layers(3)(5*32, 128), 2)>1);
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(3)(5*32, 127));
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(3)(5*32, 129));
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(3)(5*32+1, 128));
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(3)(5*32-1, 128));
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(2)(5*32, 128));
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(2)(5*32, 127));
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(2)(5*32, 129));
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(2)(5*32+1, 128));
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(2)(5*32-1, 128));
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(4)(5*32, 128));
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(4)(5*32, 127));
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(4)(5*32, 129));
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(4)(5*32+1, 128));
		BOOST_CHECK_NE(finder.get_layers(3)(5*32, 128), finder.get_layers(4)(5*32-1, 128));
	}

	BOOST_AUTO_TEST_CASE( ellipses )
	{
		OctaveFinder finder;
		OctaveFinder::Image input(256, 256);
		input.setTo(0);
		for(int i=0; i<7; ++i)
			for(int j=0; j<7; ++j)
			cv::ellipse(input, cv::Point((j+1)*32, (i+1)*32), cv::Size(4.0/(1+0.1*i),4*(1+0.1*i)), 90.0/8 *j, 90.0/8 *j, 360+90.0/8 *j, 1.0, -1);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(1))[0], 0);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(2))[0], 21);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(3))[0], 0);
		/*cv::namedWindow("truc");
		cv::imshow("truc", input);
		cv::imshow("binary2", 255*finder.get_binary(2));
		cv::imshow("binary3", 255*finder.get_binary(3));
		cv::waitKey();*/
	}

	BOOST_AUTO_TEST_CASE( local_max_speed )
	{
		OctaveFinder finder;
		OctaveFinder::Image input(256, 256);
		input.setTo(0);
		cv::circle(input, cv::Point(128, 128), 4, 1.0, -1);
		finder.preblur_and_fill(input);
		boost::progress_timer ti;
		for (size_t i=0; i<100; ++i)
			finder.initialize_binary();
		std::cout<<"100 local minima detections in ";
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( subpix )

	BOOST_AUTO_TEST_CASE( subpix_single_circle )
	{
		OctaveFinder finder;
		OctaveFinder::Image input(256, 256);
		input.setTo(0);
		cv::circle(input, cv::Point(100, 200), 4, 1.0, -1);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		BOOST_CHECK_EQUAL(finder.get_binary(2)(200, 100), 1);
		std::vector<Center2D> v;
		finder.subpix(v);
		BOOST_REQUIRE_EQUAL(v.size(), 1);
		//x
		BOOST_CHECK_GE(v[0][0], 99);
		BOOST_CHECK_LE(v[0][0], 101);
		//y
		BOOST_CHECK_GE(v[0][1], 199);
		BOOST_CHECK_LE(v[0][1], 201);
		//scale
		BOOST_CHECK_GE(v[0].r, 1);
		BOOST_CHECK_LE(v[0].r, 3);
		//gaussian response
		BOOST_CHECK_CLOSE(finder.gaussianResponse(200, 100, 2.0), finder.get_layersG(2)(200, 100), 1e-9);
		BOOST_CHECK_CLOSE(finder.gaussianResponse(200, 100, 3.0)-finder.gaussianResponse(200, 100, 2.0), finder.get_layers(2)(200, 100), 1e-9);
		BOOST_CHECK_LT(finder.gaussianResponse(100, 200, 2.0), finder.gaussianResponse(200, 100, 2.0));
		BOOST_CHECK_LT(finder.gaussianResponse(100, 200, 2.5), finder.gaussianResponse(200, 100, 2.5));
		BOOST_CHECK_LT(finder.gaussianResponse(200, 100, 3.0), finder.gaussianResponse(200, 100, 2.0));
		BOOST_CHECK_LT(finder.gaussianResponse(200, 100, 2.5), finder.gaussianResponse(200, 100, 2.0));
		BOOST_CHECK_LT(finder.gaussianResponse(200, 100, 3.5), finder.gaussianResponse(200, 100, 2.5));
		BOOST_CHECK_LT(finder.gaussianResponse(200, 100, 3.5), finder.gaussianResponse(200, 100, 3.0));
		BOOST_CHECK_GT(finder.gaussianResponse(200, 100, 2.5)-finder.gaussianResponse(200, 100, 1.5), finder.get_layers(2)(200, 100));
		BOOST_CHECK_GT(finder.gaussianResponse(200, 100, 3.5)-finder.gaussianResponse(200, 100, 2.5), finder.get_layers(2)(200, 100));
		boost::array<double,8> sublayerG;
		for(int u = 0;u < 3;++u)
			sublayerG[u] = finder.gaussianResponse(200, 100, 2 - 0.5 + u);
		boost::array<double,5> a = {{
			finder.get_layers(1)(200, 100),
			sublayerG[1] - sublayerG[0],
			finder.get_layers(2)(200, 100),
			sublayerG[2] - sublayerG[1],
			finder.get_layers(3)(200, 100)
		}};
		BOOST_CHECK_LT(a[0], 0);
		BOOST_CHECK_LT(a[2], 0);
		BOOST_CHECK_LT(a[4], 0);
		BOOST_REQUIRE_NE(a[4]-2*a[2]+a[0], 0);
		double ds = (-a[4] + 8*a[3] - 8*a[1] + a[0])/6.0 /(a[4]-2*a[2]+a[0]), z = 2-ds;
		BOOST_CHECK_GT(ds, 0);
		BOOST_CHECK_LT(ds, 0.5);
		for(size_t u = 0;u < sublayerG.size(); ++u)
			sublayerG[u] = finder.gaussianResponse(200, 100, z - 1 + 0.5*u);
		for(size_t u =0; u<a.size();++u)
			a[u] = sublayerG[u+2] - sublayerG[u];
		BOOST_CHECK_LT(a[0], 0);
		BOOST_CHECK_LT(a[2], 0);
		BOOST_CHECK_LT(a[4], 0);
		BOOST_REQUIRE_NE(a[4]-2*a[2]+a[0], 0);
		ds = (-a[4] + 8*a[3] - 8*a[1] + a[0])/6.0 /(a[4]-2*a[2]+a[0]);
		BOOST_CHECK_GT(z-ds, 1.5);
		BOOST_CHECK_LT(z-ds, 2.5);

	}

	BOOST_AUTO_TEST_CASE( subpix_rectangular )
	{
		OctaveFinder finder(503, 151);
		OctaveFinder::Image input(503, 151);
		input.setTo(0);
		cv::circle(input, cv::Point(100, 200), 4, 1.0, -1);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		BOOST_CHECK_EQUAL(finder.get_binary(2)(200, 100), 1);
		std::vector<Center2D> v;
		finder.subpix(v);
		BOOST_REQUIRE_EQUAL(v.size(), 1);
		//x
		BOOST_CHECK_GE(v[0][0], 99);
		BOOST_CHECK_LE(v[0][0], 101);
		//y
		BOOST_CHECK_GE(v[0][1], 199);
		BOOST_CHECK_LE(v[0][1], 201);
		//scale
		BOOST_CHECK_GE(v[0].r, 1);
		BOOST_CHECK_LE(v[0].r, 3);

	}

	BOOST_AUTO_TEST_CASE( subpix_relative_positions )
	{
		OctaveFinder finder(32, 32);
		//cv cannot draw circle positions better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(32*32, 32*32), small_input(32,32);
		std::vector<Center2D> v;
		for(int i=0; i<6; ++i)
		{
			input.setTo(0);
			cv::circle(input, cv::Point(32*(16+pow(0.5, i+1)), 32*16), 32*4, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			finder.preblur_and_fill(small_input);
			finder.initialize_binary();
			std::vector<Center2D> v_s;
			finder.subpix(v_s);
			std::copy(v_s.begin(), v_s.end(), std::back_inserter(v));
		}
		BOOST_REQUIRE_EQUAL(v.size(), 6);
		for(int i=0; i<5; ++i)
			BOOST_WARN_MESSAGE(
					v[i][0] > v[i+1][0],
					"spatial resolution is larger than 1/" << pow(2,i+1) <<"th of a pixel"
					);
		BOOST_CHECK_GT(v[2][0], v[3][0]);

	}

	BOOST_AUTO_TEST_CASE( subpix_relative_sizes )
	{
		OctaveFinder finder(32,32);
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(8*32, 8*32), small_input(32,32);
		std::vector<Center2D> v;
		for(int i=0; i<7; ++i)
		{
			input.setTo(0);
			cv::circle(input, cv::Point(8*16, 8*16), 8*(4+0.125*i), 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			finder.preblur_and_fill(small_input);
			finder.initialize_binary();
			std::vector<Center2D> v_s;
			finder.subpix(v_s);
			std::copy(v_s.begin(), v_s.end(), std::back_inserter(v));
		}
		BOOST_REQUIRE_EQUAL(v.size(), 7);

		for (int i=1; i<7;++i)
			BOOST_WARN_MESSAGE(v[0].r< v[7-i].r, "resolution in size is 1/"<<(8*(7-i))<< "th of a scale");

		for(size_t c=0; c<v.size(); ++c)
			finder.scale(v[c]);
		BOOST_CHECK_CLOSE(v[0].r, 4, 2);
		BOOST_CHECK_CLOSE(v[1].r, 4.125, 2);
		BOOST_CHECK_CLOSE(v[2].r, 4.25, 2);
		BOOST_CHECK_CLOSE(v[3].r, 4.325, 2);
		BOOST_CHECK_CLOSE(v[4].r, 4.5, 2);
		BOOST_CHECK_CLOSE(v[5].r, 4.625, 2);
		BOOST_CHECK_CLOSE(v[6].r, 4.75, 2);
		/*for (int i=0; i<30;++i)
		{
			input.setTo(0);
			cv::circle(input, cv::Point(8*16, 8*16), 8*(2+0.125*i), 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			finder.preblur_and_fill(small_input);
			finder.initialize_binary();
			std::vector<cv::Vec4d> v_s = finder.subpix();
			finder.scale(v_s);
			for(size_t j=0; j<v_s.size(); ++j)
				std::cout <<"["<< (2+0.125*i) << ", " << v_s[j][2] << "], ";
		}
		std::cout<<std::endl;*/
	}

	BOOST_AUTO_TEST_CASE( subpix_speed )
	{
		OctaveFinder finder;
		OctaveFinder::Image input(256, 256);
		input.setTo(0);
		cv::circle(input, cv::Point(128, 128), 4, 1.0, -1);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		boost::progress_timer ti;
		std::vector<Center2D> v;
		for (size_t i=0; i<10000; ++i)
			 finder.subpix(v);
		std::cout<<"10000 subpixel resolutions in ";
	}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( octave_limit_cases )

	BOOST_AUTO_TEST_CASE( octave_minimum_detected_size )
	{
		OctaveFinder finder(32,32);
		cv::Mat_<uchar>input(16*32, 16*32), small_input(32,32);
		int i = 0;
		std::vector<Center2D> v(1);
		while(3-0.01*i>0 && v.size()>0)
		{
			i++;
			input.setTo(0);
			cv::circle(input, cv::Point(16*16, 16*16), 16*(3-0.01*i), 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			v =	finder(small_input, true);
		}
		BOOST_WARN_MESSAGE(false, "Cannot detect sizes below "<< (3-0.01*(i-1))<<" pixels");
	}

	BOOST_AUTO_TEST_CASE( octave_minimum_detector_size )
	{
		int i = 16;
		std::vector<Center2D> v(1);
		while(i>0 && v.size()>0)
		{
			i--;
			OctaveFinder finder(i,i);
			cv::Mat_<uchar>input(16*i, 16*i), small_input(i,i);
			input.setTo(0);
			cv::circle(input, cv::Point(8*i, 8*i), 16*2.88, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			v =	finder(small_input, true);
		}
		BOOST_WARN_MESSAGE(false, "An octave detector smaller than "<< (i+1)<<" pixels cannot detect anything");
	}
	BOOST_AUTO_TEST_CASE( size_at_border )
	{
		OctaveFinder finder(32,32);
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(32*32, 32*32), small_input(32,32);
		std::ofstream f("pos_size_at_border.out");
		for(int i = 32*16; i>32*4; --i)
		{
			const double position = i/32.0;
			BOOST_TEST_CHECKPOINT("position = "<<position);
			input.setTo(0);
			cv::circle(input, cv::Point(32*16, i), 32*4, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			std::vector<Center2D> v_s = finder(small_input, true);
			BOOST_REQUIRE_MESSAGE(
				v_s.size()==1,
				""<<((v_s.size()==0)?"No center":"More than one center")<<" for input position "<<position
				);
			BOOST_CHECK_GT(finder.gaussianResponse(position, 16, 0), 0);
			BOOST_CHECK_LT(finder.gaussianResponse(position, 16, 3.5), finder.gaussianResponse(position, 16, 3));
			BOOST_CHECK_LT(finder.gaussianResponse(position, 16, 3), finder.gaussianResponse(position, 16, 2.5));
			BOOST_CHECK_LT(finder.gaussianResponse(position, 16, 2.5), finder.gaussianResponse(position, 16, 2));
			BOOST_CHECK_LT(finder.gaussianResponse(position, 16, 2), finder.gaussianResponse(position, 16, 1.5));
			BOOST_CHECK_LT(finder.gaussianResponse(position, 16, 1.5), finder.gaussianResponse(position, 16, 1));
			BOOST_CHECK_LT(finder.gaussianResponse(position, 16, 1), finder.gaussianResponse(position, 16, 0.5));
			BOOST_CHECK_LT(finder.gaussianResponse(position, 16, 0.5), finder.gaussianResponse(position, 16, 0));
			BOOST_CHECK_CLOSE(v_s[0].r, 4, 50);
			for(size_t j=0; j<v_s.size(); ++j)
				f << position << "\t" << v_s[j][1] << "\t"  << v_s[j].r << "\n";
		}
		f<<std::endl;
	}

	BOOST_AUTO_TEST_CASE( close_neighbours )
	{
		OctaveFinder finder(64,32);
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(32*64, 32*32), small_input(64,32);
		std::ofstream f("close_neighbours.out");
		for(int i = 32*32; i>7*32; --i)
		{
			const double distance = i/32.0;
			BOOST_TEST_CHECKPOINT("distance = "<<distance);
			input.setTo(0);
			cv::circle(input, cv::Point(32*16, 32*16), 32*4, 255, -1);
			cv::circle(input, cv::Point(32*16, 32*16+i), 32*4, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			std::vector<Center2D> v_s = finder(small_input, true);
			BOOST_REQUIRE_EQUAL(v_s.size(), 2);
			BOOST_CHECK_CLOSE(v_s[0][1], 16, distance<16?6:2);
			/*cv::namedWindow("truc");
			cv::imshow("truc", small_input);
			cv::imshow("binary", 255*finder.get_binary(2));
			cv::waitKey();*/
			BOOST_REQUIRE_CLOSE(v_s[1][1], 16+distance, distance<16?6:2);
			BOOST_CHECK_CLOSE(v_s[0].r, 4, distance<16?6:2);
			BOOST_CHECK_CLOSE(v_s[1][1] - v_s[0][1], distance, distance<9.06375?22:2);
			f << distance << "\t" << v_s[0][1] << "\t" << v_s[1][1] << "\t" << v_s[0].r << "\n";
		}
		f<<std::endl;
	}

	BOOST_AUTO_TEST_CASE( octave_minimum_signal_intensity )
	{
		OctaveFinder finder(32,32);
		OctaveFinder::Image  input(32,32);
		int i = 0;
		std::vector<Center2D> v(1);
		while(i<10 && v.size()>0)
		{
			i++;
			input.setTo(0);
			cv::circle(input, cv::Point(16, 16), 3, pow(0.5, i), -1);
			v =	finder(input, true);
		}
		BOOST_WARN_MESSAGE(i>10, "Cannot detect intensities below "<< pow(0.5, i-1));
		BOOST_CHECK_LT(finder.get_layers(1)(16,16), 0);
	}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END() //octave

BOOST_AUTO_TEST_SUITE( multiscale )

BOOST_AUTO_TEST_SUITE( multiscale_constructors )

	BOOST_AUTO_TEST_CASE( multiscale_constructor_square )
	{
		MultiscaleFinder2D finder;
		BOOST_REQUIRE_EQUAL(finder.get_n_octaves(), 6);
		BOOST_CHECK_EQUAL(finder.get_octave(1).get_width(), 256);
		BOOST_CHECK_EQUAL(finder.get_octave(1).get_height(), 256);
		BOOST_CHECK_EQUAL(finder.get_octave(1).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(finder.get_octave(0).get_width(), 512);
		BOOST_CHECK_EQUAL(finder.get_octave(0).get_height(), 512);
		BOOST_CHECK_EQUAL(finder.get_octave(0).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(finder.get_octave(2).get_width(), 128);
		BOOST_CHECK_EQUAL(finder.get_octave(2).get_height(), 128);
		BOOST_CHECK_EQUAL(finder.get_octave(2).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(finder.get_octave(3).get_width(), 64);
		BOOST_CHECK_EQUAL(finder.get_octave(3).get_height(), 64);
		BOOST_CHECK_EQUAL(finder.get_octave(3).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(finder.get_octave(4).get_width(), 32);
		BOOST_CHECK_EQUAL(finder.get_octave(4).get_height(), 32);
		BOOST_CHECK_EQUAL(finder.get_octave(4).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(finder.get_octave(5).get_width(), 16);
		BOOST_CHECK_EQUAL(finder.get_octave(5).get_height(), 16);
		BOOST_CHECK_EQUAL(finder.get_octave(5).get_n_layers(), 3);
	}
	BOOST_AUTO_TEST_CASE( multiscale_constructor_rectangular_flat )
    {
		MultiscaleFinder2D finder(101, 155, 1, 1.0);
        BOOST_REQUIRE_EQUAL(finder.get_n_octaves(), 5);
        BOOST_CHECK_EQUAL(finder.get_octave(1).get_width(), 101);
        BOOST_CHECK_EQUAL(finder.get_octave(1).get_height(), 155);
        BOOST_CHECK_EQUAL(finder.get_octave(0).get_width(), 202);
        BOOST_CHECK_EQUAL(finder.get_octave(0).get_height(), 310);
        BOOST_CHECK_EQUAL(finder.get_octave(2).get_width(), 50);
        BOOST_CHECK_EQUAL(finder.get_octave(2).get_height(), 77);
        BOOST_CHECK_EQUAL(finder.get_octave(3).get_width(), 25);
        BOOST_CHECK_EQUAL(finder.get_octave(3).get_height(), 38);
        BOOST_CHECK_EQUAL(finder.get_octave(4).get_width(), 12);
        BOOST_CHECK_EQUAL(finder.get_octave(4).get_height(), 19);
        BOOST_CHECK_CLOSE(finder.get_radius_preblur(), 1.0, 1e-9);
        finder.set_radius_preblur(1.6);
        BOOST_CHECK_CLOSE(finder.get_radius_preblur(), 1.6, 1e-9);
    }
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( multiscale_call )

	BOOST_AUTO_TEST_CASE( single_circle )
	{
		MultiscaleFinder2D finder;
		OctaveFinder::Image input(256, 256);
		input.setTo(0);
		cv::circle(input, cv::Point(128, 128), 4, 1.0, -1);
		std::vector<Center2D> v = finder(input);
		//the first octave's should be the same as the one of a Octavefinder of the same size
		OctaveFinder ofinder;
		ofinder.preblur_and_fill(input);
		BOOST_CHECK_CLOSE(finder.get_octave(1).get_layers(0)(128,128), ofinder.get_layers(0)(128,128), 1e-9);
		images_are_close(finder.get_octave(1).get_layers(0), ofinder.get_layers(0), 1e-4);
		//with a radius of 4, the maximum should be in layer 2 of octave 1 and nowhere else
		for(size_t o=0; o<finder.get_n_octaves(); ++o)
			for(size_t l=1; l<finder.get_n_layers()+1; ++l)
			{
				const int u = cv::sum(finder.get_octave(o).get_binary(l))[0];
				BOOST_CHECK_MESSAGE(
						u == ((o==1 && l==2)?1:0),
						"Octave "<<o<<" layer "<<l<<" has "<< u <<" center"<<((u>1)?"s":"")
				);
			}
		BOOST_REQUIRE_EQUAL(v.size(), 1);
		BOOST_CHECK_CLOSE(v[0][0], 128, 10);
		BOOST_CHECK_CLOSE(v[0][1], 128, 10);
		BOOST_CHECK_CLOSE(v[0].r, 4, 2);
	}
	BOOST_AUTO_TEST_CASE( multiscale_rectangular )
	{
		MultiscaleFinder2D finder(101, 155);
		BOOST_REQUIRE_EQUAL(finder.get_n_octaves(), 5);
		OctaveFinder::Image input(101, 155);
		input.setTo(0);
		cv::circle(input, cv::Point(28, 28), 2, 1.0, -1);
		std::vector<Center2D> v = finder(input);
		//the first octave's should be the same as the one of a OctaveFinder of the same size
		OctaveFinder ofinder(101, 155);
		ofinder.preblur_and_fill(input);
		BOOST_CHECK_CLOSE(finder.get_octave(1).get_layers(0)(28,28), ofinder.get_layers(0)(28,28), 1e-9);
		images_are_close(finder.get_octave(1).get_layers(0), ofinder.get_layers(0), 1e-4);

		//with a radius of 2, the maximum should be in layer 2 of octave 0 and nowhere else
		BOOST_CHECK_LT(finder.get_octave(0).get_layers(2)(56, 56), finder.get_octave(0).get_layers(2)(55, 56));
		BOOST_CHECK_LT(finder.get_octave(0).get_layers(2)(56, 56), finder.get_octave(0).get_layers(2)(57, 56));
		BOOST_CHECK_LT(finder.get_octave(0).get_layers(2)(56, 56), finder.get_octave(0).get_layers(2)(56, 55));
		BOOST_CHECK_LT(finder.get_octave(0).get_layers(2)(56, 56), finder.get_octave(0).get_layers(2)(56, 57));
		//cv::imshow("o0 binary2", 255*finder.get_octave(0).get_binary(2));
		//cv::waitKey();
		for(size_t o=0; o<finder.get_n_octaves(); ++o)
			for(size_t l=1; l<finder.get_n_layers()+1; ++l)
			{
				const int u = cv::sum(finder.get_octave(o).get_binary(l))[0];
				BOOST_CHECK_MESSAGE(
						u == ((o==0 && l==2)?1:0),
						"Octave "<<o<<", layer "<<l<<" has "<< u <<" center"<<((u>1)?"s":"")
				);
			}
		BOOST_REQUIRE_EQUAL(v.size(), 1);
		BOOST_CHECK_CLOSE(v[0][0], 28, 10);
		BOOST_CHECK_CLOSE(v[0][1], 28, 10);
		BOOST_CHECK_CLOSE(v[0].r, 2, 5);
	}
	BOOST_AUTO_TEST_CASE( multiscale_relative_sizes )
	{
		const int s = 5;
		MultiscaleFinder2D finder(pow(2, s+1), pow(2, s+1));
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(32*pow(2, s+1), 32*pow(2, s+1)), small_input(pow(2, s+1), pow(2, s+1));
		std::vector<Center2D> v;
		std::ofstream f("multiscale_relative_sizes.out");
		std::set<int> large_radii;
		for(int k=1; k<s-1; ++k)
			for(int i=0; i<32; ++i)
				large_radii.insert(48*pow(2, k-1)+i*pow(2, k));
		for(std::set<int>::const_iterator lr = large_radii.begin(); lr != large_radii.end(); ++lr)
		{
			input.setTo(0);
			const double radius =  *lr/32.0;
			BOOST_TEST_CHECKPOINT("radius = "<<radius<<" call");
			cv::circle(input, cv::Point(32*pow(2, s), 32*pow(2, s)), *lr, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			std::vector<Center2D> v_s = finder(small_input);
			for(size_t j=0; j<v_s.size(); ++j)
				f << radius << "\t" << v_s[j].r << "\n";
			BOOST_TEST_CHECKPOINT("radius = "<<radius<<" size");
			BOOST_REQUIRE_MESSAGE(
					v_s.size()==1,
					""<<((v_s.size()==0)?"No center":"More than one center")<<" for input radius "<<radius
					);
			BOOST_CHECK_CLOSE(v_s[0][1], pow(2, s), 2);
			/*if(radius<pow(2, s-1))
				BOOST_CHECK_CLOSE(v_s[0][2], radius, 5);*/
			BOOST_TEST_CHECKPOINT("radius = "<<radius<<" copy");
			std::copy(v_s.begin(), v_s.end(), std::back_inserter(v));

		}
		f<<std::endl;
		BOOST_REQUIRE_EQUAL(v.size(), large_radii.size());
	}
	BOOST_AUTO_TEST_CASE( multiscale_minimum_detector_size )
	{
		int i = 16;
		std::vector<Center2D> v(1);
		while(i>6 && v.size()>0)
		{
			i--;
			BOOST_CHECKPOINT("i="<<i);
			MultiscaleFinder2D finder(i,i);
			cv::Mat_<uchar>input(16*i, 16*i), small_input(i,i);
			input.setTo(0);
			cv::circle(input, cv::Point(8*i, 8*i), 16*2.88, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			v =	finder(small_input);
		}
		BOOST_WARN_MESSAGE(false, "An multiscale detector smaller than "<< (i+1)<<" pixels cannot detect anything");
	}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END() //multiscale

BOOST_AUTO_TEST_SUITE_END() //2D

BOOST_AUTO_TEST_SUITE( oneD )

BOOST_AUTO_TEST_SUITE( octave )

	BOOST_AUTO_TEST_CASE( constructor )
	{
		OctaveFinder1D finder(128);
		BOOST_CHECK_EQUAL(finder.get_width(), 1);
		BOOST_CHECK_EQUAL(finder.get_height(), 128);
		BOOST_CHECK_EQUAL(finder.get_n_layers(), 3);
	}
	BOOST_AUTO_TEST_CASE( fill )
	{
		OctaveFinder1D finder(64);
		OctaveFinder::Image input(1, 64), other(1, 64);
		//the finder should contain a copy of the input data
		input.setTo(1);
		finder.fill(input);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], 64, 1e-9);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(input)[0], 1e-9);
		images_are_close(finder.get_layersG(0), input);
		//the internal layersG[0] should not be the same object as the input, but a deep copy
		input.setTo(0);
		input(0, 32)=1.0;
		BOOST_CHECK_NE(cv::sum(finder.get_layersG(0))[0], 1.0);
		BOOST_CHECK_NE(cv::sum(finder.get_layersG(0))[0], cv::sum(input)[0]);
		finder.fill(input);
		//Gaussian blur should be normalized
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(finder.get_layersG(1))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(1))[0], cv::sum(finder.get_layersG(2))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(2))[0], cv::sum(finder.get_layersG(3))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(3))[0], cv::sum(finder.get_layersG(4))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(finder.get_layersG(4))[0], 1e-5);
		//The second to last Gaussian layer should be a good approxiation
		//to the input blurred by a two time larger radius than the preblur
		cv::GaussianBlur(input, other, cv::Size(0,0), 2*1.6);
		finder.preblur_and_fill(input);
		BOOST_REQUIRE_CLOSE(cv::sum(finder.get_layersG(3))[0], cv::sum(other)[0], 1e-4);
		BOOST_CHECK_CLOSE(finder.get_layersG(3)(0, 32), other(0, 32), 1e-2);
		images_are_close(finder.get_layersG(3), other, 1e-4);
		//Sum of layers should reconstruct blurred image
		other = finder.get_layersG(0) + finder.get_layers(0);
		BOOST_REQUIRE_CLOSE(cv::sum(finder.get_layersG(1))[0], cv::sum(other)[0], 1e-9);
		images_are_close(finder.get_layersG(1), other);
		other += finder.get_layers(1);
		images_are_close(finder.get_layersG(2), other);
		other += finder.get_layers(2);
		images_are_close(finder.get_layersG(3), other);
		other += finder.get_layers(3);
		images_are_close(finder.get_layersG(4), other);
	}
	BOOST_AUTO_TEST_CASE( single_1Dblob )
	{
		OctaveFinder1D finder(256);
		OctaveFinder::Image input(1, 256);
		input.setTo(0);
		cv::circle(input, cv::Point(128, 0), 3, 1.0, -1);
		for(int i=128-3; i<128+4; ++i)
			BOOST_REQUIRE_CLOSE(input(0, i), 1, 1e-9);
		finder.preblur_and_fill(input);
		//with a radius of 3, the maximum should be in layer 3
		BOOST_CHECK_GE(finder.get_layers(0)(0, 128), finder.get_layers(1)(0, 128));
		BOOST_CHECK_GE(finder.get_layers(1)(0, 128), finder.get_layers(2)(0, 128));
		BOOST_CHECK_GE(finder.get_layers(2)(0, 128), finder.get_layers(3)(0, 128));
		BOOST_CHECK_LE(finder.get_layers(3)(0, 128), finder.get_layers(4)(0, 128));
		finder.initialize_binary();

		//The minimum should be at the center of the circle
		double min_layers[5];
		for (int l=0; l<5; ++l)
			min_layers[l] =  *std::min_element(finder.get_layers(l).begin(), finder.get_layers(l).end());
		double global_min = *std::min_element(min_layers, min_layers+5);
		BOOST_CHECK_CLOSE(finder.get_layers(3)(0, 128), global_min, 1e-9);

		BOOST_CHECK_EQUAL(finder.get_binary(3)(0, 128), true);
		BOOST_CHECK_EQUAL(finder.get_binary(3)(0, 129), false);
		BOOST_CHECK_EQUAL(finder.get_binary(3)(0, 127), false);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(2))[0], 0);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(1))[0], 0);
		BOOST_CHECK_EQUAL(cv::sum(finder.get_binary(3))[0], 1);
		/*cv::namedWindow("truc");
		cv::imshow("truc", 255*finder.get_binary(2));
		cv::waitKey();*/
		//gaussian response
		BOOST_CHECK_LT(finder.gaussianResponse(0, 128, 2.5), finder.gaussianResponse(0, 128, 2.0));
		BOOST_CHECK_LT(finder.gaussianResponse(0, 128, 3.5), finder.gaussianResponse(0, 128, 2.5));
		BOOST_CHECK_LT(finder.gaussianResponse(0, 128, 3.5), finder.gaussianResponse(0, 128, 3.0));
		BOOST_CHECK_GT(finder.gaussianResponse(0, 128, 3.5)-finder.gaussianResponse(0, 128, 2.5), finder.get_layers(3)(0, 128));
		BOOST_CHECK_GT(finder.gaussianResponse(0, 128, 4.5)-finder.gaussianResponse(0, 128, 3.5), finder.get_layers(3)(0, 128));
		//lower bound
		BOOST_CHECK_LT(finder.gaussianResponse(0, 128, 1.5), finder.gaussianResponse(0, 128, 1.0));
		BOOST_CHECK_LT(finder.gaussianResponse(0, 128, 0.5), finder.gaussianResponse(0, 128, 0));
		//further than the top layer
		BOOST_CHECK_LT(finder.gaussianResponse(0, 128, 10.5), finder.gaussianResponse(0, 128, 0));
		BOOST_CHECK_LT(finder.gaussianResponse(0, 128, 15.5), finder.gaussianResponse(0, 128, 5));
	}
	BOOST_AUTO_TEST_CASE( asymetrical )
	{
		OctaveFinder1D finder(32);
		OctaveFinder::Image input(1, 32);
		input.setTo(0);
		std::fill_n(&input(0, 10), 7, 1);
		BOOST_REQUIRE_EQUAL(finder(input).size(), 1);
		input.setTo(0);
		std::fill_n(&input(0, 10), 6, 1);
		BOOST_REQUIRE_EQUAL(finder(input).size(), 1);
	}
	BOOST_AUTO_TEST_CASE( subpix )
	{
		OctaveFinder1D finder(256);
		OctaveFinder::Image input(1, 256);
		input.setTo(0);
		cv::circle(input, cv::Point(100, 0), 3, 1.0, -1);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		BOOST_REQUIRE_EQUAL(finder.get_binary(3)(0, 100), 1);
		std::vector<Center2D> v;
		finder.subpix(v);
		BOOST_REQUIRE_EQUAL(v.size(), 1);
		//x
		BOOST_CHECK_GE(v[0][0], 99);
		BOOST_CHECK_LE(v[0][0], 101);
		//y
		BOOST_CHECK_GE(v[0][1], 0);
		BOOST_CHECK_LT(v[0][1], 1);
		//scale
		BOOST_CHECK_GE(v[0].r, 2);
		BOOST_CHECK_LE(v[0].r, 4);
		//gaussian response
		BOOST_CHECKPOINT("gaussian response");
		BOOST_CHECK_CLOSE(finder.gaussianResponse(0, 100, 2.0), finder.get_layersG(2)(0, 100), 1e-9);
		BOOST_CHECK_CLOSE(finder.gaussianResponse(0, 100, 3.0)-finder.gaussianResponse(0, 100, 2.0), finder.get_layers(2)(0, 100), 1e-9);
		BOOST_CHECK_LT(finder.gaussianResponse(0, 100, 3.0), finder.gaussianResponse(0, 100, 2.0));
		BOOST_CHECK_LT(finder.gaussianResponse(0, 100, 2.5), finder.gaussianResponse(0, 100, 2.0));
		BOOST_CHECK_LT(finder.gaussianResponse(0, 100, 3.5), finder.gaussianResponse(0, 100, 2.5));
		BOOST_CHECK_LT(finder.gaussianResponse(0, 100, 3.5), finder.gaussianResponse(0, 100, 3.0));
		BOOST_CHECK_GT(finder.gaussianResponse(0, 100, 3.5)-finder.gaussianResponse(0, 100, 2.5), finder.get_layers(3)(0, 100));
		BOOST_CHECK_GT(finder.gaussianResponse(0, 100, 4.5)-finder.gaussianResponse(0, 100, 3.5), finder.get_layers(3)(0, 100));
		BOOST_CHECKPOINT("sublayers 1");
		boost::array<double,8> sublayerG;
		for(size_t u = 0;u < sublayerG.size(); ++u)
			sublayerG[u] = finder.gaussianResponse(0, 100, 2 + 0.5*u);
		boost::array<double,5> a;
		for(size_t u =0; u<a.size();++u)
			a[u] = sublayerG[u+2] - sublayerG[u];
		BOOST_CHECK_LT(a[0], 0);
		BOOST_CHECK_LT(a[1], 0);
		BOOST_CHECK_LT(a[2], 0);
		BOOST_CHECK_LT(a[3], 0);
		BOOST_CHECK_LT(a[4], 0);
		BOOST_REQUIRE_NE(a[4]-2*a[2]+a[0], 0);
		double ds = (-a[4] + 8*a[3] - 8*a[1] + a[0])/6.0 /(a[4]-2*a[2]+a[0]), z = 3-ds;
		BOOST_CHECK_GT(ds, -0.5);
		BOOST_CHECK_LT(ds, 0.5);
		BOOST_CHECKPOINT("sublayers 2");
		for(size_t u = 0;u < sublayerG.size(); ++u)
			sublayerG[u] = finder.gaussianResponse(0, 100, z - 1 + 0.5*u);
		for(size_t u =0; u<a.size();++u)
			a[u] = sublayerG[u+2] - sublayerG[u];
		BOOST_CHECK_LT(a[0], 0);
		BOOST_CHECK_LT(a[2], 0);
		BOOST_CHECK_LT(a[4], 0);
		BOOST_REQUIRE_NE(a[4]-2*a[2]+a[0], 0);
		ds = (-a[4] + 8*a[3] - 8*a[1] + a[0])/6.0 /(a[4]-2*a[2]+a[0]);
		BOOST_CHECK_GT(z-ds, 2.5);
		BOOST_CHECK_LT(z-ds, 3.5);

	}
	BOOST_AUTO_TEST_CASE( subpix_relative_positions )
	{
		OctaveFinder1D finder(32);
		//cv cannot draw circle positions better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(1, 32*32), small_input(1,32);
		std::vector<Center2D> v;
		for(int i=0; i<6; ++i)
		{
			input.setTo(0);
			cv::circle(input, cv::Point(32*(16+pow(0.5, i+1)), 0), 32*3, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			finder.preblur_and_fill(small_input);
			finder.initialize_binary();
			std::vector<Center2D> v_s;
			finder.subpix(v_s);
			std::copy(v_s.begin(), v_s.end(), std::back_inserter(v));
		}
		BOOST_REQUIRE_EQUAL(v.size(), 6);
		for(int i=0; i<5; ++i)
			BOOST_WARN_MESSAGE(
					v[i][0] > v[i+1][0],
					"1D spatial resolution is larger than 1/" << pow(2,i+1) <<"th of a pixel"
					);
		BOOST_CHECK_GT(v[2][0], v[3][0]);
	}
	BOOST_AUTO_TEST_CASE( subpix_relative_sizes )
	{
		OctaveFinder1D finder(64, 3);
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(1, 16*64), small_input(1,64);
		std::vector<Center2D> v;
		std::ofstream f("1d_relative_sizes.out");
		for(int i=0; i<32; ++i)
		{
			input.setTo(0);
			cv::circle(input, cv::Point(16*32, 0), 16*2+i, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			finder.preblur_and_fill(small_input);
			finder.initialize_binary();
			std::vector<Center2D> v_s;
			finder.subpix(v_s);
			BOOST_REQUIRE_EQUAL(v_s.size(), 1);
			std::copy(v_s.begin(), v_s.end(), std::back_inserter(v));
			for(size_t c=0; c<v_s.size(); ++c)
				finder.scale(v_s[c]);
			f<< (2+i/16.0) <<"\t" << v_s[0].r <<"\n";
		}

		for (int i=1; i<7;++i)
			BOOST_WARN_MESSAGE(v[0].r< v[7-i].r, "resolution in size is 1/"<<(8*(7-i))<< "th of a scale");

		for(size_t c=0; c<v.size(); ++c)
			finder.scale(v[c]);
		for (size_t i=0; i<v.size();++i)
			BOOST_REQUIRE_CLOSE(v[i].r, (2+i/16.0), 2.5);
	}

	BOOST_AUTO_TEST_CASE( octave_minimum_detected_size )
	{
		OctaveFinder1D finder(32);
		cv::Mat_<uchar>input(1, 16*32), small_input(1,32);
		int i = 0;
		std::vector<Center2D> v(1);
		while(3-0.01*i>0 && v.size()>0)
		{
			i++;
			input.setTo(0);
			cv::circle(input, cv::Point(16*16, 0), 16*(3-0.01*i), 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			v =	finder(small_input, true);
		}
		BOOST_WARN_MESSAGE(false, "Cannot detect sizes below "<< (3-0.01*(i-1))<<" pixels in 1D");
	}

	BOOST_AUTO_TEST_CASE( octave_minimum_detector_size )
	{
		int i = 32;
		std::vector<Center2D> v(1);
		while(i>0 && v.size()>0)
		{
			i--;
			OctaveFinder1D finder(i);
			cv::Mat_<uchar>input(1, 16*i), small_input(1,i);
			input.setTo(0);
			cv::circle(input, cv::Point(8*i, 0), 16*2, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			v =	finder(small_input, true);
		}
		BOOST_WARN_MESSAGE(false, "An 1D octave detector smaller than "<< (i+1)<<" pixels cannot detect anything");
	}
	BOOST_AUTO_TEST_CASE( octave_minimum_signal_intensity )
	{
		OctaveFinder1D finder(16);
		OctaveFinder::Image  input(1,16);
		int i = 0;
		std::vector<Center2D> v(1);
		while(i<10 && v.size()>0)
		{
			i++;
			input.setTo(0);
			cv::circle(input, cv::Point(8, 0), 4, pow(0.5, i), -1);
			v =	finder(input, true);
		}
		BOOST_WARN_MESSAGE(i>10, "Cannot detect intensities below "<< pow(0.5, i-1)<<" in 1D");
	}
BOOST_AUTO_TEST_SUITE_END() //octave
BOOST_AUTO_TEST_SUITE( multiscale )

	BOOST_AUTO_TEST_CASE( multiscale_constructor_1D )
	{
		MultiscaleFinder1D finder;
		BOOST_REQUIRE_EQUAL(finder.get_n_octaves(), 6);
		BOOST_CHECK_EQUAL(finder.get_octave(1).get_width(), 1);
		BOOST_CHECK_EQUAL(finder.get_octave(1).get_height(), 256);
		BOOST_CHECK_EQUAL(finder.get_octave(1).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(finder.get_octave(0).get_width(), 1);
		BOOST_CHECK_EQUAL(finder.get_octave(0).get_height(), 512);
		BOOST_CHECK_EQUAL(finder.get_octave(0).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(finder.get_octave(2).get_width(), 1);
		BOOST_CHECK_EQUAL(finder.get_octave(2).get_height(), 128);
		BOOST_CHECK_EQUAL(finder.get_octave(2).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(finder.get_octave(3).get_width(), 1);
		BOOST_CHECK_EQUAL(finder.get_octave(3).get_height(), 64);
		BOOST_CHECK_EQUAL(finder.get_octave(3).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(finder.get_octave(4).get_width(), 1);
		BOOST_CHECK_EQUAL(finder.get_octave(4).get_height(), 32);
		BOOST_CHECK_EQUAL(finder.get_octave(4).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(finder.get_octave(5).get_width(), 1);
		BOOST_CHECK_EQUAL(finder.get_octave(5).get_height(), 16);
		BOOST_CHECK_EQUAL(finder.get_octave(5).get_n_layers(), 3);
	}
	BOOST_AUTO_TEST_CASE( multiscale1D_single_blob )
	{
		MultiscaleFinder1D finder;
		OctaveFinder::Image input(1, 256);
		input.setTo(0);
		cv::circle(input, cv::Point(128, 0), 3, 1.0, -1);
		std::vector<Center2D> v = finder(input);
		//the first octave's should be the same as the one of an Octavefinder of the same size
		OctaveFinder1D ofinder;
		ofinder.preblur_and_fill(input);
		BOOST_CHECK_CLOSE(finder.get_octave(1).get_layers(0)(0,128), ofinder.get_layers(0)(0,128), 1e-9);
		images_are_close(finder.get_octave(1).get_layers(0), ofinder.get_layers(0), 1e-4);
		//with a radius of 4, the maximum should be in layer 2 of octave 1 and nowhere else
		for(size_t o=0; o<finder.get_n_octaves(); ++o)
			for(size_t l=1; l<finder.get_n_layers()+1; ++l)
			{
				const int u = cv::sum(finder.get_octave(o).get_binary(l))[0];
				BOOST_CHECK_MESSAGE(
						u == ((o==1 && l==3)?1:0),
						"Octave "<<o<<" layer "<<l<<" has "<< u <<" center"<<((u>1)?"s":"")
				);
			}
		BOOST_REQUIRE_EQUAL(v.size(), 1);
		BOOST_CHECK_CLOSE(v[0][0], 128, 10);
		BOOST_CHECK_CLOSE(v[0][1], 0, 10);
		BOOST_CHECK_CLOSE(v[0].r, 3.5, 2.5);
	}
	BOOST_AUTO_TEST_CASE( multiscale_relative_sizes1D )
	{
		const int s = 5;
		MultiscaleFinder1D finder(pow(2, s+1));
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(1, 32*pow(2, s+1)), small_input(1, pow(2, s+1));
		std::vector<Center2D> v;
		std::ofstream f("multiscale1D_relative_sizes.out");
		std::set<int> large_radii;
		for(int k=1; k<s-1; ++k)
			for(int i=0; i<32; ++i)
				large_radii.insert(48*pow(2, k-1)+i*pow(2, k));
		double radius;
		for(std::set<int>::const_iterator lr = large_radii.begin(); lr != large_radii.end(); ++lr)
		{
			input.setTo(0);
			radius =  *lr/32.0;
			BOOST_TEST_CHECKPOINT("radius = "<<radius<<" call");
			cv::circle(input, cv::Point(32*pow(2, s), 0), *lr, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			std::vector<Center2D> v_s = finder(small_input);
			BOOST_TEST_CHECKPOINT("radius = "<<radius<<" size");
			BOOST_CHECK_MESSAGE(
					v_s.size()==1,
					""<<v_s.size()<<" centers for input radius "<<radius
					);
			if(v_s.size()>1)
			{
				for(size_t o=0; o<finder.get_n_octaves(); ++o)
				{
					std::cout<<finder.get_octave(o).get_nb_centers()<<" centers in octave "<<o<<" at"<<std::endl;
					for(size_t l=1; l<finder.get_n_layers()+1; ++l)
						for(size_t i=0; i<finder.get_octave(o).get_height(); ++i)
							if(finder.get_octave(o).get_binary(l)(0,i))
							{
								std::cout<<"\tl="<<l<<"\ti="<<i<<std::endl;
							}
				}
			}
			BOOST_TEST_CHECKPOINT("radius = "<<radius<<" copy");
			std::copy(v_s.begin(), v_s.end(), std::back_inserter(v));
			for(size_t j=0; j<v_s.size(); ++j)
				f << radius << "\t" << v_s[j].r<< "\t" << v_s[j][0]<< "\n";
		}
		f<<std::endl;
		BOOST_REQUIRE_EQUAL(v.size(), large_radii.size());
	}
	BOOST_AUTO_TEST_CASE( multiscale1D_close_neighbours )
	{
		MultiscaleFinder1D finder(64);
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(1, 32*64), small_input(1, 64);
		std::vector<Center2D> v_s(2);
		std::ofstream f("close_neighbours1D.out");
		for(int i = 32*32; i>8*32 && v_s.size()==2; --i)
		{
			const double distance = i/32.0;
			BOOST_TEST_CHECKPOINT("distance = "<<distance);
			input.setTo(0);
			cv::circle(input, cv::Point(32*16, 0), 32*3, 255, -1);
			cv::circle(input, cv::Point(32*16+i, 0), 32*3, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			finder.get_centers(small_input, v_s);
			if(v_s.size()!=2)
			{
				BOOST_WARN_MESSAGE(false, "Cannot differentiate step-like 1D blobs with edges closer than "<<distance/3-2<<" radius");
				std::ofstream out("close_neighbours1D.layersG");
				for(size_t j=0; j<64; ++j)
				{
					for(size_t l=0; l<finder.get_n_layers()+3; ++l)
						out<<finder.get_octave(1).get_layersG(l)(0, j)<<"\t";
					out<<"\n";
				}
				break;
			}
			BOOST_CHECK_CLOSE(v_s[0][0], 16, distance<16?6:2);
			BOOST_REQUIRE_CLOSE(v_s[1][0], 16+distance, distance<16?6:2);
			BOOST_CHECK_CLOSE(v_s[0].r, 3, distance<16?6:2);
			BOOST_CHECK_CLOSE(v_s[1][0] - v_s[0][0], distance, distance<9.03125?21:2);
			f << distance << "\t" << v_s[0][0] << "\t" << v_s[1][0] << "\t" << v_s[0].r << "\n";
		}
	}
BOOST_AUTO_TEST_SUITE_END() //multiscale

BOOST_AUTO_TEST_SUITE_END() //1D

BOOST_AUTO_TEST_SUITE( trajectories )
	BOOST_AUTO_TEST_CASE( traj )
	{
		Traj traj(-1, 5);
		BOOST_CHECK_EQUAL(traj.get_start(), -1);
		BOOST_CHECK_EQUAL(traj.get_finish(), 0);
		BOOST_CHECK_EQUAL(traj[-1], 5);
		traj.push_back(2);
		BOOST_CHECK_EQUAL(traj.get_start(), -1);
		BOOST_CHECK_EQUAL(traj.get_finish(), 1);
		BOOST_CHECK_EQUAL(traj[-1], 5);
		BOOST_CHECK_EQUAL(traj[0], 2);
	}
	BOOST_AUTO_TEST_CASE( trajindex )
	{
		TrajIndex ti(2);
		BOOST_REQUIRE_EQUAL(ti.nbFrames(), 1);
		BOOST_REQUIRE_EQUAL(ti.size(), 2);
		BOOST_REQUIRE_EQUAL(ti.getInverse(0).size(), 2);
		for(size_t p = 0; p<ti.getInverse(0).size(); ++p)
			BOOST_CHECK_EQUAL(ti.getInverse(0)[p], p);
		for(size_t tr=0; tr<ti.size(); ++tr)
			BOOST_CHECK_EQUAL(ti[tr][0], tr);
		std::vector<size_t> from(6), to(6);
		std::fill_n(from.begin(), 3, 0);
		std::fill_n(from.rbegin(), 3, 1);
		to[0] = 0;
		to[1] = 1;
		to[2] = 2;
		to[3] = 0;
		to[4] = 1;
		to[5] = 2;
		std::vector<double> distances(6);
		distances[0] = 1,1;
		distances[1] = 1.2;
		distances[2] = 0.9;
		distances[3] = 0.1;
		distances[4] = 1.1;
		distances[5] = 1.0;
		//case where we have more new positions than old ones
		ti.add_Frame(4, distances, from, to);
		BOOST_REQUIRE_EQUAL(ti.nbFrames(), 2);
		BOOST_REQUIRE_EQUAL(ti.getInverse(1).size(), 4);
		BOOST_REQUIRE_EQUAL(ti[1].size(), 2);
		BOOST_CHECK_EQUAL(ti[1][1], 0);
		BOOST_CHECK_EQUAL(ti.getTraj(1,0), 1);
		BOOST_REQUIRE_EQUAL(ti[0].size(), 2);
		BOOST_CHECK_EQUAL(ti[0][1], 2);
		BOOST_CHECK_EQUAL(ti.getTraj(1,2), 0);
		BOOST_REQUIRE(ti[2].exist(1));
		BOOST_CHECK_EQUAL(ti[2][1], 1);
		BOOST_CHECK_EQUAL(ti.getTraj(1,1), 2);
		BOOST_REQUIRE(ti[3].exist(1));
		BOOST_CHECK_EQUAL(ti[3][1], 3);
		BOOST_CHECK_EQUAL(ti.getTraj(1,3), 3);
		BOOST_REQUIRE_EQUAL(ti.size(), 4);

		//case where we have less positions than old ones
		ti.add_Frame(3, distances, from, to);
		BOOST_REQUIRE_EQUAL(ti.nbFrames(), 3);
		BOOST_REQUIRE_EQUAL(ti.size(), 5);
		BOOST_REQUIRE_EQUAL(ti.getInverse(2).size(), 3);
		BOOST_CHECK_EQUAL(ti[0].size(), 2);
		BOOST_REQUIRE_EQUAL(ti[1].size(), 3);
		BOOST_CHECK_EQUAL(ti[1][2], 2);
		BOOST_CHECK_EQUAL(ti.getTraj(2,2), 1);
		BOOST_CHECK(ti[2].exist(2));
		BOOST_CHECK_EQUAL(ti[2][2], 0);
		BOOST_CHECK_EQUAL(ti.getTraj(2,0), 2);
		BOOST_CHECK(!ti[3].exist(2));
		BOOST_CHECK(ti[4].exist(2));
		BOOST_CHECK_EQUAL(ti[4][2], 1);
		BOOST_CHECK_EQUAL(ti.getTraj(2,1), 4);

		//load a smaller frame
		boost::array<size_t, 2> fr0 = {{2, 1}};
		ti.add_Frame(fr0.begin(), fr0.end());
		BOOST_REQUIRE_EQUAL(ti.nbFrames(), 4);
		BOOST_REQUIRE_EQUAL(ti.size(), 5);
		BOOST_REQUIRE_EQUAL(ti.getInverse(3).size(), 2);
		BOOST_CHECK_EQUAL(ti[0].size(), 2);
		BOOST_REQUIRE_EQUAL(ti[1].size(), 4);
		BOOST_CHECK_EQUAL(ti[1][3], 1);
		BOOST_CHECK_EQUAL(ti.getTraj(3,1), 1);
		BOOST_CHECK(ti[2].exist(3));
		BOOST_CHECK_EQUAL(ti[2][3], 0);
		BOOST_CHECK_EQUAL(ti.getTraj(3,0), 2);
		BOOST_CHECK(!ti[3].exist(3));
		BOOST_CHECK(!ti[4].exist(3));

		//load a larger frame
		boost::array<size_t, 3> fr1 = {{2, 1, 5}};
		ti.add_Frame(fr1.begin(), fr1.end());
		BOOST_REQUIRE_EQUAL(ti.nbFrames(), 5);
		BOOST_REQUIRE_EQUAL(ti.size(), 6);
		BOOST_REQUIRE_EQUAL(ti.getInverse(4).size(), 3);
		BOOST_REQUIRE_EQUAL(ti[1].size(), 5);
		BOOST_CHECK_EQUAL(ti[1][4], 1);
		BOOST_CHECK_EQUAL(ti.getTraj(4,1), 1);
		BOOST_REQUIRE(ti[2].exist(4));
		BOOST_CHECK_EQUAL(ti[2][4], 0);
		BOOST_CHECK_EQUAL(ti.getTraj(4,0), 2);
		BOOST_REQUIRE(ti[5].exist(4));
		BOOST_CHECK_EQUAL(ti[5][4], 2);
		BOOST_CHECK_EQUAL(ti.getTraj(4,2), 5);

		//frame with no link
		ti.add_Frame(4, std::vector<double>(), std::vector<size_t>(), std::vector<size_t>());
		BOOST_REQUIRE_EQUAL(ti.nbFrames(), 6);
		BOOST_REQUIRE_EQUAL(ti.size(), 10);
		BOOST_REQUIRE_EQUAL(ti.getInverse(5).size(), 4);
		BOOST_REQUIRE_EQUAL(ti[1].size(), 5);
		BOOST_CHECK_EQUAL(ti[1][4], 1);
		BOOST_CHECK_EQUAL(ti.getTraj(4,1), 1);
		BOOST_REQUIRE(!ti[2].exist(5));
		BOOST_CHECK_EQUAL(ti[2][4], 0);
		BOOST_CHECK_EQUAL(ti.getTraj(4,0), 2);
		BOOST_REQUIRE(!ti[5].exist(5));
		BOOST_CHECK_EQUAL(ti[5][4], 2);
		BOOST_CHECK_EQUAL(ti.getTraj(4,2), 5);

	}
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE( Overlap )
	BOOST_AUTO_TEST_CASE(Overlap)
	{
		typedef RStarTree<size_t, 2, 4, 32, double> RTree;
		std::vector<Center2D> centers;
		std::auto_ptr<RTree> tree = removeOverlapping(centers);
		BOOST_REQUIRE(centers.empty());
		Center2D c(0, 1, 0);
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 1);
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 1);
		c.intensity = -1.0;
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 1);
		BOOST_CHECK_CLOSE(centers[0].intensity, -1, 1e-9);
		c[0] = 3;
		c.intensity = -2.0;
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 2);
		BOOST_CHECK_CLOSE(centers[0].intensity, -2, 1e-9);
		BOOST_CHECK_CLOSE(centers[1].intensity, -1, 1e-9);
		c[0] = 2;
		c.intensity = 0;
		c.r = 1.1;
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 2);
		BOOST_CHECK_CLOSE(centers[0].intensity, -2, 1e-9);
		BOOST_CHECK_CLOSE(centers[1].intensity, -1, 1e-9);
		c.intensity = -1.5;
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 2);
		BOOST_CHECK_CLOSE(centers[0].intensity, -2, 1e-9);
		BOOST_CHECK_CLOSE(centers[1].intensity, -1, 1e-9);
		c.intensity = -3.0;
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 1);
		BOOST_CHECK_CLOSE(centers[0].intensity, -3, 1e-9);
	}
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE( Reconstruction )
	BOOST_AUTO_TEST_CASE( add_frame )
	{
		Reconstructor rec;
		BOOST_REQUIRE(rec.empty());
		Reconstructor::Frame centers(1);
		centers.back().r=1;
		rec.push_back(centers);
		BOOST_REQUIRE_EQUAL(rec.size(), 1);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
		BOOST_CHECK_CLOSE(rec.get_clusters().front().back().r, 1.0, 1e-9);
		rec.clear();
		BOOST_REQUIRE(rec.empty());
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 0);
		rec.push_back(centers);
		BOOST_REQUIRE_EQUAL(rec.size(), 1);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
		BOOST_CHECK_CLOSE(rec.get_clusters().front().back().r, 1.0, 1e-9);
		rec.push_back(centers);
		BOOST_REQUIRE_EQUAL(rec.size(), 2);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
		BOOST_CHECK_CLOSE(rec.get_clusters().front().back().r, 1.0, 1e-9);
		Center2D c(10.0, 0.5);
		centers.push_back(c);
		rec.push_back(centers);
		BOOST_REQUIRE_EQUAL(rec.size(), 3);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 2);
		BOOST_CHECK_CLOSE(rec.get_clusters().front().back().r, 1.0, 1e-9);
		BOOST_CHECK_CLOSE(rec.get_clusters().back().back().r, 0.5, 1e-9);
		centers[0].r = 0.5;
		centers[1].r = 1.0;
		rec.push_back(centers);
		BOOST_REQUIRE_EQUAL(rec.size(), 4);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 2);
		BOOST_CHECK_CLOSE(rec.get_clusters().front().back().r, 0.5, 1e-9);
		BOOST_CHECK_CLOSE(rec.get_clusters().back().back().r, 1.0, 1e-9);
		//More centers in a frame
		centers.push_back(Center2D(4,3));
		rec.push_back(centers);
		BOOST_REQUIRE_EQUAL(rec.size(), 5);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 3);
		BOOST_CHECK_CLOSE(rec.get_clusters().front().back().r, 0.5, 1e-9);
		BOOST_CHECK_CLOSE(rec.get_clusters()[1].back().r, 1.0, 1e-9);
		BOOST_CHECK_CLOSE(rec.get_clusters().back().back().r, 3.0, 1e-9);
		//less centers in a frame
		centers.erase(centers.begin());
		rec.push_back(centers);
		BOOST_REQUIRE_EQUAL(rec.size(), 6);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 3);
		BOOST_REQUIRE_CLOSE(rec.get_clusters().front().back()[2], 4, 1e-9);
		BOOST_CHECK_CLOSE(rec.get_clusters().front().back().r, 0.5, 1e-9);
		BOOST_CHECK_CLOSE(rec.get_clusters()[1].back().r, 1.0, 1e-9);
		BOOST_CHECK_CLOSE(rec.get_clusters().back().back().r, 3.0, 1e-9);

	}
	BOOST_AUTO_TEST_CASE( far )
	{
		//particles too far away
		Reconstructor rec;
		Reconstructor::Frame centers(1);
		centers.back().r=1;
		rec.push_back(centers);
		Center2D c(10.0, 0.5);
		centers.assign(1, c);
		rec.push_back(centers);
		BOOST_REQUIRE_EQUAL(rec.size(), 2);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 2);
		BOOST_CHECK_CLOSE(rec.get_clusters().front().back().r, 1.0, 1e-9);
		BOOST_CHECK_CLOSE(rec.get_clusters().back().back().r, 0.5, 1e-9);
		//particles too far away, with tolerance
		rec.clear();
		centers = Reconstructor::Frame(1);
		centers.back().r=1;
		rec.push_back(centers);
		centers.assign(1, c);
		rec.push_back(centers, 10);
		BOOST_REQUIRE_EQUAL(rec.size(), 2);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
		BOOST_CHECK_CLOSE(rec.get_clusters().front().front().r, 1.0, 1e-9);
		BOOST_CHECK_CLOSE(rec.get_clusters().front().back().r, 0.5, 1e-9);
	}
	BOOST_AUTO_TEST_CASE( cluster_split )
	{
		//no need to split
		Reconstructor rec;
		Reconstructor::Frame centers(1, Center2D(0, 1));
		for(size_t i=0; i<16; ++i)
			rec.push_back(centers);
		BOOST_REQUIRE_EQUAL(rec.size(), 16);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
		rec.split_clusters();
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
		//need to split
		rec.clear();
		for(size_t i=0; i<6; ++i)
			rec.push_back(centers);
		for(size_t i=0; i<7; ++i)
		{
			centers[0][0] += 0.1;
			rec.push_back(centers);
		}
		for(size_t i=0; i<7; ++i)
			rec.push_back(centers);
		BOOST_REQUIRE_EQUAL(rec.size(), 20);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
		rec.split_clusters();
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 2);
		BOOST_CHECK_GE(rec.get_clusters().front().back()[2], 6);
		BOOST_CHECK_LE(rec.get_clusters().back().front()[2], 10);
		//need to split, with magnitude contrast
		rec.clear();
		for(size_t i=0; i<6; ++i)
			rec.push_back(centers);
		for(size_t i=0; i<7; ++i)
		{
			centers[0][0] += 0.05;
			rec.push_back(centers);
		}
		for(size_t i=0; i<6; ++i)
			rec.push_back(centers);
		for(size_t i=0; i<7; ++i)
		{
			centers[0][0] += 0.1;
			rec.push_back(centers);
		}
		for(size_t i=0; i<7; ++i)
			rec.push_back(centers);
		BOOST_REQUIRE_EQUAL(rec.size(), 33);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
		rec.split_clusters();
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 3);
		BOOST_CHECK_GE(rec.get_clusters().front().back()[2], 6);
		BOOST_CHECK_LE(rec.get_clusters().back().front()[2], 13);
		BOOST_CHECK_GE(rec.get_clusters().back().back()[2], 19);
		BOOST_CHECK_LE(rec.get_clusters()[1].front()[2], 26);
	}
	BOOST_AUTO_TEST_CASE( one_sphere )
	{
		Reconstructor rec;
		std::ofstream f("one_sphere.out");
		for(int z0=0; z0<10; ++z0)
		{
			//slice a sphere of radius 4 centeres on 4+z0/10
			const double pos = 4.0 + z0 / 10.0;
			for(int z=0; z<10; ++z)
			{
				const double radsq =  4*4-pow(z - pos, 2);
				if(radsq >= 0)
					rec.push_back(Reconstructor::Frame(1, Center2D(0, sqrt(radsq))));
				else
					rec.push_back(Reconstructor::Frame());
			}
			BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
			std::deque<Center3D> centers;
			rec.get_blobs(centers);
			BOOST_REQUIRE_EQUAL(centers.size(), 1);
			BOOST_CHECK_CLOSE(centers.front()[0], 0, 1e-9);
			BOOST_CHECK_CLOSE(centers.front()[1], 0, 1e-9);
			BOOST_CHECK_CLOSE(centers.front()[2], pos, 2);
			BOOST_CHECK_CLOSE(centers.front().r, 4, 2);
			f<<pos<<"\t"<<centers.front()[2]<<"\n";
			rec.clear();
		}
	}
	/*BOOST_AUTO_TEST_CASE( two_identical_spheres )
	{
		Reconstructor rec;
		//std::ofstream f("one_sphere.out");
		for(int z0=1; z0<10; ++z0)
		{
			//slice a sphere of radius 4 centeres on 4
			for(int z=0; z<8; ++z)
			{
				const double radsq =  4*4-pow(z - 4, 2);
				rec.push_back(Reconstructor::Frame(1, Center2D(0, sqrt(radsq))));
			}
			rec.push_back(Reconstructor::Frame(1, Center2D(0, 1)));
			//slice a sphere of radius 4 centeres on 8+z0/10
			const double pos = 4.0 + z0 / 10.0;
			for(int z=0; z<10; ++z)
			{
				const double radsq =  4*4-pow(z - pos, 2);
				if(radsq >= 0)
					rec.push_back(Reconstructor::Frame(1, Center2D(0, sqrt(radsq))));
				else
					rec.push_back(Reconstructor::Frame(1, Center2D(0, 1)));
			}
			BOOST_REQUIRE_MESSAGE(rec.nb_cluster()==1, ""<<rec.nb_cluster()<<" clusters at pos="<<pos);
			std::deque<Center3D> centers;
			rec.get_blobs(centers);
			BOOST_REQUIRE_MESSAGE(centers.size()==2, ""<<centers.size()<<" centers at pos="<<pos);
			BOOST_CHECK_CLOSE(centers.front()[0], 0, 1e-9);
			BOOST_CHECK_CLOSE(centers.front()[1], 0, 1e-9);
			BOOST_CHECK_CLOSE(std::min(centers.front()[2], centers.back()[2]), 4, 2);
			BOOST_CHECK_CLOSE(centers.front().r, 4, 2);
			BOOST_CHECK_CLOSE(centers.back()[0], 0, 1e-9);
			BOOST_CHECK_CLOSE(centers.back()[1], 0, 1e-9);
			BOOST_CHECK_CLOSE(std::max(centers.front()[2], centers.back()[2]), pos+7, 2);
			BOOST_CHECK_CLOSE(centers.back().r, 4, 2);
			//f<<pos<<"\t"<<centers.front()[2]<<"\n";
			rec.clear();
		}
	}*/
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE( Lif )
	BOOST_AUTO_TEST_CASE( export_z_scan )
	{
		LifReader reader("/home/mathieu/Code_data/liftest/Tsuru11dm_phi=52.53_J36.lif");
		LifSerie &serie = reader.getSerie(0);
		std::vector<size_t> dims = serie.getSpatialDimensions();
		BOOST_CHECK_GE(dims.size(), 2);
		cv::Mat_<unsigned char> slice(dims[0], dims[1]);
		cv::Mat color_slice(slice.size(), CV_8UC3);
		BOOST_CHECK_GT(dims.size(), 2);
		const size_t total_z = dims.size()>2 ? dims[2] : 1;
		BOOST_CHECK_EQUAL(total_z, 256);
		cv::VideoWriter w("z_scan.avi", CV_FOURCC('D', 'I', 'V', 'X'), 16, slice.size(), true);
		for(size_t z=0; z<total_z; ++z)
		{
			serie.fill2DBuffer(static_cast<void*>(slice.data), 0, z);
			unsigned char *c = color_slice.data;
			cv::Mat_<unsigned char>::const_iterator b = slice.begin();
			while(b!=slice.end())
			{
				*c++ = *b;
				*c++ = *b;
				*c++ = *b++;
			}
			w << color_slice;
		}
	}
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE( LifTrack )
	BOOST_AUTO_TEST_CASE( fill_one_slice )
	{
		LifReader reader("/home/mathieu/Code_data/liftest/Tsuru11dm_phi=52.53_J36.lif");
		LocatorFromLif locator(&reader.getSerie(0));
		BOOST_CHECK_EQUAL(cv::sum(locator.get_slice())[0], 0);
		locator.fill_next_slice();
		BOOST_CHECK_EQUAL(cv::sum(locator.get_slice())[0], 5357342);
		BOOST_CHECK_EQUAL(locator.get_reconstructor().size(), 1);
		locator.clear();
		BOOST_CHECK_EQUAL(cv::sum(locator.get_slice())[0], 0);
		BOOST_CHECK_EQUAL(locator.get_z(), 0);
	}
	BOOST_AUTO_TEST_CASE( fill_one_stack )
	{
		LifReader reader("/home/mathieu/Code_data/liftest/Tsuru11dm_phi=52.53_J36.lif");
		LocatorFromLif locator(&reader.getSerie(0));
		{
			std::cout<<"fill_time_step ";
			boost::progress_timer ti;
			locator.fill_time_step();
		}
		BOOST_CHECK_EQUAL(locator.get_z(), 256);

		//export to VTK format for human-eye comparison
		typedef Reconstructor::Cluster Cluster;
		typedef std::deque<Cluster> Clusters;
		const Clusters & clusters = locator.get_reconstructor().get_clusters();
		Clusters::const_iterator cl_min=clusters.begin();
		int count_file =0;
		//their may be too many clusters than the maximum of int type, so we have to divide in several files
		while(cl_min!=clusters.end())
		{
			Clusters::const_iterator cl_max = cl_min;
			size_t nb = 0, nb_cl=0;
			while(cl_max!=clusters.end() && (nb + cl_max->size() + nb_cl +1 < (size_t)std::numeric_limits<int>::max()))
			{
				nb += cl_max->size();
				nb_cl++;
				cl_max++;
			}
			std::cout<<"file "<<count_file<<" has "<<nb_cl<<" clusters and "<<nb<<" particles"<<std::endl;
			std::ostringstream os;
			os<<"fill_one_stack_clusters"<< count_file++ <<".vtk";
			std::ofstream f_cl(os.str().c_str());
			f_cl<<"# vtk DataFile Version 3.0\n"
							"fill_one_stack\n"
							"ASCII\n"
							"DATASET POLYDATA\n"
							"POINTS "<<nb<<" double\n";
			for(Clusters::const_iterator cl=cl_min;cl!=cl_max;++cl)
				for(Cluster::const_iterator p = cl->begin(); p!=cl->end(); ++p)
				{
					for(size_t d=0;d<3;++d)
						f_cl<<(*p)[d]<<" ";
					f_cl<<"\n";
				}
			f_cl<< "LINES "<<nb_cl<<" "<<nb + nb_cl<<"\n";
			size_t l=0;
			for(Clusters::const_iterator cl=cl_min;cl!=cl_max;++cl)
			{
				f_cl<<cl->size()<<" ";
				for(size_t p=0; p<cl->size();++p)
					f_cl<< l++ <<" ";
				f_cl<<"\n";
			}
			f_cl<<"POINT_DATA "<<nb<<"\n"
					"SCALARS r double\n"
					"LOOKUP_TABLE default\n";
			for(Clusters::const_iterator cl=cl_min;cl!=cl_max;++cl)
				for(Cluster::const_iterator p = cl->begin(); p!=cl->end(); ++p)
					f_cl<< p->r <<"\n";
			f_cl.close();
			cl_min = cl_max;
		}

		LocatorFromLif::Centers centers;
		{
			std::cout<<"get_centers ";
			boost::progress_timer ti;
			locator.get_centers(centers);
		}
		std::cout<<std::endl;
		BOOST_CHECK_EQUAL(locator.get_t(), 1);
		BOOST_REQUIRE(!centers.empty());
		BOOST_REQUIRE_EQUAL(centers.size(), 7946);
		//export to VTK format for human-eye comparison
		std::ofstream f_rec("fill_one_stack_reconstructed.vtk");
		f_rec<<"# vtk DataFile Version 3.0\n"
				"fill_one_stack\n"
				"ASCII\n"
				"DATASET POLYDATA\n"
				"POINTS "<<centers.size()<<" double\n";
		for(LocatorFromLif::Centers::const_iterator p=centers.begin();p!=centers.end();++p)
		{
			for(size_t d=0;d<3;++d)
				f_rec<<(*p)[d]<<" ";
			f_rec<<"\n";
		}
		f_rec<<"POINT_DATA "<<centers.size()<<"\n"
				"SCALARS r double\n"
				"LOOKUP_TABLE default\n";
		for(LocatorFromLif::Centers::const_iterator p=centers.begin();p!=centers.end();++p)
			f_rec<< p->r <<"\n";
		f_rec.close();
	}
BOOST_AUTO_TEST_SUITE_END()
