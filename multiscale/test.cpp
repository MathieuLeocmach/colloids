#define BOOST_TEST_MODULE multiscale test
#define BOOST_TEST_DYN_LINK

#include "multiscalefinder.hpp"
#include <fstream>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
#include <boost/array.hpp>

using namespace Colloids;

void images_are_close(const cv::Mat &a, const cv::Mat &b, double precision=1e-9)
{
	cv::Mat_<double> M = cv::abs(a)+cv::abs(b);
	double peak = *std::max_element(M.begin(), M.end());
	BOOST_REQUIRE_CLOSE(cv::sum(a)[0], cv::sum(b)[0], precision/peak);
	cv::Mat_<double> diff = cv::abs(a-b) / peak;
	cv::MatConstIterator_<double> u = diff.begin();
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
		BOOST_REQUIRE_CLOSE(cv::sum(finder.get_layersG(3))[0], cv::sum(other)[0], 1e-9);
		BOOST_CHECK_CLOSE(finder.get_layersG(3)(32,32), other(32,32), 1e-2);
		images_are_close(finder.get_layersG(3), other, 1e-4);
		//Sum of layers should reconstruct blurred image
		other = finder.get_layersG(0) + finder.get_layers(0);
		BOOST_REQUIRE_CLOSE(cv::sum(finder.get_layersG(1))[0], cv::sum(other)[0], 1e-9);
		images_are_close(finder.get_layersG(1), other, 1e-9);
		other += finder.get_layers(1);
		images_are_close(finder.get_layersG(2), other, 1e-9);
		other += finder.get_layers(2);
		images_are_close(finder.get_layersG(3), other, 1e-9);
		other += finder.get_layers(3);
		images_are_close(finder.get_layersG(4), other, 1e-9);
	}
	BOOST_AUTO_TEST_CASE( fill_rectangular )
	{
		OctaveFinder finder(51, 103);
		cv::Mat_<double>input(51, 103), other(51, 103);
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
		cv::Mat_<double>input(256, 256);
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
		cv::Mat_<double>input(256, 256);
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
		cv::Mat_<double>input(128, 512);
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
		cv::Mat_<double>input(151, 103), other(151, 103);
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
		cv::Mat_<double>input(256, 256);
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
		std::vector<cv::Vec4d> v = finder.subpix();
	}

	BOOST_AUTO_TEST_CASE( multiple_circles_various_sizes )
	{
		OctaveFinder finder;
		cv::Mat_<double>input(256, 256);
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
		cv::Mat_<double>input(256, 256);
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
		cv::Mat_<double>input(256, 256);
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
		cv::Mat_<double>input(256, 256);
		input.setTo(0);
		cv::circle(input, cv::Point(100, 200), 4, 1.0, -1);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		BOOST_CHECK_EQUAL(finder.get_binary(2)(200, 100), 1);
		std::vector<cv::Vec4d> v = finder.subpix();
		BOOST_REQUIRE_EQUAL(v.size(), 1);
		//x
		BOOST_CHECK_GE(v[0][0], 99);
		BOOST_CHECK_LE(v[0][0], 101);
		//y
		BOOST_CHECK_GE(v[0][1], 199);
		BOOST_CHECK_LE(v[0][1], 201);
		//scale
		BOOST_CHECK_GE(v[0][2], 1);
		BOOST_CHECK_LE(v[0][2], 3);
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
		cv::Mat_<double>input(503, 151);
		input.setTo(0);
		cv::circle(input, cv::Point(100, 200), 4, 1.0, -1);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		BOOST_CHECK_EQUAL(finder.get_binary(2)(200, 100), 1);
		std::vector<cv::Vec4d> v = finder.subpix();
		BOOST_REQUIRE_EQUAL(v.size(), 1);
		//x
		BOOST_CHECK_GE(v[0][0], 99);
		BOOST_CHECK_LE(v[0][0], 101);
		//y
		BOOST_CHECK_GE(v[0][1], 199);
		BOOST_CHECK_LE(v[0][1], 201);
		//scale
		BOOST_CHECK_GE(v[0][2], 1);
		BOOST_CHECK_LE(v[0][2], 3);

	}

	BOOST_AUTO_TEST_CASE( subpix_relative_positions )
	{
		OctaveFinder finder(32, 32);
		//cv cannot draw circle positions better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(8*32, 8*32), small_input(32,32);
		std::vector<cv::Vec4d> v;
		for(int i=0; i<7; ++i)
		{
			input.setTo(0);
			cv::circle(input, cv::Point(8*(16+pow(0.5, i+1)), 8*16), 8*4, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			finder.preblur_and_fill(small_input);
			finder.initialize_binary();
			std::vector<cv::Vec4d> v_s = finder.subpix();
			std::copy(v_s.begin(), v_s.end(), std::back_inserter(v));
		}
		BOOST_REQUIRE_EQUAL(v.size(), 7);
		for(int i=0; i<6; ++i)
			BOOST_WARN_MESSAGE(
					v[i][0] > v[i+1][0],
					"spatial resolution is larger than 1/" << pow(2,i+1) <<"th of a pixel"
					);
		BOOST_CHECK_GT(v[2][0], v[3][0]);
		//can differentiate between 1/8 and 1/16, but not between 1/16 and 1/32

	}

	BOOST_AUTO_TEST_CASE( subpix_relative_sizes )
	{
		OctaveFinder finder(32,32);
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(8*32, 8*32), small_input(32,32);
		std::vector<cv::Vec4d> v;
		for(int i=0; i<7; ++i)
		{
			input.setTo(0);
			cv::circle(input, cv::Point(8*16, 8*16), 8*(4+0.125*i), 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			finder.preblur_and_fill(small_input);
			finder.initialize_binary();
			std::vector<cv::Vec4d> v_s = finder.subpix();
			std::copy(v_s.begin(), v_s.end(), std::back_inserter(v));
		}
		BOOST_REQUIRE_EQUAL(v.size(), 7);

		for (int i=1; i<7;++i)
			BOOST_WARN_MESSAGE(v[0][2]< v[7-i][2], "resolution in size is 1/"<<(8*(7-i))<< "th of a scale");

		finder.scale(v);
		BOOST_CHECK_CLOSE(v[0][2], 4, 2);
		BOOST_CHECK_CLOSE(v[1][2], 4.125, 2);
		BOOST_CHECK_CLOSE(v[2][2], 4.25, 2);
		BOOST_CHECK_CLOSE(v[3][2], 4.325, 2);
		BOOST_CHECK_CLOSE(v[4][2], 4.5, 2);
		BOOST_CHECK_CLOSE(v[5][2], 4.625, 2);
		BOOST_CHECK_CLOSE(v[6][2], 4.75, 2);
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
		cv::Mat_<double>input(256, 256);
		input.setTo(0);
		cv::circle(input, cv::Point(128, 128), 4, 1.0, -1);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		boost::progress_timer ti;
		for (size_t i=0; i<100; ++i)
			std::vector<cv::Vec4d> v = finder.subpix();
		std::cout<<"100 subpixel resolutions in ";
	}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( octave_limit_cases )

	BOOST_AUTO_TEST_CASE( octave_minimum_detected_size )
	{
		OctaveFinder finder(32,32);
		cv::Mat_<uchar>input(16*32, 16*32), small_input(32,32);
		int i = 0;
		std::vector<cv::Vec4d> v(1);
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
		int i = 32;
		std::vector<cv::Vec4d> v(1);
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
			std::vector<cv::Vec4d> v_s = finder(small_input, true);
			BOOST_CHECK_MESSAGE(
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
			BOOST_CHECK_CLOSE(v_s[0][2], 4, 50);
			for(size_t j=0; j<v_s.size(); ++j)
				f << position << "\t" << v_s[j][1] << "\t"  << v_s[j][2] << "\n";
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
			std::vector<cv::Vec4d> v_s = finder(small_input, true);
			BOOST_REQUIRE_EQUAL(v_s.size(), 2);
			BOOST_CHECK_CLOSE(v_s[0][1], 16, distance<16?6:2);
			/*cv::namedWindow("truc");
			cv::imshow("truc", small_input);
			cv::imshow("binary", 255*finder.get_binary(2));
			cv::waitKey();*/
			BOOST_REQUIRE_CLOSE(v_s[1][1], 16+distance, distance<16?6:2);
			BOOST_CHECK_CLOSE(v_s[0][2], 4, distance<16?6:2);
			BOOST_CHECK_CLOSE(v_s[1][1] - v_s[0][1], distance, distance<9.03125?21:2);
			f << distance << "\t" << v_s[0][1] << "\t" << v_s[1][1] << "\t" << v_s[0][2] << "\n";
		}
		f<<std::endl;
	}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END() //octave

BOOST_AUTO_TEST_SUITE( multiscale )

BOOST_AUTO_TEST_SUITE( multiscale_constructors )

	BOOST_AUTO_TEST_CASE( multiscale_constructor_square )
	{
		MultiscaleFinder finder;
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
		MultiscaleFinder finder(101, 155, 1, 1.0);
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
		MultiscaleFinder finder;
		cv::Mat_<double>input(256, 256);
		input.setTo(0);
		cv::circle(input, cv::Point(128, 128), 4, 1.0, -1);
		std::vector<cv::Vec4d> v = finder(input);
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
		BOOST_CHECK_CLOSE(v[0][2], 4, 2);
	}
	BOOST_AUTO_TEST_CASE( multiscale_rectangular )
	{
		MultiscaleFinder finder(101, 155);
		BOOST_REQUIRE_EQUAL(finder.get_n_octaves(), 5);
		cv::Mat_<double>input(101, 155);
		input.setTo(0);
		cv::circle(input, cv::Point(28, 28), 2, 1.0, -1);
		std::vector<cv::Vec4d> v = finder(input);
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
		BOOST_CHECK_CLOSE(v[0][2], 2, 5);
	}
	BOOST_AUTO_TEST_CASE( multiscale_relative_sizes )
	{
		const int s = 5;
		MultiscaleFinder finder(pow(2, s+1), pow(2, s+1));
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(32*pow(2, s+1), 32*pow(2, s+1)), small_input(pow(2, s+1), pow(2, s+1));
		std::vector<cv::Vec4d> v;
		std::ofstream f("multiscale_relative_sizes.out");
		for(int k=1; k<s; ++k)
			for(int i=0; i<32; ++i)
			{
				input.setTo(0);
				const int large_radius = 32*1.5*pow(2, k-1)+i*pow(2, k);
				const double radius =  large_radius/32.0;
				BOOST_TEST_CHECKPOINT("radius = "<<radius<<" call");
				cv::circle(input, cv::Point(32*pow(2, s), 32*pow(2, s)), large_radius, 255, -1);
				cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
				std::vector<cv::Vec4d> v_s = finder(small_input);
				BOOST_TEST_CHECKPOINT("radius = "<<radius<<" size");
				BOOST_CHECK_MESSAGE(
						v_s.size()==1,
						""<<((v_s.size()==0)?"No center":"More than one center")<<" for input radius "<<radius
						);
				/*if(radius<pow(2, s-1))
					BOOST_CHECK_CLOSE(v_s[0][2], radius, 5);*/
				BOOST_TEST_CHECKPOINT("radius = "<<radius<<" copy");
				std::copy(v_s.begin(), v_s.end(), std::back_inserter(v));
				for(size_t j=0; j<v_s.size(); ++j)
					f << radius << "\t" << v_s[j][2] << "\n";
			}
		f<<std::endl;
		BOOST_REQUIRE_EQUAL(v.size(), 32*(s-1));
	}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END() //multiscale
