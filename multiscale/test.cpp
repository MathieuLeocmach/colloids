#define BOOST_TEST_MODULE multiscale test
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
//#include <boost/lambda/lambda.hpp>

#include "octavefinder.hpp"

using namespace Colloids;
//using namespace boost::lambda;

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
			BOOST_CHECK_MESSAGE(*u<precision, "at x=" << i <<" y=" << j);
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

	}

	BOOST_AUTO_TEST_CASE( subpix_relative_positions )
	{
		OctaveFinder finder;
		//cv cannot draw circle positions better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(8*256, 8*256), small_input(256,256);
		input.setTo(0);
		for(int i=0; i<7; ++i)
			cv::circle(input, cv::Point(8*(128+pow(0.5, i+1)), 8*((i+1)*32)), 8*4, 255, -1);
		//reduce the resolution of the input
		cv::resize(input, small_input, small_input.size());
		finder.preblur_and_fill(small_input);
		finder.initialize_binary();
		std::vector<cv::Vec4d> v = finder.subpix();
		BOOST_REQUIRE_EQUAL(v.size(), 7);
		std::vector<cv::Vec4d> v_by_y = v;
		std::sort(v_by_y.begin(), v_by_y.end(), by_coordinate<cv::Vec4d>(1));
		BOOST_CHECK_GT(v_by_y[0][0], v_by_y[1][0]);
		BOOST_CHECK_GT(v_by_y[1][0], v_by_y[2][0]);
		BOOST_CHECK_GT(v_by_y[2][0], v_by_y[3][0]);
		//can differentiate between 1/8 and 1/16, but not between 1/16 and 1/32
		BOOST_WARN_GT(v_by_y[3][0], v_by_y[4][0]);
		BOOST_WARN_GT(v_by_y[4][0], v_by_y[5][0]);
		BOOST_WARN_GT(v_by_y[5][0], v_by_y[6][0]);

	}

	BOOST_AUTO_TEST_CASE( subpix_relative_sizes )
	{
		OctaveFinder finder;
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		cv::Mat_<uchar>input(8*256, 8*256), small_input(256,256);
		input.setTo(0);
		for(int i=0; i<7; ++i)
			cv::circle(input, cv::Point(8*128, 8*((i+1)*32)), 8*(4+0.125*i), 255, -1);
		//reduce the resolution of the input
		cv::resize(input, small_input, small_input.size());
		finder.preblur_and_fill(small_input);
		finder.initialize_binary();
		std::vector<cv::Vec4d> v = finder.subpix();
		BOOST_REQUIRE_EQUAL(v.size(), 7);
		std::vector<cv::Vec4d> v_by_y = v;
		std::sort(v_by_y.begin(), v_by_y.end(), by_coordinate<cv::Vec4d>(1));
		BOOST_CHECK_LT(v_by_y[0][2], v_by_y[6][2]);
		BOOST_CHECK_LT(v_by_y[0][2], v_by_y[5][2]);
		BOOST_CHECK_LT(v_by_y[0][2], v_by_y[4][2]);
		BOOST_CHECK_LT(v_by_y[0][2], v_by_y[3][2]);
		BOOST_CHECK_LT(v_by_y[0][2], v_by_y[2][2]);
		//resolution in size is 1/4 of a scale
		BOOST_WARN_LT(v_by_y[0][2], v_by_y[1][2]);

	}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
