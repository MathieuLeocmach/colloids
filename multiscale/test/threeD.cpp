#define BOOST_TEST_DYN_LINK

#include "../src/multiscalefinder.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>
#include <boost/array.hpp>
#include <set>

using namespace Colloids;
using namespace boost::posix_time;

void images_are_close(const cv::Mat &a, const cv::Mat &b, float precision=1e-5);

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

BOOST_AUTO_TEST_SUITE( threeD )

BOOST_AUTO_TEST_SUITE( circular )

	BOOST_AUTO_TEST_CASE( one_pixel )
	{
		int dims[3] = {32, 32, 32};
		OctaveFinder::Image input(3, dims, (OctaveFinder::PixelType)0);
		input(4, 16, 16) = 255;
		CircularZ4D circ(6, 32, 32);
		circ.loadplanes(&input(0,0,0), 1, -3, 6);
		BOOST_CHECK_CLOSE(circ.getG(1, 1, 16, 16), 255, 1e-4);
		BOOST_CHECK_CLOSE(circ.getDoG(0, 1, 16, 16), 255, 1e-4);
		BOOST_CHECK_CLOSE(circ.getDoG(1, 1, 16, 16), -255, 1e-4);
		//is the minimum where we put it ?
		int ml, mk, mj, mi;
		OctaveFinder::PixelType value;
		circ.blockmin(1, 16, 16, ml, mk, mj, mi, value);
		BOOST_CHECK_EQUAL(ml, 1);
		BOOST_CHECK_EQUAL(mk, 1);
		BOOST_CHECK_EQUAL(mj, 16);
		BOOST_CHECK_EQUAL(mi, 16);
		BOOST_CHECK_CLOSE(value, -255, 1e-4);
		BOOST_CHECK(circ.is_localmin(1, 16, 16, ml, mk, mj, mi, value));
		//focus goes on the next plane of blocks, the first minimum should be out of focus
		++circ;
		BOOST_CHECK_CLOSE(circ.getG(1, -1, 16, 16), 255, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, 1, 16, 16), 0, 1e-4);
		circ.blockmin(1, 16, 16, ml, mk, mj, mi, value);
		BOOST_CHECK_GT(value, -128);
		BOOST_CHECK(!circ.is_localmin(1, 16, 16, ml, mk, mj, mi, value));
		//load two more planes that have a bright pixel, they should be out of focus
		input(4, 16, 16) = 0;
		input(0, 16, 16) = 128;
		circ.loadplanes(&input(0,0,0), 1);
		BOOST_CHECK_CLOSE(circ.getG(1, -1, 16, 16), 255, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, 1, 16, 16), 0, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, 3, 16, 16), 128, 1e-4);
		//focus goes on the next plane of blocks, the second minimum should be in focus
		++circ;
		BOOST_CHECK_CLOSE(circ.getG(1, -3, 16, 16), 255, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, -1, 16, 16), 0, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, 1, 16, 16), 128, 1e-4);
		circ.blockmin(1, 16, 16, ml, mk, mj, mi, value);
		BOOST_CHECK_EQUAL(ml, 1);
		BOOST_CHECK_EQUAL(mk, 1);
		BOOST_CHECK_EQUAL(mi, 16);
		BOOST_CHECK_CLOSE(value, -128, 1e-4);
		BOOST_CHECK(circ.is_localmin(1, 16, 16, ml, mk, mj, mi, value));
		//load two more planes that have a out-of-block brighter pixel neighbour
		input(0, 15, 16) = 256;
		circ.loadplanes(&input(0,0,0), 1);
		BOOST_CHECK_CLOSE(circ.getG(1, -3, 16, 16), 255, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, -1, 16, 16), 0, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, 1, 16, 16), 128, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, 3, 16, 16), 128, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, 3, 15, 16), 256, 1e-4);
		//focus goes on the third block minimum, which is not a local minimum
		++circ;
		BOOST_CHECK_CLOSE(circ.getG(1, -3, 16, 16), 0, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, -1, 16, 16), 128, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, 1, 16, 16), 128, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, 1, 15, 16), 256, 1e-4);
		circ.blockmin(1, 16, 16, ml, mk, mj, mi, value);
		BOOST_CHECK_EQUAL(ml, 1);
		BOOST_CHECK_EQUAL(mk, 1);
		BOOST_CHECK_EQUAL(mi, 16);
		BOOST_CHECK_CLOSE(value, -128, 1e-4);
		BOOST_CHECK(!circ.is_localmin(1, 16, 16, ml, mk, mj, mi, value));
		circ.blockmin(1, 14, 16, ml, mk, mj, mi, value);
		BOOST_CHECK_EQUAL(ml, 1);
		BOOST_CHECK_EQUAL(mk, 1);
		BOOST_CHECK_EQUAL(mj, 15);
		BOOST_CHECK_EQUAL(mi, 16);
		BOOST_CHECK(circ.is_localmin(1, 14, 16, ml, mk, mj, mi, value));
		circ.blockmin(1, 15, 16, ml, mk, mj, mi, value);
		BOOST_CHECK_EQUAL(ml, 1);
		BOOST_CHECK_EQUAL(mk, 1);
		BOOST_CHECK_EQUAL(mj, 15);
		BOOST_CHECK_EQUAL(mi, 16);
		BOOST_CHECK(circ.is_localmin(1, 15, 16, ml, mk, mj, mi, value));
		//insert 2 new planes where the pixel is at odd position in z
		input(0, 16, 16) = 0;
		input(0, 15, 16) = 0;
		input(1, 16, 16) = 300;
		circ.loadplanes(&input(0,0,0), 1);
		++circ;
		BOOST_CHECK_CLOSE(circ.getG(1, -3, 16, 16), 128, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, -1, 16, 16), 128, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, -1, 15, 16), 256, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, 2, 16, 16), 300, 1e-4);
		++circ;
		BOOST_CHECK_CLOSE(circ.getG(1, -3, 16, 16), 128, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, -3, 15, 16), 256, 1e-4);
		BOOST_CHECK_CLOSE(circ.getG(1, 0, 16, 16), 300, 1e-4);
		circ.blockmin(1, 16, 16, ml, mk, mj, mi, value);
		BOOST_CHECK_EQUAL(ml, 1);
		BOOST_CHECK_EQUAL(mk, 0);
		BOOST_CHECK_EQUAL(mi, 16);
		BOOST_CHECK_CLOSE(value, -300, 1e-4);
		BOOST_CHECK(circ.is_localmin(1, 16, 16, ml, mk, mj, mi, value));
	}

BOOST_AUTO_TEST_SUITE_END() //circular

BOOST_AUTO_TEST_SUITE( octave3D )

BOOST_AUTO_TEST_SUITE( octave3D_constructors )

	BOOST_AUTO_TEST_CASE( octave3D_constructor_square )
    {
        OctaveFinder3D finder;
        BOOST_CHECK_EQUAL(finder.get_depth(), 256);
        BOOST_CHECK_EQUAL(finder.get_width(), 256);
        BOOST_CHECK_EQUAL(finder.get_height(), 256);
        BOOST_CHECK_EQUAL(finder.get_n_layers(), 3);

    }
	BOOST_AUTO_TEST_CASE( octave3D_constructor_rectangular_flat )
    {
        OctaveFinder3D finder(12, 128, 512, 1, 1.0);
        BOOST_CHECK_EQUAL(finder.get_depth(), 12);
        BOOST_CHECK_EQUAL(finder.get_width(), 128);
        BOOST_CHECK_EQUAL(finder.get_height(), 512);
        BOOST_CHECK_EQUAL(finder.get_n_layers(), 1);
    }
	BOOST_AUTO_TEST_CASE( octave3D_constructor_rectangular_ten_layers )
    {
		OctaveFinder3D finder(12, 128, 512, 10, 2.8);
        BOOST_CHECK_EQUAL(finder.get_depth(), 12);
		BOOST_CHECK_EQUAL(finder.get_width(), 128);
		BOOST_CHECK_EQUAL(finder.get_height(), 512);
		BOOST_CHECK_EQUAL(finder.get_n_layers(), 10);
	}
	BOOST_AUTO_TEST_CASE( octave3D_constructor_incore )
	{
		OctaveFinder3D finder(256, 256, 256, 3, 1.6, true);
		BOOST_CHECK_EQUAL(finder.get_depth(), 256);
		BOOST_CHECK_EQUAL(finder.get_width(), 256);
		BOOST_CHECK_EQUAL(finder.get_height(), 256);
		BOOST_CHECK_EQUAL(finder.get_n_layers(), 3);

	}

BOOST_AUTO_TEST_SUITE_END() //constructors

BOOST_AUTO_TEST_SUITE( octave3D_fill )

	BOOST_AUTO_TEST_CASE( fill3D_test )
	{
		OctaveFinder3D finder(64, 64, 64);
		int dims[3] = {64,64,64};
		OctaveFinder::Image input(3, dims);
		//the finder should contain a copy of the input data
		input.setTo(1);
		finder.fill(input);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], 64*64*64, 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(input)[0], 1e-4);
		images_are_close(finder.get_layersG(0), input, 1e-4);
		//the internal layersG[0] should not be the same object as the input, but a deep copy
		input.setTo(0);
		input.at<float>(32, 32, 32) = 1.0;
		BOOST_CHECK_NE(cv::sum(finder.get_layersG(0))[0], 1.0);
		BOOST_CHECK_NE(cv::sum(finder.get_layersG(0))[0], cv::sum(input)[0]);
		finder.fill(input);
		//Gaussian blur should be normalized
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(finder.get_layersG(1))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(1))[0], cv::sum(finder.get_layersG(2))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(2))[0], cv::sum(finder.get_layersG(3))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(3))[0], cv::sum(finder.get_layersG(4))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(finder.get_layersG(4))[0], 1e-5);
		//The blur should have acted in all directions
		BOOST_CHECK_GT(finder.get_layersG(1).at<float>(32, 32, 31), 0);
		BOOST_CHECK_GT(finder.get_layersG(1).at<float>(32, 31, 32), 0);
		BOOST_CHECK_GT(finder.get_layersG(1).at<float>(31, 32, 32), 0);
		BOOST_CHECK_CLOSE(finder.get_layersG(1).at<float>(31, 32, 32), finder.get_layersG(1).at<float>(32, 31, 32), 1e-5);
		BOOST_CHECK_CLOSE(finder.get_layersG(1).at<float>(31, 32, 32), finder.get_layersG(1).at<float>(32, 32, 31), 1e-5);
	}
	BOOST_AUTO_TEST_CASE( fill3D_rectangular )
	{
		OctaveFinder3D finder(60, 63, 65);
		int dims[3] = {60,63,65};
		OctaveFinder::Image input(3, dims);
		//the finder should contain a copy of the input data
		input.setTo(1);
		finder.fill(input);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], 60*63*65, 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(input)[0], 1e-4);
		images_are_close(finder.get_layersG(0), input, 1e-4);
		//the internal layersG[0] should not be the same object as the input, but a deep copy
		input.setTo(0);
		input.at<float>(32, 33, 31) = 1.0;
		BOOST_CHECK_NE(cv::sum(finder.get_layersG(0))[0], 1.0);
		BOOST_CHECK_NE(cv::sum(finder.get_layersG(0))[0], cv::sum(input)[0]);
		finder.fill(input);
		//Gaussian blur should be normalized
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(finder.get_layersG(1))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(1))[0], cv::sum(finder.get_layersG(2))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(2))[0], cv::sum(finder.get_layersG(3))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(3))[0], cv::sum(finder.get_layersG(4))[0], 1e-5);
		BOOST_CHECK_CLOSE(cv::sum(finder.get_layersG(0))[0], cv::sum(finder.get_layersG(4))[0], 1e-5);
		//The blur should have acted in all directions
		BOOST_CHECK_GT(finder.get_layersG(1).at<float>(32, 33, 30), 0);
		BOOST_CHECK_GT(finder.get_layersG(1).at<float>(32, 32, 31), 0);
		BOOST_CHECK_GT(finder.get_layersG(1).at<float>(31, 33, 31), 0);
		BOOST_CHECK_CLOSE(finder.get_layersG(1).at<float>(31, 33, 31), finder.get_layersG(1).at<float>(32, 32, 31), 1e-5);
		BOOST_CHECK_CLOSE(finder.get_layersG(1).at<float>(31, 33, 31), finder.get_layersG(1).at<float>(32, 33, 30), 1e-5);
	}
BOOST_AUTO_TEST_SUITE_END() //fill

BOOST_AUTO_TEST_SUITE( local_max )

	BOOST_AUTO_TEST_CASE( single_sphere )
	{
		OctaveFinder3D finder(64, 64, 64);
		int dims[3] = {64,64,64};
		OctaveFinder::Image input(3, dims);
		input.setTo(0);
		//draw a sphere
		drawsphere(input, 32, 32, 32, 4.0, (OctaveFinder::PixelType)1.0);
		finder.preblur_and_fill(input);
		finder.initialize_binary();
		//there should be only one center
		BOOST_REQUIRE_EQUAL(finder.get_nb_centers(), 1);
		std::vector<int> ci = finder.get_center_pixel(0);
		//with a radius of 4, the maximum should be in layer 1
		BOOST_CHECK_EQUAL(ci.back(), 1);
		//The minimum should be at the center of the circle
		BOOST_CHECK_EQUAL(ci[0], 32);
		BOOST_CHECK_EQUAL(ci[1], 32);
		BOOST_CHECK_EQUAL(ci[2], 32);
	}
BOOST_AUTO_TEST_SUITE_END() //local_max3D

BOOST_AUTO_TEST_SUITE( subpix )
	BOOST_AUTO_TEST_CASE( subpix_relative_sizes3D )
	{
		OctaveFinder3D finder(32, 32, 32);
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		int dims[3] = {32,32,32}, ldims[3] = {8*32, 8*32, 8*32};
		cv::Mat_<uchar>input(3, ldims, (unsigned char)0), small_input(3, dims, (unsigned char)0);
		std::vector<Center3D> v;
		for(int i=0; i<24; ++i)
		{
			input.setTo(0);
			drawsphere(input, 8*16, 8*16, 8*16, 8*4+i);
			volume_shrink(input, small_input, 8);
			finder.preblur_and_fill(small_input);
			finder.initialize_binary();
			std::vector<Center3D> v_s;
			finder.subpix(v_s);
			std::copy(v_s.begin(), v_s.end(), std::back_inserter(v));
		}
		BOOST_REQUIRE_EQUAL(v.size(), 24);

		for (int i=1; i<7;++i)
			BOOST_WARN_MESSAGE(v[0].r< v[7-i].r, "resolution in size is 1/"<<(8*(7-i))<< "th of a scale");

		for(size_t c=0; c<v.size(); ++c)
			finder.scale(v[c]);
		std::ofstream out("test_output/subpix_relative_sizes3D");
		for(size_t c=0; c<v.size(); ++c)
		{
			BOOST_CHECK_CLOSE(v[c].r, 4+0.125*c, 2);
			out<<4+0.125*c<<"\t"<<v[c].r<<"\n";
		}
	}
	BOOST_AUTO_TEST_CASE( subpix_positions )
	{
		OctaveFinder3D finder(32, 32, 32);
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		int dims[3] = {32,32,32}, ldims[3] = {8*32, 8*32, 8*32};
		cv::Mat_<uchar>input(3, ldims, (unsigned char)0), small_input(3, dims, (unsigned char)0);
		std::vector<Center3D> v;
		for(int x=-4; x<4; ++x)
		{
			for(int i=0; i<24; ++i)
			{
				input.setTo(0);
				drawsphere(input, 8*16, 8*16, 8*16+x, 8*4+i);
				volume_shrink(input, small_input, 8);
				finder.preblur_and_fill(small_input);
				finder.initialize_binary();
				std::vector<Center3D> v_s;
				finder.subpix(v_s);
				BOOST_CHECK_MESSAGE(v_s.size()==1, "x="<<x/8.<<" r="<< 4+0.125*i<<"\t"<<v_s.size()<<" centers");
				std::copy(v_s.begin(), v_s.end(), std::back_inserter(v));
			}
		}
		BOOST_REQUIRE_EQUAL(v.size(), 8*24);

		for(size_t c=0; c<v.size(); ++c)
			finder.scale(v[c]);
		std::ofstream out("test_output/subpix_positions3D");
		for(size_t c=0; c<24; ++c)
		{
			const double r = 4.0+0.125*c,
					s = log(r/finder.get_prefactor()/finder.get_radius_preblur())/log(2.0)*finder.get_n_layers()-1;
			out<<r<<"\t"<<s<<"\t";
			for(int x=0; x<8; ++x)
			{
				BOOST_CHECK_CLOSE(v[c+24*x].r, 4+0.125*c, 2);
				out<<v[c+24*x][0]<<"\t"<<v[c+24*x].r<<"\t";
			}
			out<<"\n";
		}
	}
BOOST_AUTO_TEST_SUITE_END() //subpix 3D

BOOST_AUTO_TEST_SUITE( octave_limit_cases )

	BOOST_AUTO_TEST_CASE( octave3D_minimum_detected_size )
	{
		OctaveFinder3D finder(32, 32, 32);
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		int dims[3] = {32,32,32}, ldims[3] = {16*32, 16*32, 16*32};
		cv::Mat_<uchar>input(3, ldims, (unsigned char)0), small_input(3, dims, (unsigned char)0);
		int i = 0;
		std::vector<Center3D> v(1);
		while(4-0.01*i>0 && v.size()>0)
		{
			i++;
			input.setTo(0);
			drawsphere(input, 16*16, 16*16, 16*16, 16*(4-0.01*i), (unsigned char)255);
			volume_shrink(input, small_input, 16);
			v =	finder.get_centers<3>(small_input, true);
		}
		BOOST_WARN_MESSAGE(false, "Cannot detect sizes below "<< (4-0.01*(i-1))<<" pixels in 3D");
	}

	BOOST_AUTO_TEST_CASE( octave3D_minimum_detector_size )
	{
		int i = 13;
		std::vector<Center3D> v(1);
		while(i>0 && v.size()>0)
		{
			i--;
			OctaveFinder3D finder(i,i,i);
			int dims[3] = {i, i, i};
			cv::Mat_<uchar> input(3, dims);
			input.setTo(0);
			drawsphere(input, i/2, i/2, i/2, 4, (unsigned char)255);
			v =	finder.get_centers<3>(input, true);
		}
		BOOST_WARN_MESSAGE(false, "An octave detector smaller than "<< (i+1)<<" pixels cannot detect anything in 3D");
	}
	/*BOOST_AUTO_TEST_CASE( size_at_border3D )
	{
		OctaveFinder3D finder(24,24,24);
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		int dims[3] = {24,24,24}, ldims[3] = {16*24, 16*24, 16*24};
		cv::Mat_<uchar>input(3, ldims, (unsigned char)0), small_input(3, dims, (unsigned char)0);
		std::ofstream f("test_output/pos_size_at_border3D.out");
		for(int i = 16*12; i>16*5; --i)
		{
			const double position = i/16.0;
			BOOST_TEST_CHECKPOINT("position = "<<position);
			input.setTo(0);
			drawsphere(input, 16*12, 16*12, i, 16*5, 255);
			volume_shrink(input, small_input, 16);
			std::vector<Center3D> v_s = finder.get_centers<3>(small_input, true);
			BOOST_REQUIRE_MESSAGE(
				v_s.size()==1,
				""<<((v_s.size()==0)?"No center":"More than one center")<<" for input position "<<position
				);
			std::vector<int> ci(3, 12);
			ci[0] = position;
			BOOST_CHECK_CLOSE(v_s[0].r, 5, 50);
			for(size_t j=0; j<v_s.size(); ++j)
				f << position << "\t" << v_s[j][0] << "\t"  << v_s[j].r << "\n";
		}
		f<<std::endl;
	}*/

	BOOST_AUTO_TEST_CASE( close_neighbours3D )
	{
		OctaveFinder3D finder(40,24,24);
		//cv cannot draw circle sizes better than a pixel, so the input image is drawn in high resolution
		int dims[3] = {40,24,24}, ldims[3] = {8*40, 8*24, 8*24};
		cv::Mat_<unsigned char>input(3, ldims, (unsigned char)0), small_input(3, dims, (unsigned char)0);
		std::ofstream f("test_output/close_neighbours3D.out");
		for(int i = 8*20; i>8*6; --i)
		{
			const double distance = i/8.0;
			BOOST_TEST_CHECKPOINT("distance = "<<distance);
			input.setTo(0);
			drawsphere(input, 8*14, 8*12, 8*12, 8*4, (unsigned char)255);
			drawsphere(input, 8*14+i, 8*12, 8*12, 8*4, (unsigned char)255);
			volume_shrink(input, small_input, 8);
			std::vector<Center3D> v_s = finder.get_centers<3>(small_input, true);
			BOOST_REQUIRE_EQUAL(v_s.size(), 2);
			BOOST_CHECK_CLOSE(v_s[0][2], 14, distance<16?6:2);
			//BOOST_CHECK_CLOSE(v_s[1][2], 12+distance, distance<16?6:2);
			BOOST_CHECK_CLOSE(v_s[0].r, 4, distance<16?6:2);
			BOOST_CHECK_CLOSE(v_s[1][2] - v_s[0][2], distance, distance<9.06375?22:2);
			f << distance << "\t" << v_s[0][2] << "\t" << v_s[1][2] << "\t" << v_s[0].r << "\t" << v_s[1].r << "\n";
		}
		f<<std::endl;
	}

BOOST_AUTO_TEST_SUITE_END() //octave limit cases 3D

BOOST_AUTO_TEST_SUITE_END() //octave 3D

BOOST_AUTO_TEST_SUITE( multiscale )

BOOST_AUTO_TEST_SUITE( multiscale3D_constructors )

	BOOST_AUTO_TEST_CASE( multiscale3D_constructor_square )
	{
		MultiscaleFinder3D finder;
		BOOST_REQUIRE_EQUAL(finder.get_n_octaves(), 6);
		BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(1)).get_depth(), 256);
		BOOST_CHECK_EQUAL(finder.get_octave(1).get_width(), 256);
		BOOST_CHECK_EQUAL(finder.get_octave(1).get_height(), 256);
		BOOST_CHECK_EQUAL(finder.get_octave(1).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(0)).get_depth(), 512);
		BOOST_CHECK_EQUAL(finder.get_octave(0).get_width(), 512);
		BOOST_CHECK_EQUAL(finder.get_octave(0).get_height(), 512);
		BOOST_CHECK_EQUAL(finder.get_octave(0).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(2)).get_depth(), 128);
		BOOST_CHECK_EQUAL(finder.get_octave(2).get_width(), 128);
		BOOST_CHECK_EQUAL(finder.get_octave(2).get_height(), 128);
		BOOST_CHECK_EQUAL(finder.get_octave(2).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(3)).get_depth(), 64);
		BOOST_CHECK_EQUAL(finder.get_octave(3).get_width(), 64);
		BOOST_CHECK_EQUAL(finder.get_octave(3).get_height(), 64);
		BOOST_CHECK_EQUAL(finder.get_octave(3).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(4)).get_depth(), 32);
		BOOST_CHECK_EQUAL(finder.get_octave(4).get_width(), 32);
		BOOST_CHECK_EQUAL(finder.get_octave(4).get_height(), 32);
		BOOST_CHECK_EQUAL(finder.get_octave(4).get_n_layers(), 3);
		BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(5)).get_depth(), 16);
		BOOST_CHECK_EQUAL(finder.get_octave(5).get_width(), 16);
		BOOST_CHECK_EQUAL(finder.get_octave(5).get_height(), 16);
		BOOST_CHECK_EQUAL(finder.get_octave(5).get_n_layers(), 3);
	}
	BOOST_AUTO_TEST_CASE( multiscale3D_constructor_rectangular_flat )
    {
		MultiscaleFinder3D finder(156, 101, 155, 1, 1.0);
        BOOST_REQUIRE_EQUAL(finder.get_n_octaves(), 5);
        BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(1)).get_depth(), 156);
        BOOST_CHECK_EQUAL(finder.get_octave(1).get_width(), 101);
        BOOST_CHECK_EQUAL(finder.get_octave(1).get_height(), 155);
        BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(0)).get_depth(), 312);
        BOOST_CHECK_EQUAL(finder.get_octave(0).get_width(), 202);
        BOOST_CHECK_EQUAL(finder.get_octave(0).get_height(), 310);
        BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(2)).get_depth(), 78);
        BOOST_CHECK_EQUAL(finder.get_octave(2).get_width(), 50);
        BOOST_CHECK_EQUAL(finder.get_octave(2).get_height(), 77);
        BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(3)).get_depth(), 39);
        BOOST_CHECK_EQUAL(finder.get_octave(3).get_width(), 25);
        BOOST_CHECK_EQUAL(finder.get_octave(3).get_height(), 38);
        BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(4)).get_depth(), 19);
        BOOST_CHECK_EQUAL(finder.get_octave(4).get_width(), 12);
        BOOST_CHECK_EQUAL(finder.get_octave(4).get_height(), 19);
        BOOST_CHECK_CLOSE(finder.get_radius_preblur(), 1.0, 1e-9);
        finder.set_radius_preblur(1.6);
        BOOST_CHECK_CLOSE(finder.get_radius_preblur(), 1.6, 1e-9);
    }
	BOOST_AUTO_TEST_CASE( multiscale3D_constructor_incore )
		{
			MultiscaleFinder3D finder(256,256,256,3,1.6,true);
			BOOST_REQUIRE_EQUAL(finder.get_n_octaves(), 6);
			BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(1)).get_depth(), 256);
			BOOST_CHECK_EQUAL(finder.get_octave(1).get_width(), 256);
			BOOST_CHECK_EQUAL(finder.get_octave(1).get_height(), 256);
			BOOST_CHECK_EQUAL(finder.get_octave(1).get_n_layers(), 3);
			BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(0)).get_depth(), 512);
			BOOST_CHECK_EQUAL(finder.get_octave(0).get_width(), 512);
			BOOST_CHECK_EQUAL(finder.get_octave(0).get_height(), 512);
			BOOST_CHECK_EQUAL(finder.get_octave(0).get_n_layers(), 3);
			BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(2)).get_depth(), 128);
			BOOST_CHECK_EQUAL(finder.get_octave(2).get_width(), 128);
			BOOST_CHECK_EQUAL(finder.get_octave(2).get_height(), 128);
			BOOST_CHECK_EQUAL(finder.get_octave(2).get_n_layers(), 3);
			BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(3)).get_depth(), 64);
			BOOST_CHECK_EQUAL(finder.get_octave(3).get_width(), 64);
			BOOST_CHECK_EQUAL(finder.get_octave(3).get_height(), 64);
			BOOST_CHECK_EQUAL(finder.get_octave(3).get_n_layers(), 3);
			BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(4)).get_depth(), 32);
			BOOST_CHECK_EQUAL(finder.get_octave(4).get_width(), 32);
			BOOST_CHECK_EQUAL(finder.get_octave(4).get_height(), 32);
			BOOST_CHECK_EQUAL(finder.get_octave(4).get_n_layers(), 3);
			BOOST_CHECK_EQUAL(dynamic_cast<const OctaveFinder3D&>(finder.get_octave(5)).get_depth(), 16);
			BOOST_CHECK_EQUAL(finder.get_octave(5).get_width(), 16);
			BOOST_CHECK_EQUAL(finder.get_octave(5).get_height(), 16);
			BOOST_CHECK_EQUAL(finder.get_octave(5).get_n_layers(), 3);
		}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( multiscale3D_call )

	BOOST_AUTO_TEST_CASE( single_sphere )
	{
		MultiscaleFinder3D finder(32, 32, 32);
		int dims[3] = {32,32,32};
		cv::Mat_<uchar>input(3, dims, (unsigned char)0);
		input.setTo(0);
		drawsphere(input, 16, 16, 16, 5, (unsigned char)255);
		std::ofstream out("test_output/multiscale_single_sphere.raw", std::ios_base::binary);
		out.write((char*)input.data, 32*32*32);
		std::vector<Center3D> v;
		finder.get_centers(input, v);
		std::ofstream out0("test_output/multiscale_single_sphere_o0lg0.raw", std::ios_base::binary);
		out0.write((char*)finder.get_octave(0).get_layersG(0).data, 64*64*64*sizeof(OctaveFinder::PixelType));
		std::ofstream out1("test_output/multiscale_single_sphere_o1lg0.raw", std::ios_base::binary);
		out1.write((char*)finder.get_octave(1).get_layersG(0).data, 32*32*32*sizeof(OctaveFinder::PixelType));
		//with a radius of 5, the maximum should be in layer 2 of octave 1 and nowhere else
		for(size_t o=0; o<finder.get_n_octaves(); ++o)
		{
			const size_t u = finder.get_octave(o).get_nb_centers();
			BOOST_CHECK_MESSAGE(u==(o==1), "Octave "<<o<<" has "<< u <<" center"<<((u>1)?"s":"")<<" in layers ");
			if(u!=(o==1))
				for(size_t i=0; i<u; ++i)
					BOOST_CHECK_MESSAGE(false, finder.get_octave(o).get_center_pixel(i).back());
		}
		std::ofstream f("test_output/multiscale_single_sphere.out");
		for(size_t c=0; c<v.size(); ++c)
			f<<v[c][0]<<"\t"<<v[c][1]<<"\t"<<v[c][2]<<"\t"<<v[c].r<<"\t"<<v[c].intensity<<"\n";
		BOOST_REQUIRE_EQUAL(v.size(), 1);
		BOOST_CHECK_CLOSE(v[0][0], 16, 10);
		BOOST_CHECK_CLOSE(v[0][1], 16, 10);
		BOOST_CHECK_CLOSE(v[0][2], 16, 10);
		BOOST_CHECK_CLOSE(v[0].r, 5, 2);
	}

	BOOST_AUTO_TEST_CASE( multiscale3D_minimum_detector_size )
	{
		int i = 13;
		std::vector<Center3D> v(1);
		while(i>0 && v.size()>0)
		{
			i--;
			MultiscaleFinder3D finder(i,i,i);
			int dims[3] = {i, i, i};
			cv::Mat_<uchar> input(3, dims);
			input.setTo(0);
			drawsphere(input, i/2, i/2, i/2, 3, (unsigned char)255);
			finder.get_centers(input, v);
			std::ofstream f("test_output/multiscale3D_minimum_detector_size.raw", std::ios_base::binary);
			f.write((const char*)finder.get_octave(0).get_layersG(0).data, i*i*i*8*sizeof(OctaveFinder::PixelType));
		}
		BOOST_WARN_MESSAGE(false, "An multiscale detector smaller than "<< (i+1)<<" pixels cannot detect anything in 3D");
	}
BOOST_AUTO_TEST_SUITE_END() //multiscale 3D call

BOOST_AUTO_TEST_SUITE( john )

	BOOST_AUTO_TEST_CASE( exact_mono )
	{
		//read simulation output and select 1/8th of the volume
		std::vector<Center3D> simulation;
		Center3D c(0.0, 5.0);
		std::ifstream sim("test_input/poly00_phi05.dat");
		BOOST_REQUIRE_MESSAGE(sim.good(), "could not find test_input/poly00_phi05.dat");
		size_t nb=0;
		sim >> nb;
		sim >> nb;
		double boxsize;
		sim >> boxsize;
		sim >> boxsize;
		sim >> boxsize;
		for(size_t p=0; p<nb; ++p)
		{
			for(size_t d=0; d<3; ++d)
			{
				sim >> c[d];
				//periodic boundary conditions
				if(c[d]<0)
					c[d] += boxsize;
				if(c[d] >= boxsize)
					c[d] -= boxsize;
			}
			//select 1/8th of the volume
			//if(std::max(c[0], std::max(c[1], c[2])) > boxsize*0.5)
					//continue;
			//rescale
			for(size_t d=0; d<3; ++d)
				c[d] = (c[d]+0.5) * 128.0 / (2.0 + boxsize);
			simulation.push_back(c);
		}
		sim.close();
		//create an empty, black input picture as big as the simulation box + margins
		const double radius = 0.5 * 128.0 / (2.0 + boxsize);
		int dims[3] = {128, 128, 128}, ldims[3] = {512, 512, 512};
		cv::Mat_<uchar> input(3, dims);
		{
			cv::Mat_<uchar> linput(3, ldims);
			linput.setTo(0);
			//paint spheres at the positions and sizes given by the simulation output
			for(size_t p=0; p<simulation.size(); ++p)
			{
				drawsphere(linput,
						4*simulation[p][2],
						4*simulation[p][1],
						4*simulation[p][0],
						4*radius, (unsigned char)255);
			}
			volume_shrink(linput, input, 4);
		}
		std::ofstream f("test_output/john_exact_mono.raw", std::ios_base::binary);
		f.write((const char*)input.data, dims[0]*dims[1]*dims[2]);
		//track in 3D
		MultiscaleFinder3D finder(128,128,128);
		std::vector<Center3D> centers;
		finder.get_centers(input, centers);
		//compare track output and simulation output
		BOOST_REQUIRE(!centers.empty());
		BOOST_CHECK_GT(centers.size(), simulation.size());
		removeOverlapping(centers);
		BOOST_CHECK_EQUAL(centers.size(), simulation.size());
		std::ofstream out("test_output/john_exact_mono.csv");
		out<<"x y z r i\n";
		double meanr = 0.0, minr = centers.front().r , maxr = minr;
		for(size_t p=0; p<centers.size(); ++p)
		{
			for(size_t d=0; d<3; ++d)
				out<< centers[p][d] - radius << " ";
			out<<centers[p].r<<" "<<centers[p].intensity<<"\n";
			meanr += centers[p].r;
			if(centers[p].r < minr)
				minr = centers[p].r;
			if(maxr < centers[p].r)
				maxr = centers[p].r;
		}
		meanr /= centers.size();
		BOOST_CHECK_CLOSE(meanr, radius, 0.3);
		BOOST_CHECK_GT(minr, radius*0.9);
		BOOST_CHECK_LT(maxr, radius*0.11);
	}
	BOOST_AUTO_TEST_CASE( exact_mono_dilute )
	{
		//read simulation output
		std::vector<Center3D> simulation;
		std::vector<Center3D> simulationL;
		Center3D c(0.0, 5.0), cL=c;
		std::ifstream sim("test_input/poly00_phi05.dat");
		BOOST_REQUIRE_MESSAGE(sim.good(), "could not find test_input/poly00_phi05.dat");
		size_t nb=0;
		sim >> nb;
		sim >> nb;
		double boxsize;
		sim >> boxsize;
		sim >> boxsize;
		sim >> boxsize;
		for(size_t p=0; p<nb; ++p)
		{
			for(size_t d=0; d<3; ++d)
			{
				sim >> c[d];
				//periodic boundary conditions
				if(c[d]<0)
					c[d] += boxsize;
				if(c[d] >= boxsize)
					c[d] -= boxsize;
			}
			//rescale
			for(size_t d=0; d<3; ++d)
			{
				cL[d] = (c[d]+0.5) * 256.0 / (2.0 + boxsize);
				c[d] = (c[d]+0.5) * 128.0 / (2.0 + boxsize);
			}
			simulation.push_back(c);
			simulationL.push_back(cL);
		}
		sim.close();
		const double radius = 0.5 * 128.0 / (2.0 + boxsize),
				radiusL = 0.5 * 256.0 / (2.0 + boxsize);
		//track each particle independently
		std::ofstream out("test_output/john_exact_mono_dilute.csv");
		std::ofstream outL("test_output/john_exact_mono_large_dilute.csv");
		out<<"x y z r i\n";
		outL<<"x y z r i\n";
		MultiscaleFinder3D finder(24,24,24);
		std::vector<Center3D> centers;
		double meanr = 0.0, minr = 2*radius , maxr = 0.5*radius;
		double meanrL = 0.0, minrL = 2*radiusL , maxrL = 0.5*radiusL;
		int dims[3] = {24, 24, 24}, ldims[3] = {24*4, 24*4, 24*4};
		cv::Mat_<uchar> input(3, dims), linput(3, ldims);
		for(size_t p=0; p<simulation.size(); ++p)
		{
			linput.setTo(0);
			//draw a sphere about at the center, but keeping subpixel dilute position
			drawsphere(linput,
				4*(simulation[p][2]-floor(simulation[p][2])+12),
				4*(simulation[p][1]-floor(simulation[p][1])+12),
				4*(simulation[p][0]-floor(simulation[p][0])+12),
				4*radius, (unsigned char)255);
			volume_shrink(linput, input, 4);
			finder.get_centers(input, centers);
			removeOverlapping(centers);
			BOOST_REQUIRE_EQUAL(centers.size(), 1);
			for(size_t d=0; d<3; ++d)
				out<< centers[0][d] - radius << " ";
			out<<centers.front().r<<" "<<centers[0].intensity<<"\n";
			meanr += centers[0].r;
			if(centers[0].r < minr)
				minr = centers[0].r;
			if(maxr < centers[0].r)
				maxr = centers[0].r;
			linput.setTo(0);
			//draw a sphere about at the center, but keeping subpixel dilute position
			drawsphere(linput,
				4*(simulationL[p][2]-floor(simulationL[p][2])+12),
				4*(simulationL[p][1]-floor(simulationL[p][1])+12),
				4*(simulationL[p][0]-floor(simulationL[p][0])+12),
				4*radiusL, (unsigned char)255);
			volume_shrink(linput, input, 4);
			finder.get_centers(input, centers);
			removeOverlapping(centers);
			BOOST_REQUIRE_EQUAL(centers.size(), 1);
			for(size_t d=0; d<3; ++d)
				outL<< centers[0][d] - radius << " ";
			outL<<centers.front().r<<" "<<centers[0].intensity<<"\n";
			meanrL += centers[0].r;
			if(centers[0].r < minrL)
				minrL = centers[0].r;
			if(maxrL < centers[0].r)
				maxrL = centers[0].r;
		}
		meanr /= simulation.size();
		BOOST_CHECK_CLOSE(meanr, radius, 0.3);
		BOOST_CHECK_GT(minr, radius*0.9);
		BOOST_CHECK_LT(maxr, radius*0.11);
		meanrL /= simulationL.size();
		BOOST_CHECK_CLOSE(meanrL, radiusL, 0.3);
		BOOST_CHECK_GT(minrL, radiusL*0.9);
		BOOST_CHECK_LT(maxrL, radiusL*0.11);
	}
	BOOST_AUTO_TEST_CASE( exact_mono_large )
	{
		//read simulation output
		std::vector<Center3D> simulation;
		Center3D c(0.0, 5.0);
		std::ifstream sim("test_input/poly00_phi05.dat");
		BOOST_REQUIRE_MESSAGE(sim.good(), "could not find test_input/poly00_phi05.dat");
		size_t nb=0;
		sim >> nb;
		sim >> nb;
		double boxsize;
		sim >> boxsize;
		sim >> boxsize;
		sim >> boxsize;
		for(size_t p=0; p<nb; ++p)
		{
			for(size_t d=0; d<3; ++d)
			{
				sim >> c[d];
				//periodic boundary conditions
				if(c[d]<0)
					c[d] += boxsize;
				if(c[d] >= boxsize)
					c[d] -= boxsize;
			}
			//select 1/8th of the volume
			//if(std::max(c[0], std::max(c[1], c[2])) > boxsize*0.5)
					//continue;
			//rescale
			for(size_t d=0; d<3; ++d)
				c[d] = (c[d]+0.5) * 192.0 / (2.0 + boxsize);
			simulation.push_back(c);
		}
		sim.close();
		//create an empty, black input picture as big as the simulation box + margins
		const double radius = 0.5 * 192.0 / (2.0 + boxsize);
		BOOST_WARN_MESSAGE(false, "radius is "<<radius);
		int dims[3] = {192, 192, 192}, ldims[3] = {768, 768, 768};
		cv::Mat_<uchar> input(3, dims);
		{
			cv::Mat_<uchar> linput(3, ldims);
			linput.setTo(0);
			//paint spheres at the positions and sizes given by the simulation output
			for(size_t p=0; p<simulation.size(); ++p)
			{
				drawsphere(linput,
						4*simulation[p][2],
						4*simulation[p][1],
						4*simulation[p][0],
						4*radius, (unsigned char)255);
			}
			volume_shrink(linput, input, 4);
		}
		//track in 3D
		MultiscaleFinder3D finder(192,192,192);
		std::vector<Center3D> centers;
		finder.get_centers(input, centers);
		//compare track output and simulation output
		BOOST_REQUIRE(!centers.empty());
		BOOST_CHECK_GT(centers.size(), simulation.size());
		removeOverlapping(centers);
		BOOST_CHECK_EQUAL(centers.size(), simulation.size());
		std::ofstream out("test_output/john_exact_mono_large.csv");
		out<<"x y z r i\n";
		double meanr = 0.0, minr = centers.front().r , maxr = minr;
		for(size_t p=0; p<centers.size(); ++p)
		{
			for(size_t d=0; d<3; ++d)
				out<< centers[p][d] - radius << " ";
			out<<centers[p].r<<" "<<centers[p].intensity<<"\n";
			meanr += centers[p].r;
			if(centers[p].r < minr)
				minr = centers[p].r;
			if(maxr < centers[p].r)
				maxr = centers[p].r;
		}
		meanr /= centers.size();
		BOOST_CHECK_CLOSE(meanr, radius, 0.3);
		BOOST_CHECK_GT(minr, radius*0.9);
		BOOST_CHECK_LT(maxr, radius*0.11);
	}

BOOST_AUTO_TEST_SUITE_END() //John

BOOST_AUTO_TEST_SUITE( real )

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
		//track in 3D
		MultiscaleFinder3D finder(31, 64, 64);
		std::vector<Center3D> centers;
		finder.get_centers(image, centers);
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
		std::ofstream out("test_output/gel.csv");
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
		//now add the ZX ratio
		finder.set_ZXratio(1.0437917621692017);
		finder.get_centers(image, centers);
		//check track output
		BOOST_REQUIRE(!centers.empty());
		contains = 0;
		for(size_t p=0; p<centers.size(); ++p)
		{
			if(centers[p][0]>30 && centers[p][0]<39 &&
					centers[p][1]>15 && centers[p][1]<25 &&
					centers[p][2]>11 && centers[p][2]<25)
				++contains;
		}
		BOOST_WARN_MESSAGE(contains>0, "Bridge particle not detected with ZX ratio");
		BOOST_CHECK_MESSAGE(contains<2, ""<<contains<<" bridge particles detected with ZX ratio");
	}

	BOOST_AUTO_TEST_CASE( Z_preblur )
	{
		int dims[3] = {30, 25, 25};
		cv::Mat_<uchar> image(3, dims);
		image.setTo(0);
		//read test data from disk
		std::ifstream imf("test_input/Z_preblur.raw");
		BOOST_REQUIRE_MESSAGE(imf.good(), "could not find test_input/Z_preblur.raw");
		imf.read((char*)image.data, 30*25*25);
		imf.close();
		//track in 3D
		MultiscaleFinder3D finder(30, 25, 25);
		std::vector<Center3D> centers;
		finder.get_centers(image, centers);
		//check track output
		BOOST_REQUIRE(!centers.empty());
		BOOST_CHECK_EQUAL(centers.size(), 3);
		finder.set_halfZpreblur(true);
		finder.get_centers(image, centers);
		//check track output
		BOOST_REQUIRE(!centers.empty());
		BOOST_CHECK_EQUAL(centers.size(), 3);
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
		MultiscaleFinder3D finder(55, 50, 50);
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
		BOOST_CHECK_EQUAL(tooclose, 1);
		//preblur only in XY, not in Z
		OctaveFinder3D fiXY(55, 50, 50);
		OctaveFinder::Image input;
		image.convertTo(input, input.type());
		inplace_blurXY(input, 1.6);
		fiXY.fill(input);
		fiXY.initialize_binary();
		centers.clear();
		fiXY.subpix(centers);
		BOOST_REQUIRE(!centers.empty());
		BOOST_CHECK_EQUAL(centers.size(), 4);
		for(size_t c=0; c< centers.size(); ++c)
			fiXY.scale(centers[c]);
		//two centers should be closer than the sum of their radii
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
	}

BOOST_AUTO_TEST_SUITE_END() //real

BOOST_AUTO_TEST_SUITE_END() //multiscale 3D

BOOST_AUTO_TEST_SUITE_END() //threeD
