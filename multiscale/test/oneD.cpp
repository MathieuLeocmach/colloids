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
		std::vector<int> ci(1, 128);
		BOOST_CHECK_LT(finder.gaussianResponse(ci, 2.5), finder.gaussianResponse(ci, 2.0));
		BOOST_CHECK_LT(finder.gaussianResponse(ci, 3.5), finder.gaussianResponse(ci, 2.5));
		BOOST_CHECK_LT(finder.gaussianResponse(ci, 3.5), finder.gaussianResponse(ci, 3.0));
		BOOST_CHECK_GT(finder.gaussianResponse(ci, 3.5)-finder.gaussianResponse(ci, 2.5), finder.get_layers(3)(0, 128));
		BOOST_CHECK_GT(finder.gaussianResponse(ci, 4.5)-finder.gaussianResponse(ci, 3.5), finder.get_layers(3)(0, 128));
		//lower bound
		BOOST_CHECK_LT(finder.gaussianResponse(ci, 1.5), finder.gaussianResponse(ci, 1.0));
		BOOST_CHECK_LT(finder.gaussianResponse(ci, 0.5), finder.gaussianResponse(ci, 0));
		//further than the top layer
		BOOST_CHECK_LT(finder.gaussianResponse(ci, 10.5), finder.gaussianResponse(ci, 0));
		BOOST_CHECK_LT(finder.gaussianResponse(ci, 15.5), finder.gaussianResponse(ci, 5));
	}
	BOOST_AUTO_TEST_CASE( asymetrical )
	{
		OctaveFinder1D finder(32);
		OctaveFinder::Image input(1, 32);
		input.setTo(0);
		std::fill_n(&input(0, 10), 7, 1);
		BOOST_REQUIRE_EQUAL(finder.get_centers<2>(input).size(), 1);
		input.setTo(0);
		std::fill_n(&input(0, 10), 6, 1);
		BOOST_REQUIRE_EQUAL(finder.get_centers<2>(input).size(), 1);
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
		std::vector<int> ci(2, 3);
		ci[0] = 100;
		BOOST_CHECK_GE(finder.scale_subpix(ci),2);
		BOOST_CHECK_LE(finder.scale_subpix(ci),4);
		Center1D c;
		finder.single_subpix(ci, c);
		BOOST_CHECK_GE(c.r, 2);
		BOOST_CHECK_LE(c.r, 4);
		std::vector<Center1D> v;
		finder.subpix(v);
		BOOST_REQUIRE_EQUAL(v.size(), 1);
		//x
		BOOST_CHECK_GE(v[0][0], 99);
		BOOST_CHECK_LE(v[0][0], 101);
		//scale
		BOOST_CHECK_GE(v[0].r, 2);
		BOOST_CHECK_LE(v[0].r, 4);
		//gaussian response
		BOOST_CHECKPOINT("gaussian response");
		BOOST_CHECK_CLOSE(finder.gaussianResponse(ci, 2.0), finder.get_layersG(2)(0, 100), 1e-9);
		BOOST_CHECK_CLOSE(finder.gaussianResponse(ci, 3.0)-finder.gaussianResponse(ci, 2.0), finder.get_layers(2)(0, 100), 1e-9);
		BOOST_CHECK_LT(finder.gaussianResponse(ci, 3.0), finder.gaussianResponse(ci, 2.0));
		BOOST_CHECK_LT(finder.gaussianResponse(ci, 2.5), finder.gaussianResponse(ci, 2.0));
		BOOST_CHECK_LT(finder.gaussianResponse(ci, 3.5), finder.gaussianResponse(ci, 2.5));
		BOOST_CHECK_LT(finder.gaussianResponse(ci, 3.5), finder.gaussianResponse(ci, 3.0));
		BOOST_CHECK_GT(finder.gaussianResponse(ci, 3.5)-finder.gaussianResponse(ci, 2.5), finder.get_layers(3)(0, 100));
		BOOST_CHECK_GT(finder.gaussianResponse(ci, 4.5)-finder.gaussianResponse(ci, 3.5), finder.get_layers(3)(0, 100));
		BOOST_CHECKPOINT("sublayers 1");
		boost::array<double,8> sublayerG;
		for(size_t u = 0;u < sublayerG.size(); ++u)
			sublayerG[u] = finder.gaussianResponse(ci, 2 + 0.5*u);
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
			sublayerG[u] = finder.gaussianResponse(ci, z - 1 + 0.5*u);
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
		std::vector<Center1D> v;
		std::ofstream f("test_output/1d_relative_sizes.out");
		for(int i=0; i<32; ++i)
		{
			input.setTo(0);
			cv::circle(input, cv::Point(16*32, 0), 16*2+i, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			finder.preblur_and_fill(small_input);
			finder.initialize_binary();
			std::vector<Center1D> v_s;
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
		std::vector<Center1D> v(1);
		while(3-0.01*i>0 && v.size()>0)
		{
			i++;
			input.setTo(0);
			cv::circle(input, cv::Point(16*16, 0), 16*(3-0.01*i), 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			v =	finder.get_centers<1>(small_input, true);
		}
		BOOST_WARN_MESSAGE(false, "Cannot detect sizes below "<< (3-0.01*(i-1))<<" pixels in 1D");
	}

	BOOST_AUTO_TEST_CASE( octave_minimum_detector_size )
	{
		int i = 32;
		std::vector<Center1D> v(1);
		while(i>0 && v.size()>0)
		{
			i--;
			OctaveFinder1D finder(i);
			cv::Mat_<uchar>input(1, 16*i), small_input(1,i);
			input.setTo(0);
			cv::circle(input, cv::Point(8*i, 0), 16*2, 255, -1);
			cv::resize(input, small_input, small_input.size(), 0, 0, cv::INTER_AREA);
			v =	finder.get_centers<1>(small_input, true);
		}
		BOOST_WARN_MESSAGE(false, "An 1D octave detector smaller than "<< (i+1)<<" pixels cannot detect anything");
	}
	BOOST_AUTO_TEST_CASE( octave_minimum_signal_intensity )
	{
		OctaveFinder1D finder(16);
		OctaveFinder::Image  input(1,16);
		int i = 0;
		std::vector<Center1D> v(1);
		while(i<10 && v.size()>0)
		{
			i++;
			input.setTo(0);
			cv::circle(input, cv::Point(8, 0), 4, pow(0.5, i), -1);
			v =	finder.get_centers<1>(input, true);
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
		std::vector<Center2D> v;
		finder.get_centers(input, v);
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
		std::vector<Center2D> v, v_s;
		std::ofstream f("test_output/multiscale1D_relative_sizes.out");
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
			finder.get_centers(small_input, v_s);
			removeOverlapping(v_s);
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
						for(int i=0; i<finder.get_octave(o).get_height(); ++i)
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
		std::ofstream f("test_output/close_neighbours1D.out");
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
				std::ofstream out("test_output/close_neighbours1D.layersG");
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
	BOOST_AUTO_TEST_CASE( multiscale1D_real_data )
	{
		boost::array<float, 10> signal = {{3.1579 ,  3.57596,  3.59357,  3.73922,  3.72126,  3.66554, 3.67778,  3.57596,  2.9114 ,  3.207}};
		cv::Mat_<float> input_margin(1, signal.size()*3, 0.0f);
		std::copy(signal.begin(), signal.end(), input_margin.begin()+signal.size());
		MultiscaleFinder1D finder_margin(input_margin.cols);
		std::vector<Center2D> blobs_margin;
		finder_margin.get_centers(input_margin, blobs_margin);
		std::ofstream f_margin("test_output/multiscale1D_real_data_margin.layersG");
		BOOST_REQUIRE_GE(finder_margin.get_n_octaves(), 2);
		for(size_t l=0; l<finder_margin.get_n_layers()+2; ++l)
		{
			std::copy(
					finder_margin.get_octave(1).get_layersG(l).begin(),
					finder_margin.get_octave(1).get_layersG(l).end(),
					std::ostream_iterator<float>(f_margin, " "));
			f_margin<<"\n";
		}
		BOOST_REQUIRE_EQUAL(blobs_margin.size(), 1);
		BOOST_CHECK_CLOSE(blobs_margin[0][0], 14.63085423, 2);
	}
BOOST_AUTO_TEST_SUITE_END() //multiscale

BOOST_AUTO_TEST_SUITE_END() //1D
