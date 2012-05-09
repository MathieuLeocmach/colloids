#define BOOST_TEST_DYN_LINK

#include "../src/multiscalefinder.hpp"
#include "../src/reconstructor.hpp"
#include <boost/test/unit_test.hpp>

using namespace Colloids;

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
		rec.push_back(centers, 15);
		BOOST_REQUIRE_EQUAL(rec.size(), 2);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
		BOOST_CHECK_CLOSE(rec.get_clusters().front().front().r, 1.0, 1e-9);
		BOOST_CHECK_CLOSE(rec.get_clusters().front().back().r, 0.5, 1e-9);
	}
	BOOST_AUTO_TEST_CASE( real_data )
	{
		Reconstructor rec;
		Reconstructor::Frame centers(1);
		centers.back()[0] = 112.91;
		centers.back()[1] = 120.444;
		centers.back().r = 1.70649;
		rec.push_back(centers);
		centers.back()[0] = 112.358;
		centers.back()[1] = 120.663;
		centers.back().r = 1.71538;
		rec.push_back(centers);
		centers.back()[0] = 112.431;
		centers.back()[1] = 121.25;
		centers.back().r = 1.59051;
		rec.push_back(centers);
		centers.back()[0] = 112.502;
		centers.back()[1] = 120.567;
		centers.back().r = 1.58156;
		rec.push_back(centers);
		BOOST_REQUIRE_EQUAL(rec.size(), 4);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
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
	BOOST_AUTO_TEST_CASE( split_real )
	{
		Reconstructor rec;
		Reconstructor::Frame centers(1);
		centers.back()[0] = 109.843;
		centers.back()[1] = 127.541;
		centers.back().r = 2.83824;
		centers.back().intensity = -5.8995;
		rec.push_back(centers);
		centers.back()[0] = 110.112;
		centers.back()[1] = 127.957;
		centers.back().r = 2.85067;
		centers.back().intensity = -6.68224;
		rec.push_back(centers);
		centers.back()[0] = 109.975;
		centers.back()[1] = 127.733;
		centers.back().r = 2.73108;
		centers.back().intensity = -6.93194;
		rec.push_back(centers);
		centers.back()[0] = 110.168;
		centers.back()[1] = 127.88;
		centers.back().r = 2.54598;
		centers.back().intensity = -7.28323;
		rec.push_back(centers);
		centers.back()[0] = 110.117;
		centers.back()[1] = 127.984;
		centers.back().r = 2.63205;
		centers.back().intensity = -7.36591;
		rec.push_back(centers);
		centers.back()[0] = 110.781;
		centers.back()[1] = 127.935;
		centers.back().r = 2.78161;
		centers.back().intensity = -6.05116;
		rec.push_back(centers);
		centers.back()[0] = 110.714;
		centers.back()[1] = 127.818;
		centers.back().r = 2.63343;
		centers.back().intensity = -6.07131;
		rec.push_back(centers);
		centers.back()[0] = 111.183;
		centers.back()[1] = 127.836;
		centers.back().r = 3.11998;
		centers.back().intensity = -7.8848;
		rec.push_back(centers);
		centers.back()[0] = 111.281;
		centers.back()[1] = 127.663;
		centers.back().r = 3.44455;
		centers.back().intensity = -9.67106;
		rec.push_back(centers);
		centers.back()[0] = 111.304;
		centers.back()[1] = 127.587;
		centers.back().r = 3.58571;
		centers.back().intensity = -12.0548;
		rec.push_back(centers);
		centers.back()[0] = 111.714;
		centers.back()[1] = 127.506;
		centers.back().r = 3.66983;
		centers.back().intensity = -12.5797;
		rec.push_back(centers);
		centers.back()[0] = 111.663;
		centers.back()[1] = 127.511;
		centers.back().r = 3.72784;
		centers.back().intensity = -12.6989;
		rec.push_back(centers);
		centers.back()[0] = 111.488;
		centers.back()[1] = 127.865;
		centers.back().r = 3.54346;
		centers.back().intensity = -12.9572;
		rec.push_back(centers);
		centers.back()[0] = 111.433;
		centers.back()[1] = 127.599;
		centers.back().r = 3.31208;
		centers.back().intensity = -11.081;
		rec.push_back(centers);
		centers.back()[0] = 111.525;
		centers.back()[1] = 127.506;
		centers.back().r = 3.1168 ;
		centers.back().intensity = -9.94293;
		rec.push_back(centers);
		centers.back()[0] = 111.267;
		centers.back()[1] = 127.671;
		centers.back().r = 2.72046;
		centers.back().intensity = -6.34706;
		rec.push_back(centers);
		BOOST_REQUIRE_EQUAL(rec.size(), 16);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
		OctaveFinder::Image signal(1, 16*3);
		//const double minrad = std::min_element(rec.get_clusters()[0].begin(), rec.get_clusters()[0].end(), compare_radii<3>())->r;
		std::fill(&signal(0,0), &signal(0,16), 0.0);
		OctaveFinder::PixelType * s = &signal(0,16);
		for(std::list<Center3D>::const_iterator it = rec.get_clusters()[0].begin(); it!=rec.get_clusters()[0].end(); ++it)
			*s++ = it->intensity;
		std::fill(signal.begin()+32, signal.end(), 0.0);
		MultiscaleFinder1D finder(16*3);
		std::ofstream f0("test_output/split_real.layersG0");
		for(size_t l=0; l<finder.get_n_layers()+3; ++l)
		{
			std::copy(
					finder.get_octave(0).get_layersG(l).begin(),
					finder.get_octave(0).get_layersG(l).end(),
					std::ostream_iterator<float>(f0, "\t"));
			f0<<"\n";
		}
		std::ofstream f("test_output/split_real.layersG");
		for(size_t l=0; l<finder.get_n_layers()+3; ++l)
		{
			std::copy(
					finder.get_octave(1).get_layersG(l).begin(),
					finder.get_octave(1).get_layersG(l).end(),
					std::ostream_iterator<float>(f, "\t"));
			f<<"\n";
		}
		rec.split_clusters();
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 2);
	}
	BOOST_AUTO_TEST_CASE( real_long )
	{
		Reconstructor rec;
		Reconstructor::Frame frame(1);
		std::ifstream input("test_input/real.blob");
		while(input.good())
		{
			input >> frame.back()[0] >> frame.back()[1];
			input >> frame.back().r >> frame.back().r;
			input >> frame.back().intensity;
			if(input.good())
				rec.push_back(frame);
		}
		BOOST_REQUIRE_EQUAL(rec.size(), 74);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
		std::deque<Center3D> cen;
		rec.get_blobs(cen);
		BOOST_CHECK_EQUAL(cen.size(), 7);
		std::vector<Center3D> centers(cen.begin(), cen.end());
		removeOverlapping_brute_force<3>(centers);
		BOOST_CHECK_EQUAL(centers.size(), cen.size());
	}
	/*BOOST_AUTO_TEST_CASE( real_3particles )
	{
		Reconstructor rec;
		Reconstructor::Frame frame(1);
		std::ifstream input("test_input/real2.blob");
		while(input.good())
		{
			input >> frame.back()[0] >> frame.back()[1];
			input >> frame.back().r >> frame.back().r;
			input >> frame.back().intensity;
			if(input.good())
				rec.push_back(frame);
		}
		BOOST_REQUIRE_EQUAL(rec.size(), 21);
		BOOST_REQUIRE_EQUAL(rec.nb_cluster(), 1);
		std::deque<Center3D> cen;
		rec.get_blobs(cen);
		BOOST_CHECK_EQUAL(cen.size(), 3);
		std::vector<Center3D> centers(cen.begin(), cen.end());
		removeOverlapping_brute_force<3>(centers);
		BOOST_CHECK_EQUAL(centers.size(), cen.size());
	}
	BOOST_AUTO_TEST_CASE( one_sphere )
	{
		Reconstructor rec;
		std::ofstream f("test_output/one_sphere.out");
		for(int z0=0; z0<10; ++z0)
		{
			//slice a sphere of radius 4 centeres on 4+z0/10
			const double pos = 4.01 + z0 / 10.0;
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
	}*/
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
