#define BOOST_TEST_DYN_LINK

#include "../src/multiscalefinder.hpp"
#include "../src/traj.hpp"
#include "../src/reconstructor.hpp"
#include "../src/locatorfromlif.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>
#include <boost/array.hpp>
#include <fstream>
#include <set>
#include <limits>
#include <numeric>

using namespace Colloids;
using namespace boost::posix_time;

BOOST_AUTO_TEST_SUITE( LifTrack )
	/*BOOST_AUTO_TEST_CASE( fill_one_slice )
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
		MultiscaleFinder2D finder;
		cv::Mat_<unsigned char> slice(256, 256, (unsigned char)0);
		reader.getSerie(0).fill2DBuffer(static_cast<void*>(slice.data), 0, 114);
		std::vector<Center2D> centers = finder(slice);
		for(size_t l=1; l<finder.get_n_layers()+1; ++l)
		{
			std::ostringstream os;
			os << "fill_one_slice_bin_l" <<l<<".raw";
			std::ofstream out(os.str().c_str());
			out.write((char*)finder.get_octave(1).get_binary(l).data, 256*256);
		}
		for(size_t l=0; l<finder.get_n_layers()+3; ++l)
		{
			std::ostringstream os;
			os << "fill_one_slice_layerG_l" <<l<<".raw";
			std::ofstream out(os.str().c_str());
			out.write((char*)finder.get_octave(1).get_layersG(l).data, 256*256*finder.get_octave(1).get_layersG(l).elemSize());
		}
		BOOST_CHECK_EQUAL(finder.get_octave(1).get_binary(1)(110, 111), true);
	}*/
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
			os<<"test_output/fill_one_stack_clusters"<< count_file++ <<".vtk";
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
			f_cl<<"SCALARS response double\n"
					"LOOKUP_TABLE default\n";
			for(Clusters::const_iterator cl=cl_min;cl!=cl_max;++cl)
				for(Cluster::const_iterator p = cl->begin(); p!=cl->end(); ++p)
					f_cl<< p->intensity <<"\n";
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
		BOOST_CHECK_EQUAL(centers.size(), 7951);
		//export to VTK format for human-eye comparison
		std::ofstream f_rec("test_output/fill_one_stack_reconstructed.vtk");
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
