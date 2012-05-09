#define BOOST_TEST_DYN_LINK

#include "../src/locatorfromlif.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>
#include <boost/array.hpp>


using namespace Colloids;
using namespace boost::posix_time;

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
		cv::VideoWriter w("test_output/z_scan.avi", CV_FOURCC('D', 'I', 'V', 'X'), 16, slice.size(), true);
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
	BOOST_AUTO_TEST_CASE( locate_one_stack )
	{
		LifReader reader("/home/mathieu/Code_data/liftest/Tsuru11dm_phi=52.53_J36.lif");
		LifSerie &serie = reader.getSerie(0);
		boost::iostreams::mapped_file_source file("/home/mathieu/Code_data/liftest/Tsuru11dm_phi=52.53_J36.lif");
		const std::vector<size_t> dims = serie.getSpatialDimensions();
		int dimsint[3] = {dims[2], dims[1], dims[0]};
		cv::Mat_<uchar> input = cv::Mat(3, dimsint, CV_8UC1, (unsigned char*)(file.data() + serie.getOffset(0)));
		int sdims[3] = {dims[2], dims[1], dims[0]};
		int smin[3] = {0, 0, 0};
		cv::Mat_<uchar> small_input = input;/*(3, sdims);
		for(int k=0; k<sdims[0]; ++k)
			for(int j=0; j<sdims[1]; ++j)
			{
				const unsigned char * a = &input(k+smin[0], j+smin[1], smin[2]);
				unsigned char * b = &small_input(k, j, 0);
				for(int i=0; i<sdims[2]; ++i)
					*b++ = *a++;
			}*/
		std::ofstream f("test_output/locate_one_stack_small_input.raw", std::ios_base::binary);
		f.write((const char*)small_input.data, sdims[0]*sdims[1]*sdims[2]);
		MultiscaleFinder3D finder(sdims[0], sdims[1], sdims[2]);
		std::vector<Center3D> centers;
		{
			std::cout<<"fill ";
			ptime past = microsec_clock::local_time();
			boost::progress_timer ti;
			finder.fill(small_input);
			std::cout<< microsec_clock::local_time()-past <<" including CPU ";
		}
		{
			std::cout<<"initialize ";
			ptime past = microsec_clock::local_time();
			boost::progress_timer ti;
			finder.initialize_binary();
			std::cout<< microsec_clock::local_time()-past <<" including CPU ";
		}
		{
			std::cout<<"subpix ";
			ptime past = microsec_clock::local_time();
			boost::progress_timer ti;
			finder.subpix(centers);
			std::cout<< microsec_clock::local_time()-past <<" including CPU ";
		}
		BOOST_REQUIRE(!centers.empty());
		BOOST_CHECK_EQUAL(centers.size(), 13768);
		BOOST_CHECK_LT((*std::max_element(centers.begin(), centers.end(), compare_coord<0>()))[0], sdims[2]);
		BOOST_CHECK_GT((*std::min_element(centers.begin(), centers.end(), compare_coord<0>()))[0], 0);
		BOOST_CHECK_LT((*std::max_element(centers.begin(), centers.end(), compare_coord<1>()))[1], sdims[1]);
		BOOST_CHECK_GT((*std::min_element(centers.begin(), centers.end(), compare_coord<1>()))[1], 0);
		std::vector<Center3D>::const_iterator max_z = std::max_element(centers.begin(), centers.end(), compare_coord<2>()),
				min_z = std::min_element(centers.begin(), centers.end(), compare_coord<2>());
		BOOST_CHECK_LT((*max_z)[2], sdims[0]);
		/*std::vector<int> max_zi = finder.get_octave(0).get_center_pixel(max_z-centers.begin());
		std::copy(max_zi.begin(), max_zi.end(), std::ostream_iterator<int>(std::cout, "\t"));
		for(int u=-1;u<2;++u)
			std::cout<<finder.get_octave(0).get_layersG(max_zi.back())(max_zi[2]+u, max_zi[1], max_zi[0])<<"\t";
		std::cout<<std::endl;*/

		BOOST_CHECK_GT((*min_z)[2], 0);
		/*std::vector<int> min_zi = finder.get_octave(0).get_center_pixel(min_z-centers.begin());
		std::copy(min_zi.begin(), min_zi.end(), std::ostream_iterator<int>(std::cout, "\t"));
		for(int u=-1;u<2;++u)
			std::cout<<finder.get_octave(0).get_layersG(min_zi.back())(min_zi[2]+u, min_zi[1], min_zi[0])<<"\t";
		std::cout<<std::endl;*/
		{
			std::cout<<"remove overlap ";
			ptime past = microsec_clock::local_time();
			boost::progress_timer ti;
			removeOverlapping(centers);
			std::cout<< microsec_clock::local_time()-past <<" including CPU ";
		}
		BOOST_REQUIRE(!centers.empty());
		BOOST_CHECK_EQUAL(centers.size(), 13721);

		//export to VTK format for human-eye comparison
		std::ofstream f_rec("test_output/locate_one_stack_nooverlap.vtk");
		f_rec<<"# vtk DataFile Version 3.0\n"
				"fill_one_stack\n"
				"ASCII\n"
				"DATASET POLYDATA\n"
				"POINTS "<<centers.size()<<" double\n";
		for(std::vector<Center3D>::const_iterator p=centers.begin();p!=centers.end();++p)
		{
			for(size_t d=0;d<3;++d)
				f_rec<<(*p)[d]<<" ";
			f_rec<<"\n";
		}
		f_rec<<"POINT_DATA "<<centers.size()<<"\n"
				"SCALARS r double\n"
				"LOOKUP_TABLE default\n";
		for(std::vector<Center3D>::const_iterator p=centers.begin();p!=centers.end();++p)
			f_rec<< p->r <<"\n";
		f_rec<<"SCALARS intensity double\n"
				"LOOKUP_TABLE default\n";
		for(std::vector<Center3D>::const_iterator p=centers.begin();p!=centers.end();++p)
			f_rec<< p->intensity <<"\n";
		f_rec.close();
	}

BOOST_AUTO_TEST_SUITE_END()
