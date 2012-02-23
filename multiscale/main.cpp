#include "src/multiscalefinder.hpp"
#include "src/lifFile.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>
#include <boost/program_options.hpp>
//#include "H5Cpp.h"

namespace po = boost::program_options;
using namespace Colloids;
using namespace boost::posix_time;

int main(int ac, char* av[]){
	try {
		std::string input,output;
		int ser = -1;
		// Declare a group of options that will be
		// allowed only on command line
		po::options_description cmdline_options("Generic options");
		cmdline_options.add_options()
			("input,i", po::value<std::string>(&input), "Leica file to read as input")
			("output,o", po::value<std::string>(&output)->default_value("./output"), "file to output to")
			("series,s", po::value<int>(&ser), "Dataset number")
			("start", po::value<int>()->default_value(0), "Starting time step. The output file numbering will also start at that time.")
			("removeOverlap", po::value<double>()->default_value(0.5), "When two centers are closer than the sum of their radius x parameter, remove the weaker")
			("Octave0",
					"Enables upsampling of the image to get two times smaller particles (between 4 and 8 pixels in diameter)\n"
					"NOT resilient to noise"
					)
			("incore",
					"Do not use memory mapped file to store the data. If you have enough RAM this will speed up the calculation. If you don't it will trigger swapping and global slow down or even crash.\n"
					"If the Octave0 option in ON, the 0th octave still uses memory mapped file.")
			("verbose,v", "Output debugging information")
			("help", "Show this help and exit")
			;
		//Input file as positional option
		po::positional_options_description p_o;
		p_o.add("input", -1);
		po::variables_map vm;
		po::store(po::command_line_parser(ac, av).
		          options(cmdline_options).positional(p_o).run(), vm);
		//po::store(po::parse_command_line(ac, av, cmdline_options), vm);
		po::notify(vm);

		if(!!vm.count("help"))
		{
			std::cout << cmdline_options << std::endl;
			return EXIT_SUCCESS;
		}
		if(!vm.count("input"))
		{
			std::cerr<<cmdline_options << std::endl;
			std::cerr<<"input file needed" << std::endl;
			return EXIT_FAILURE;
		}
		//choose series
		LifReader reader(vm["input"].as<std::string>());
		if(!vm.count("series") || !(!!vm.count("series") && reader.contains(vm["series"].as<int>())))
		{
			ser = reader.chooseSerieNumber();
		}
		LifSerie &serie = reader.getSerie(ser);
		std::cout<<"Tracking "<<serie.getName()<<std::endl;
		//create hdf5 file erase the content
		/*H5::H5File h5file(output, H5F_ACC_TRUNC);
		//Create a group that contains the sample
		H5::Group sample( h5file->createGroup( "/"+serie.getName() ));*/

		//get the spatial dimensions and switch on their number (dimensionality of the dataset)
		const std::vector<size_t> dims = serie.getSpatialDimensions();
		switch (dims.size())
		{
		case 3:
		{
			//initialize the finder
			int dimsint[3] = {dims[2], dims[1], dims[0]};
			MultiscaleFinder3D finder(dimsint[0], dimsint[1], dimsint[2], 3, 1.6, vm.count("incore"));
			//set the voxel size ratio (sampling in Z is often poorer than in X and Y)
			//finder.set_ZXratio(serie.getZXratio()); //DISABLED Particles get lost
			if(!vm.count("Octave0"))
				finder.disable_Octave0();
			//do not preblur in Z to compensate for the point spread function
			finder.set_Zpreblur(false);
			//re-open the LIF as a memory mapped file
			boost::iostreams::mapped_file_source file(vm["input"].as<std::string>());
			//container for tracked particles
			std::vector<Center3D> centers;
			ptime past_total = microsec_clock::local_time();
			boost::progress_timer ti;

			//enumerate time steps
			std::auto_ptr<boost::progress_display> progress;
			if(!vm.count("verbose"))
				progress.reset(new boost::progress_display(serie.getNbTimeSteps()));
			for(size_t t=vm["start"].as<int>(); t<serie.getNbTimeSteps(); ++t)
			{
				//create the input image header pointing to the right portion of the memory mapped file
				cv::Mat_<uchar> image = cv::Mat(3, dimsint, CV_8UC1, (unsigned char*)(file.data() + serie.getOffset(t)));
				//multiscale tracking
				if(!vm.count("verbose"))
					finder.get_centers(image, centers);
				else
				{
					std::cout<<"t = "<<t<<std::endl;
					{
						std::cout<<"fill ";
						ptime past = microsec_clock::local_time();
						boost::progress_timer ti;
						finder.fill(image);
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
				}
				//scale z according to the Z/X ratio of the image voxel
				const double ZXratio = serie.getZXratio();
				for(size_t c=0; c<centers.size(); ++c)
					centers[c][2] *=  ZXratio;
				//remove overlap
				removeOverlapping(centers, vm["removeOverlap"].as<double>());
				//output
				std::ostringstream os;
				os << output <<"_t"<< std::setfill('0') << std::setw(3) << t;
				if(!!vm.count("verbose"))
					std::cout << "output to "<< os.str() <<std::endl;
				std::ofstream out(os.str().c_str());
				for(size_t c=0; c<centers.size(); ++c)
				{
					for(int d=0; d<3; ++d)
						out << centers[c][d] << "\t";
					out << centers[c].r << "\t" << centers[c].intensity <<"\n";
				}
				out.close();
				if(progress.get())
					++(*progress.get());
			}
			std::cout<< "total time" << microsec_clock::local_time()-past_total <<" including CPU ";
		}
			break;
		default:
			throw std::invalid_argument("The dimensionality of the dataset is not yet supported");
		}

	}
	catch(std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
