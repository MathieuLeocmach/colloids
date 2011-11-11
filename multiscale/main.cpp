#include "src/multiscalefinder.hpp"
#include "src/lifFile.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>
#include <boost/program_options.hpp>

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
			("output,o", po::value<std::string>(&output)->default_value("./"), "folder to output to")
			("series,s", po::value<int>(&ser), "Dataset number")
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
		LifReader reader(vm["input"].as<std::string>());
		if(!vm.count("series") || !(!!vm.count("series") && reader.contains(vm["series"].as<int>())))
		{
			ser = reader.chooseSerieNumber();
		}
		LifSerie &serie = reader.getSerie(ser);
		std::cout<<"Tracking "<<serie.getName()<<std::endl;
	}
	catch(std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
