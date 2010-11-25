/**
    Copyright 2008,2009 Mathieu Leocmach

    This file is part of Colloids.

    Colloids is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Colloids is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Colloids.  If not, see <http://www.gnu.org/licenses/>.
**/

#include "particles.hpp"
#include "files_series.hpp"
#include <boost/progress.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace std;
using namespace Colloids;

int main(int argc, char ** argv)
{
    try
    {
        string filename, token, outputPostfix;
        double minSep;
        size_t t_offset;
        po::options_description
            compulsory_options("Compulsory options"),
            additional_options("Additional options"),
            cmdline_options("Command-line options");
        compulsory_options.add_options()
            ("input", po::value<string>(&filename), "input file or pattern")
            ("minSep", po::value<double>(&minSep), "Minimum separation between particles (in pixels)")
            ;
        po::positional_options_description pd;
        pd.add("input", 1);
        pd.add("minSep", 2);
        additional_options.add_options()
            ("help", "produce help message")
            ("token", po::value<string>(&token)->default_value("_t"), "Token delimiting time step number (for time series only)")
            ("span", po::value<size_t>(), "Number of time steps to process (compulsory for time series)")
            ("offset", po::value<size_t>(&t_offset)->default_value(0), "Starting time step")
            ("both", "Removes both particles if they are closer than minSep (removes only the second one by default, ie the dimmer one)")
            ("outputPostfix,o", po::value<string>(&outputPostfix)->default_value(""), "If empty output prefix is given, overwrites the files.")
            ;

        cmdline_options.add(compulsory_options).add(additional_options);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(cmdline_options).positional(pd).run(), vm);

        if (vm.count("help") || (!vm.count("input") && !vm.count("minSep")))
        {
            cout << "cutter input minSep [options]\n";
            cout << cmdline_options << "\n";
            return EXIT_SUCCESS;
        }
        po::notify(vm);

        const string inputPath = filename.substr(0,filename.find_last_of("."));
        const string ext = filename.substr(filename.find_last_of(".")+1);

        if(vm.count("span"))
        {
            cout<<"file serie"<<endl;
            FileSerie datSerie(filename, token, vm["span"].as<size_t>(), t_offset),
                outSerie = datSerie.addPostfix(outputPostfix, ".dat");
            boost::progress_display show_progress(vm["span"].as<size_t>());
            if(vm.count("both"))
                for(size_t t=0; t<vm["span"].as<size_t>(); ++t)
                {
                    Particles parts(datSerie%t);
                    parts.makeRTreeIndex();
                    parts.removeShortRange(minSep).exportToFile(outSerie%t);
                    ++show_progress;
                }
            else
                for(size_t t=0; t<vm["span"].as<size_t>(); ++t)
                {
                    Particles(datSerie%t).cut(minSep).exportToFile(outSerie%t);
                    ++show_progress;
                }
        }
        else
            if(vm.count("both"))
            {
                Particles parts(filename);
                parts.makeRTreeIndex();
                parts.removeShortRange(minSep).exportToFile(inputPath+outputPostfix+".dat");
            }
            else
                Particles(filename).cut(minSep).exportToFile(inputPath+outputPostfix+".dat");


    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
