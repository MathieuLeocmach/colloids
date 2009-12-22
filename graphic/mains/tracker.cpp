/**
 * \file mainTracker.cpp
 * \brief main tracking program
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 7 april 2009
 *
 * Object oriented implementation of the code elaborated in 2008.
 * Old auto-configuration routine are commented at the end but not implemented this time
 * Based on Crocker & Grier algorithm
 *
 *
 * Program options are powered by boost::program_options (http://www.boost.org)
 * Need a link to boost_program_options-mt (libboost_program_options-mt.lib renamed to libboost_program_options-mt.a in MINGW)
 * A good tutorial (in French) for compiling BOOST under MinGW : http://devtricks.wordpress.com/installer-boost-sous-windows-avec-mingw/
 *
 * LIF files header parsing is powered by TinyXML (http://www.grinninglizard.com/tinyxml/)
 * You must define the preprocessor macro "TIXML_USE_STL"
 *
 * The use of FFTW, (http://www.fftw.org/)the quickest implementation of n-dimension FFT, is strongly recommended.
 * If the image sizes (x,y or z) are not a power of 2, it is even compulsory
 * Add the pre-processor definition "cimg_use_fftw3", link to fftw3-3 and put libfftw3-3.dll where the program can find it.
 *
 * To Enables CImg to use a native support for Tiff files, add the compilation definition "cimg_use_tiff"
 * File loading becomes ~2.5 time faster
 * You will need LibTIFF (http://www.libtiff.org) in your computer.
 * The library must be linked in the compiler configuration (include path + link to "tiff" in the gnuwin32/lib directory)
 * The "bin" path of libtiff (containing libtiff3.dll) must be in the windows PATH (reboot necessary)
 */

#include <boost/program_options.hpp>
#include "../graphicParticles.hpp"
#include "../tracker.hpp"
#include <boost/progress.hpp>

namespace po = boost::program_options;
using namespace std;

#ifndef INSTAL_PATH
#define INSTAL_PATH "c:/bin/"
#endif

int main(int ac, char* av[])
{
    try {
        string inputFile,outputPath;
        double Zratio=1.0,displayRadius,radiusMin, radiusMax, zradiusMin, zradiusMax,threshold;
        size_t serie;

        // Declare a group of options that will be
        // allowed only on command line
        po::options_description generic("Generic options");
        generic.add_options()
            ("help", "produce help message")
            ("modeLIF,L", "Read from a Leica LIF file")
            ("modeSerie,S", "Read from a 2D image files serie")
            ("view", "Display intermediate images")
            ("quiet","Display only a progression bar")
            ;

        // Declare a group of options that will be
        // allowed both on command line and in
        // config file
        po::options_description config("Tracking configuration");
        config.add_options()
            ("input,i", po::value<string>(&inputFile), "input file")
            ("outputPath,o", po::value<string>(&outputPath), "Path to output coordinate files")
            ("serie", po::value<size_t>(&serie), "Serie to track (modeLIF only, optional)")
            ("displayRadius",po::value<double>(&displayRadius)->default_value(5.0), "real size of a particle in pixels")
            ("radiusMin", po::value<double>(&radiusMin), "Minimum radius of the band pass filter in xy")
            ("radiusMax", po::value<double>(&radiusMax), "Maximum radius of the band pass filter in xy")
            ("zradiusMin", po::value<double>(&zradiusMin), "Minimum radius of the band pass filter in z")
            ("zradiusMax", po::value<double>(&zradiusMax), "Maximum radius of the band pass filter in z")
            ("threshold", po::value<double>(&threshold)->default_value(0.0), "Minimum intensity of a center after band passing (0,255)")
            ;

        // Options specific to file serie mode, accessible for both config file and command line
        po::options_description serieOp("Specific to file serie mode");
        serieOp.add_options()
            ("channel", po::value<size_t>()->default_value(0), "channel")
            ("xsize", po::value<size_t>()->default_value(255), "number of colums in one image")
            ("ysize", po::value<size_t>()->default_value(255), "number of rows in one image")
            ("zsize", po::value<size_t>()->default_value(255), "number of planes to track")
            ("zoffset", po::value<size_t>()->default_value(0), "first plane to track")
            ("tsize", po::value<size_t>()->default_value(1), "number of time steps to track")
            ("toffset", po::value<size_t>()->default_value(0), "first time step to track")
            ("voxelWidth", po::value<double>()->default_value(1.0), "real size of a pixel in a xy plane")
            ("voxelDepth", po::value<double>()->default_value(1.0), "real size of a pixel in z")
            ;


        po::options_description cmdline_options("Allowed options");
        cmdline_options.add(generic).add(config).add(serieOp);

        po::options_description config_file_options;
        config_file_options.add(config).add(serieOp);

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, cmdline_options), vm);

        if (vm.count("help") || !vm.count("modeLIF") && !vm.count("modeSerie"))
        {
            cout << cmdline_options << "\n";
            return EXIT_SUCCESS;
        }

        ifstream ifs(INSTAL_PATH "tracker.ini");
        if(!ifs)
        {
            cout << "Cannot open " INSTAL_PATH "tracker.ini" <<endl;
            return EXIT_FAILURE;
        }
        po::store(parse_config_file(ifs, config_file_options), vm);
        notify(vm);
        ifs.close();

        Tracker::loadWisdom(INSTAL_PATH "wisdom.fftw");

        boost::progress_display *show_progress;

        if (vm.count("modeSerie"))
        {
            Zratio = vm["voxelDepth"].as<double>() / vm["voxelWidth"].as<double>();
            //initialize the tracker
            vector<string> tokens(2);
            tokens[0]="_t";
            tokens[1]="_z";
            TokenTree input(tokens,inputFile);
            cout << input<<endl;
            FileSerieTracker track(input,Zratio,
                vm["xsize"].as<size_t>(),
                vm["ysize"].as<size_t>(),
                vm["zsize"].as<size_t>(),
                vm["channel"].as<size_t>()
                );
            Tracker::saveWisdom(INSTAL_PATH "wisdom.fftw");
            track.view = !!vm.count("view");
            track.quiet = !!vm.count("quiet");
            track.displayRadius = displayRadius;
            track.makeBandPassMask(radiusMin, radiusMax, zradiusMin, zradiusMax);

            const size_t digits = max(track.fileSerie->nbdigit,(size_t)1);
            ostringstream os;
            os<<outputPath<<"t%|0"<<digits<<"d|.dat";
            boost::format outputFileName(os.str());

            if(track.quiet) show_progress = new boost::progress_display(vm["tsize"].as<size_t>());
            for(size_t t=vm["toffset"].as<size_t>();t<vm["tsize"].as<size_t>()+vm["toffset"].as<size_t>();++t)
            {
                track.load3D(t);
                GraphicParticles Centers = track.trackXYZ(threshold);

                //export into a file
                Centers.exportToFile((outputFileName % t).str());
                if(track.quiet) ++(*show_progress);
            }
        }

        if(vm.count("modeLIF"))
        {
            auto_ptr<LifFile> lif(new LifFile(inputFile));
            //saving xml header
            FILE * header = fopen((inputFile.substr(0,inputFile.find_last_of("."))+"_header.xml").c_str(),"w");
            if(!header)
            {
                cerr<<"Create the output directory first !"<<endl;
                return EXIT_FAILURE;
            }
            lif->Header.Print(header);
            fclose(header);

            if(!vm.count("serie"))
                serie = lif->chooseSerie();

            LifTracker track(lif,serie,lif->chooseChannel(serie));
            Tracker::saveWisdom(INSTAL_PATH "wisdom.fftw");
            cout << "tracker ok"<<endl;
            track.view = !!vm.count("view");
            track.quiet = !!vm.count("quiet");
            track.displayRadius = displayRadius;
            track.makeBandPassMask(radiusMin, radiusMax, zradiusMin, zradiusMax);

            ostringstream nbframes;
            nbframes << track.lif->getNbFrames(serie);
            const size_t digits = max(nbframes.str().size(),(size_t)1);
            ostringstream os;
            os<<outputPath<<"t%|0"<<digits<<"d|.dat";
            boost::format outputFileName(os.str());

            if(track.quiet) show_progress = new boost::progress_display(track.lif->getNbFrames(serie));
            boost::progress_timer ptimer;
            const size_t tsize = track.lif->getNbFrames(serie);
            for(size_t t=0;t<tsize;++t)
            {
                track.load3D(t);
                GraphicParticles Centers = track.trackXYZ(threshold);

                Centers.exportToFile((outputFileName % t).str());
                if(track.quiet) ++(*show_progress);
            }
        }
        Tracker::saveWisdom(INSTAL_PATH "wisdom.fftw");
    }
    catch(exception& e)
    {
        cout << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
