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

 * \file mainTracker.cpp
 * \brief main tracking program
 * \author Mathieu Leocmach
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
#include "lifTracker.hpp"
#include "serieTracker.hpp"
#include "radiiTracker.hpp"
#include <boost/progress.hpp>
#include <boost/format.hpp>

namespace po = boost::program_options;
using namespace std;
using namespace Colloids;

/** \brief test the existence of tracker.ini, create it if necessary  */
string testTrackerIni()
{
    const string dir = string(getenv("HOME"))+"/.colloids",
        file = dir + "/tracker.ini";

    ifstream ifs(file.c_str());
    if(!ifs.good())
    {
        cout << "create default "<<file <<endl;
        if(!system( ("MKDIR \""+dir+"\"").c_str() ))
            throw(invalid_argument(("Cannot create the directory "+dir+". Check permissions").c_str()
                ));
        ofstream ofs(file.c_str());
        if(!ofs.good())
            throw(invalid_argument(("Cannot create "+file+ "\n"+
                "Create the directory "+dir+" or check its permissions").c_str()
                ));

        ofs <<"# Tracking configuration"<<endl;
        ofs <<"input="<<getenv("HOME")<<"/Code_data/liftest/Tsuru11dm_phi=52.53_J36.lif"<<endl;
        ofs <<"outputPath = "<<getenv("HOME")<<"/Code_output/liftest/test_"<<endl;
        ofs <<"radiusMin = 3.48"<<endl;
        ofs <<"radiusMax = 64"<<endl;
        ofs <<"threshold = 90"<<endl;
        ofs <<""<<endl;
        ofs <<"# Specific to file serie mode"<<endl;
        ofs <<"channel = 1"<<endl;
        ofs <<"xsize = 256"<<endl;
        ofs <<"ysize = 256"<<endl;
        ofs <<"zsize = 128"<<endl;
        ofs <<"tsize = 1"<<endl;
        ofs <<"voxelWidth = 1"<<endl;
        ofs <<"voxelDepth = 1"<<endl;
    }
    return dir;
}

/** \brief Functor to save a serie of particles to file */
struct ParticlesExporter : public unary_function<const list<Particles>&, void>
{
    size_t t;
    ofstream * nbs;
    FileSerie *outputFileName;
    boost::progress_display *show_progress;

    /** \brief  Constructor. Need the size of the serie only to get the right number of digits  */
    ParticlesExporter(const string& outputPath, size_t size, bool quiet)
    {
        t=0;
        outputFileName = new FileSerie(FileSerie::get0th(outputPath, size)+".dat", "_t", size, 0);
        if(quiet)
            show_progress = new boost::progress_display(size);
        else
            show_progress = 0;
        nbs = new ofstream((outputPath+".nb").c_str(), ios::out | ios::trunc);
    }
    void close()
    {
        if(nbs) delete nbs;
        if(outputFileName) delete outputFileName;
        if(show_progress) delete show_progress;
    }

    /** \brief  Export the Particles object to the next file of the serie.
        For example outputDir/30s_t025.dat
    */
    void operator()(const list<Particles>& parts)
    {
        if(!show_progress) cout<<"export to "<<(*outputFileName % t)<<endl;
        parts.front().exportToFile(*outputFileName % (t++));
        (*nbs) << parts.front().size() << endl;
        if(show_progress)
            ++(*show_progress);
    }
};

/** \brief Functor to save a serie of radii to file */
struct RadiiExporter : public unary_function<const Particles&, void>
{
    size_t t;
    FileSerie *outputFileName;
    boost::progress_display *show_progress;

    /** \brief  Constructor. Need the size of the serie only to get the right number of digits  */
    RadiiExporter(const string& outputPath, size_t size, bool quiet)
    {
        t=0;
        outputFileName = new FileSerie(FileSerie::get0th(outputPath, size)+".radii", "_t", size, 0);
        if(quiet)
            show_progress = new boost::progress_display(size-1);
        else
            show_progress = 0;
    }

    /** \brief  Export the Particles object to the next file of the serie.
        For example outputDir/30s_t025.dat
    */
    void operator()(const vector<double>& radii)
    {
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
            ofstream outputFile((*outputFileName % (t++)).c_str());
            copy(radii.begin(), radii.end(), ostream_iterator<double>(outputFile, "\n"));
            }
            if(show_progress)
                ++(*show_progress);

        }
    }
};

int main(int ac, char* av[])
{
    try {
        string inputFile,outputPath;
        //double Zratio=1.0,displayRadius,radiusMin, radiusMax, zradiusMin, zradiusMax,threshold;
        double Zratio=1.0,radiusMin, radiusMax, threshold;
        size_t serie, cores,onlyTimeStep;
        vector<double> thresholds;

        // Declare a group of options that will be
        // allowed only on command line
        po::options_description generic("Generic options");
        generic.add_options()
            ("help", "produce help message")
            ("modeLIF,L", "Read from a Leica LIF file")
            ("modeSerie,S", "Read from a 2D image files serie")
            ("extractRadii,R", "Extract radii from image data and already tracked coordinates (from a previous run without the -R option).")
            ("view", "Display intermediate images")
            ("quiet", "Display only a progression bar")
#if (TRACKER_N_THREADS>1)
            ("cores", po::value<size_t>(&cores)->default_value(TRACKER_N_THREADS), "Number of Cores to use")
            ("onlyTimeStep", po::value<size_t>(&onlyTimeStep)->default_value(-1), "Track only the given time step (default: track all time steps)")
#endif
            ;

        // Declare a group of options that will be
        // allowed both on command line and in
        // config file
        po::options_description config("Tracking configuration");
        config.add_options()
            ("input,i", po::value<string>(&inputFile), "input file")
            ("outputPath,o", po::value<string>(&outputPath), "Path to output coordinate files")
            ("serie", po::value<size_t>(&serie), "Serie to track (modeLIF only, optional)")
            //("displayRadius",po::value<double>(&displayRadius)->default_value(5.0), "real size of a particle in pixels")
            ("radiusMin", po::value<double>(&radiusMin), "Minimum radius of the band pass filter in xy")
            ("radiusMax", po::value<double>(&radiusMax), "Maximum radius of the band pass filter in xy")
            //("zradiusMin", po::value<double>(&zradiusMin), "Minimum radius of the band pass filter in z")
            //("zradiusMax", po::value<double>(&zradiusMax), "Maximum radius of the band pass filter in z")
            ("threshold", po::value<double>(&threshold)->default_value(0.0), "Minimum intensity of a center after band passing (0,255)")
            ("thresholds", po::value< vector<double> >(&thresholds), "Minimum intensity of a center after band passing (0,255)")
            ;

        // Options specific to file serie mode, accessible for both config file and command line
        po::options_description serieOp("Specific to file serie mode");
        serieOp.add_options()
            ("channel", po::value<size_t>()->default_value(0), "channel")
            ("xsize", po::value<size_t>()->default_value(255), "number of colums in one image")
            ("ysize", po::value<size_t>()->default_value(255), "number of rows in one image")
            ("zsize", po::value<size_t>()->default_value(255), "number of planes to track")
            //("zoffset", po::value<size_t>()->default_value(0), "first plane to track")
            ("tsize", po::value<size_t>()->default_value(1), "number of time steps to track")
            //("toffset", po::value<size_t>()->default_value(0), "first time step to track")
            ("voxelWidth", po::value<double>()->default_value(1.0), "real size of a pixel in a xy plane")
            ("voxelDepth", po::value<double>()->default_value(1.0), "real size of a pixel in z")
            ;


        po::options_description cmdline_options("Allowed options");
        cmdline_options.add(generic).add(config).add(serieOp);

        po::options_description config_file_options;
        config_file_options.add(config).add(serieOp);

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, cmdline_options), vm);

        if (vm.count("help") || (!vm.count("modeLIF") && !vm.count("modeSerie")))
        {
            cout << cmdline_options << "\n";
            return EXIT_SUCCESS;
        }

        const string config_path = testTrackerIni();
        ifstream ifs((config_path+"/tracker.ini").c_str());
        po::store(parse_config_file(ifs, config_file_options), vm);
        notify(vm);
        ifs.close();

        if(thresholds.empty())
            thresholds.push_back(threshold);

        //makes sure the last character of outputpath is not an underscore, because we will add one anyway
        //because the last 't' of "seriet000.dat" is not equivalent to the last '_t' of "serie_t000.dat"
        //so '_t' is a much better token than 't'
        if(outputPath.substr(outputPath.size()-1) == "_")
            outputPath.erase(outputPath.end()-1);

        //test writting on the given outputpath
        {
            ofstream output_test((outputPath+"test.test").c_str(), ios::out | ios::trunc);
            if(output_test.good())
                remove((outputPath+"test.test").c_str());
            else
                throw invalid_argument("Cannot write on specified output path :\t"+outputPath);
        }

        //configure FTTW to use multiple CPU
#if (TRACKER_N_THREADS>1)
        fftwf_init_threads();
        fftwf_plan_with_nthreads(cores);
#endif

        Tracker::loadWisdom(config_path+"/wisdom.fftw");


        if (vm.count("modeSerie"))
        {
            Zratio = vm["voxelDepth"].as<double>() / vm["voxelWidth"].as<double>();
            //initialize the tracker
            boost::array<size_t, 4> xyzt = {{
                vm["xsize"].as<size_t>(),
                vm["ysize"].as<size_t>(),
                vm["zsize"].as<size_t>(),
                vm["tsize"].as<size_t>()
            }};
            SerieTracker track(inputFile, xyzt, Zratio, vm["channel"].as<size_t>());
            cout << "tracker ok"<<endl;
            Tracker::saveWisdom(config_path+"/wisdom.fftw");
            track.setView(!!vm.count("view"));
            track.setQuiet(!!vm.count("quiet"));
            track.setThreshold(threshold);
            track.setIsotropicBandPass(radiusMin, radiusMax); //not using zradius for the moment

            boost::progress_timer ptimer;
            for_each(track, track.end(), ParticlesExporter(outputPath, vm["tsize"].as<size_t>(), track.quiet()));

        }

        if(vm.count("modeLIF"))
        {
            LifReader reader(inputFile);
            //saving xml header
            FILE * headerFile = fopen((inputFile.substr(0,inputFile.find_last_of("."))+"_header.xml").c_str(),"w");
            reader.getXMLHeader().Print(headerFile);
            fclose(headerFile);
            cout<<"header exported"<<endl;

            if(!vm.count("serie") || serie >= reader.getNbSeries())
                serie = reader.chooseSerieNumber();
            cout<<"You choose "<<reader.getSerie(serie).getName()<<endl;

            if(vm.count("extractRadii"))
            {
                cout<<"Extract radius"<<endl;
                LifRadiiTracker track(reader.getSerie(serie), outputPath);
                Tracker::saveWisdom(config_path+"/wisdom.fftw");
                cout << "tracker ok"<<endl;

                track.setView(!!vm.count("view"));
                track.setQuiet(!!vm.count("quiet"));

                boost::progress_timer ptimer;
                for_each(track, track.end(), RadiiExporter(outputPath, track.getLif().getNbTimeSteps(), track.quiet()));
                cout<<"radii extracted at all time steps in ";
            }
            else
            {
                cout<<"track particles"<<endl;
                LifTracker track(reader.getSerie(serie));
                cout << "tracker ok"<<endl;
                /*ofstream out((outputPath+"img.raw").c_str(), ios_base::out | ios_base::trunc);
                track.getTracker().copyImage(ostream_iterator<int>(out,"\t"));
                out.close();*/
                Tracker::saveWisdom(config_path+"/wisdom.fftw");
                track.setView(!!vm.count("view"));
                track.setQuiet(!!vm.count("quiet"));
                track.setThresholds(thresholds.begin(), thresholds.end());
                track.setIsotropicBandPass(radiusMin, radiusMax); //not using zradius for the moment
                //string maskoutput = outputPath+"mask.txt";
                //track.getTracker().maskToFile(maskoutput);

                if(onlyTimeStep==-1)
                {
                    boost::progress_timer ptimer;
                    ParticlesExporter pe(outputPath, track.getLif().getNbTimeSteps(), track.quiet());
                    for_each(track, track.end(), pe);
                    pe.close();
                }
                else
                {
                    track.setTimeStep(onlyTimeStep);
                    string nothrName = FileSerie
                    (
                        FileSerie::get0th(outputPath, track.getLif().getNbTimeSteps())+".dat",
                        "_t", track.getLif().getNbTimeSteps(), 0
                    ) % onlyTimeStep;
                    nothrName.insert(outputPath.size(), "_thr%g");
                    boost::format outputFileName(nothrName);

                    const list<Particles> & parts = *track;
                    vector<double>::const_iterator thr = thresholds.begin();
                    for(list<Particles>::const_iterator p=parts.begin(); p!=parts.end(); ++p)
                        p->exportToFile((outputFileName % (*thr++)).str( ));

                }

                //Destroy the inner data
                track.close();

                //display radius no more in use
                //track.displayRadius = displayRadius;
            }
        }
        Tracker::saveWisdom(config_path+"/wisdom.fftw");
    }
    catch(exception& e)
    {
        cout << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
