#define BOOST_TEST_DYN_LINK

#include "../graphic/serieTracker.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>
#include <boost/array.hpp>
#include <set>

using namespace Colloids;
using namespace boost::posix_time;

BOOST_AUTO_TEST_SUITE( serieTracker )

BOOST_AUTO_TEST_SUITE( serieTrackerConstructor )

    BOOST_AUTO_TEST_CASE( single_image )
    {
        boost::array<size_t, 4> xyzt = {{16, 16, 1, 1}};
        std::string inputFile = "test_input/single_image.tif";
        SerieTracker track(inputFile, xyzt, 1.0, 0, FFTW_ESTIMATE);
        BOOST_CHECK_EQUAL(track.getPattern(), inputFile);
        BOOST_CHECK(!track.has_depth());
        BOOST_CHECK(!track.has_time());
        xyzt[2] = 2;
        BOOST_CHECK_THROW(SerieTracker track2(inputFile, xyzt, 1.0, 0, FFTW_ESTIMATE), std::invalid_argument);
        xyzt[2] = 1;
        xyzt[3] = 2;
        BOOST_CHECK_THROW(SerieTracker track2(inputFile, xyzt, 1.0, 0, FFTW_ESTIMATE), std::invalid_argument);
    }
    
    BOOST_AUTO_TEST_CASE( zSeries )
    {
        boost::array<size_t, 4> xyzt = {{16, 16, 16, 1}};
        std::string inputFile = "test_input/ZSerie_z00.tif";
        SerieTracker track(inputFile, xyzt, 1.0, 0, FFTW_ESTIMATE);
        BOOST_CHECK(track.has_depth());
        BOOST_CHECK(!track.has_time());
        BOOST_CHECK_EQUAL(track.getPattern(), "test_input/ZSerie_z%|02d|.tif");
    }
    BOOST_AUTO_TEST_CASE( tSeries )
    {
        boost::array<size_t, 4> xyzt = {{16, 16, 1, 3}};
        std::string inputFile = "test_input/TSerie_t000.tif";
        SerieTracker track(inputFile, xyzt, 1.0, 0, FFTW_ESTIMATE);
        BOOST_CHECK(!track.has_depth());
        BOOST_CHECK(track.has_time());
        BOOST_CHECK_EQUAL(track.getPattern(), "test_input/TSerie_t%|03d|.tif");
    }
    BOOST_AUTO_TEST_CASE( ztSeries )
    {
        boost::array<size_t, 4> xyzt = {{16, 16, 1, 3}};
        std::string inputFile = "test_input/ZTSerie_z00_t000.tif";
        SerieTracker track(inputFile, xyzt, 1.0, 0, FFTW_ESTIMATE);
        BOOST_CHECK(track.has_depth());
        BOOST_CHECK(track.has_time());
        BOOST_CHECK_EQUAL(track.getPattern(), "test_input/ZTSerie_z%|02d|_t%|03d|.tif");
    }

BOOST_AUTO_TEST_SUITE_END() //serieTrackerConstructor

BOOST_AUTO_TEST_SUITE_END() //serieTracker
