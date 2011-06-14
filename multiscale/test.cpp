#define BOOST_TEST_MODULE multiscale test
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "octavefinder.hpp"

using namespace Colloids;

BOOST_AUTO_TEST_CASE( octave_constructors_test )
{
    {
        OctaveFinder finder;
        BOOST_CHECK_EQUAL(finder.width, 256);
        BOOST_CHECK_EQUAL(finder.height, 256);
        BOOST_CHECK_EQUAL(finder.n_layers, 3);
    }
    {
        OctaveFinder finder(128, 512, 1);
        BOOST_CHECK_EQUAL(finder.width, 128);
        BOOST_CHECK_EQUAL(finder.height, 512);
        BOOST_CHECK_EQUAL(finder.n_layers, 1);
    }
}
