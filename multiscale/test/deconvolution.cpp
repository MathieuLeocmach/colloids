/*
 * deconvolution.cpp
 *
 *  Created on: 4 juin 2012
 *      Author: mathieu
 */

#define BOOST_TEST_DYN_LINK

#include "../src/deconvolution.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>
#include <boost/array.hpp>


using namespace Colloids;
using namespace boost::posix_time;

BOOST_AUTO_TEST_SUITE( Deconvolution )
	BOOST_AUTO_TEST_CASE( convolver )
	{
		Convolver co(10, 1);
		BOOST_REQUIRE_EQUAL(co.size(), 10);
		BOOST_REQUIRE_EQUAL(co.fourier_size(), 6);
	}
BOOST_AUTO_TEST_SUITE_END()

