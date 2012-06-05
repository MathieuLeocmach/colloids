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
		//tests both even and odd sizes
		boost::array<int,2> sizes ={{11,10}};
		boost::array<float,6> spectrum;
		for (boost::array<int,2>::const_iterator s=sizes.begin(); s!=sizes.end(); ++s)
		{
			Convolver co(*s);
			BOOST_REQUIRE_EQUAL(co.size(), *s);
			BOOST_REQUIRE_EQUAL(co.fourier_size(), 6);
			std::vector<float> input(*s);
			//spectrum from empty input
			std::fill(input.begin(), input.end(), 0.0f);
			co.spectrum(&input[0], 1, &spectrum[0]);
			BOOST_CHECK_EQUAL(*std::max_element(spectrum.begin(), spectrum.end()), 0);
			BOOST_CHECK_EQUAL(*std::min_element(spectrum.begin(), spectrum.end()), 0);
			//step input
			input[0]=1;
			input[1]=1;
			co.spectrum(&input[0], 1, &spectrum[0]);
			BOOST_CHECK_EQUAL(spectrum[0], 4.0);
			//striding input
			std::vector<float> input2(2*(*s));
			std::fill(input2.begin(), input2.end(), 0.0f);
			input2[0]=1;
			input2[1]=1;
			input2[2]=1;
			co.spectrum(&input2[0], 2, &spectrum[0]);
			BOOST_CHECK_EQUAL(spectrum[0], 4.0);
			//simplest convolution: remove of the DC coefficient = subtract the average
			boost::array<float,6> kernel ={{0,1,1,1,1,1}};
			co(&input[0], 1, &kernel[0]);
			BOOST_CHECK_CLOSE(input[0], 1-2.0 / *s, 1e-5);
			BOOST_CHECK_CLOSE(input[2], -2.0 / *s, 1e-5);
		}
		//for a size of 10 the result is exactly null
		BOOST_CHECK_EQUAL(spectrum[5], 0.0);
		fftwf_cleanup();
	}
BOOST_AUTO_TEST_SUITE_END()

