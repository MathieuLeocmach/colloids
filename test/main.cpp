#include <iostream> 
#define BOOST_TEST_MODULE tracker test
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "../graphic/tracker.hpp"

using namespace Colloids;

/*void images_are_close(const cv::Mat &a, const cv::Mat &b, float precision)
{
	OctaveFinder::Image  M = cv::abs(a)+cv::abs(b);
	double peak = *std::max_element(M.begin(), M.end());
	BOOST_REQUIRE_CLOSE(cv::sum(a)[0], cv::sum(b)[0], precision/peak);
	cv::Mat_<float>  diff = cv::abs(a-b) / peak;
	cv::MatConstIterator_<float> u = diff.begin();
	for(int i=0; i<a.rows; i++)
		for(int j=0; j<a.cols; j++)
		{
			BOOST_REQUIRE_MESSAGE(*u<precision, "at x=" << i <<" y=" << j <<"\t"<<diff(i,j)*peak<<" > "<<precision);
			BOOST_CHECK_SMALL(*u, precision);
			u++;
		}
}*/

