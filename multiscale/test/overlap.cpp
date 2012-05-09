#define BOOST_TEST_DYN_LINK

#include "../src/center.hpp"
#include <boost/test/unit_test.hpp>

using namespace Colloids;

BOOST_AUTO_TEST_SUITE( Overlap )
	BOOST_AUTO_TEST_CASE(Overlap_RTree)
	{
		typedef RStarTree<size_t, 2, 4, 32, double> RTree;
		std::vector<Center2D> centers;
		std::auto_ptr<RTree> tree = removeOverlapping(centers);
		BOOST_REQUIRE(centers.empty());
		Center2D c(0, 1, 0);
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 1);
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 1);
		c.intensity = -1.0;
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 1);
		BOOST_CHECK_CLOSE(centers[0].intensity, -1, 1e-9);
		c[0] = 3;
		c.intensity = -2.0;
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 2);
		BOOST_CHECK_CLOSE(centers[0].intensity, -2, 1e-9);
		BOOST_CHECK_CLOSE(centers[1].intensity, -1, 1e-9);
		c[0] = 2;
		c.intensity = 0;
		c.r = 1.1;
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 2);
		BOOST_CHECK_CLOSE(centers[0].intensity, -2, 1e-9);
		BOOST_CHECK_CLOSE(centers[1].intensity, -1, 1e-9);
		c.intensity = -1.5;
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 2);
		BOOST_CHECK_CLOSE(centers[0].intensity, -2, 1e-9);
		BOOST_CHECK_CLOSE(centers[1].intensity, -1, 1e-9);
		c.intensity = -3.0;
		centers.push_back(c);
		tree = removeOverlapping(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 1);
		BOOST_CHECK_CLOSE(centers[0].intensity, -3, 1e-9);
	}
	BOOST_AUTO_TEST_CASE(Overlap_Brute)
	{
		std::vector<Center2D> centers;
		removeOverlapping_brute_force(centers);
		BOOST_REQUIRE(centers.empty());
		Center2D c(0, 1, 0);
		centers.push_back(c);
		removeOverlapping_brute_force(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 1);
		centers.push_back(c);
		removeOverlapping_brute_force(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 1);
		c.intensity = -1.0;
		centers.push_back(c);
		removeOverlapping_brute_force(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 1);
		BOOST_CHECK_CLOSE(centers[0].intensity, -1, 1e-9);
		c[0] = 3;
		c.intensity = -2.0;
		centers.push_back(c);
		removeOverlapping_brute_force(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 2);
		BOOST_CHECK_CLOSE(centers[0].intensity, -2, 1e-9);
		BOOST_CHECK_CLOSE(centers[1].intensity, -1, 1e-9);
		c[0] = 2;
		c.intensity = 0;
		c.r = 1.1;
		centers.push_back(c);
		removeOverlapping_brute_force(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 2);
		BOOST_CHECK_CLOSE(centers[0].intensity, -2, 1e-9);
		BOOST_CHECK_CLOSE(centers[1].intensity, -1, 1e-9);
		c.intensity = -1.5;
		centers.push_back(c);
		removeOverlapping_brute_force(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 2);
		BOOST_CHECK_CLOSE(centers[0].intensity, -2, 1e-9);
		BOOST_CHECK_CLOSE(centers[1].intensity, -1, 1e-9);
		c.intensity = -3.0;
		centers.push_back(c);
		removeOverlapping_brute_force(centers);
		BOOST_REQUIRE_EQUAL(centers.size(), 1);
		BOOST_CHECK_CLOSE(centers[0].intensity, -3, 1e-9);
	}
BOOST_AUTO_TEST_SUITE_END()
