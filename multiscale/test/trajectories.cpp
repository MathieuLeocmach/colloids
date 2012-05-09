#define BOOST_TEST_DYN_LINK

#include "../src/traj.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/array.hpp>

using namespace Colloids;

BOOST_AUTO_TEST_SUITE( trajectories )
	BOOST_AUTO_TEST_CASE( traj )
	{
		Traj traj(-1, 5);
		BOOST_CHECK_EQUAL(traj.get_start(), -1);
		BOOST_CHECK_EQUAL(traj.get_finish(), 0);
		BOOST_CHECK_EQUAL(traj[-1], 5);
		traj.push_back(2);
		BOOST_CHECK_EQUAL(traj.get_start(), -1);
		BOOST_CHECK_EQUAL(traj.get_finish(), 1);
		BOOST_CHECK_EQUAL(traj[-1], 5);
		BOOST_CHECK_EQUAL(traj[0], 2);
	}
	BOOST_AUTO_TEST_CASE( trajindex )
	{
		TrajIndex ti(2);
		BOOST_REQUIRE_EQUAL(ti.nbFrames(), 1);
		BOOST_REQUIRE_EQUAL(ti.size(), 2);
		BOOST_REQUIRE_EQUAL(ti.getInverse(0).size(), 2);
		for(size_t p = 0; p<ti.getInverse(0).size(); ++p)
			BOOST_CHECK_EQUAL(ti.getInverse(0)[p], p);
		for(size_t tr=0; tr<ti.size(); ++tr)
			BOOST_CHECK_EQUAL(ti[tr][0], tr);
		std::vector<size_t> from(6), to(6);
		std::fill_n(from.begin(), 3, 0);
		std::fill_n(from.rbegin(), 3, 1);
		to[0] = 0;
		to[1] = 1;
		to[2] = 2;
		to[3] = 0;
		to[4] = 1;
		to[5] = 2;
		std::vector<double> distances(6);
		distances[0] = 1,1;
		distances[1] = 1.2;
		distances[2] = 0.9;
		distances[3] = 0.1;
		distances[4] = 1.1;
		distances[5] = 1.0;
		//case where we have more new positions than old ones
		ti.add_Frame(4, distances, from, to);
		BOOST_REQUIRE_EQUAL(ti.nbFrames(), 2);
		BOOST_REQUIRE_EQUAL(ti.getInverse(1).size(), 4);
		BOOST_REQUIRE_EQUAL(ti[1].size(), 2);
		BOOST_CHECK_EQUAL(ti[1][1], 0);
		BOOST_CHECK_EQUAL(ti.getTraj(1,0), 1);
		BOOST_REQUIRE_EQUAL(ti[0].size(), 2);
		BOOST_CHECK_EQUAL(ti[0][1], 2);
		BOOST_CHECK_EQUAL(ti.getTraj(1,2), 0);
		BOOST_REQUIRE(ti[2].exist(1));
		BOOST_CHECK_EQUAL(ti[2][1], 1);
		BOOST_CHECK_EQUAL(ti.getTraj(1,1), 2);
		BOOST_REQUIRE(ti[3].exist(1));
		BOOST_CHECK_EQUAL(ti[3][1], 3);
		BOOST_CHECK_EQUAL(ti.getTraj(1,3), 3);
		BOOST_REQUIRE_EQUAL(ti.size(), 4);

		//case where we have less positions than old ones
		ti.add_Frame(3, distances, from, to);
		BOOST_REQUIRE_EQUAL(ti.nbFrames(), 3);
		BOOST_REQUIRE_EQUAL(ti.size(), 5);
		BOOST_REQUIRE_EQUAL(ti.getInverse(2).size(), 3);
		BOOST_CHECK_EQUAL(ti[0].size(), 2);
		BOOST_REQUIRE_EQUAL(ti[1].size(), 3);
		BOOST_CHECK_EQUAL(ti[1][2], 2);
		BOOST_CHECK_EQUAL(ti.getTraj(2,2), 1);
		BOOST_CHECK(ti[2].exist(2));
		BOOST_CHECK_EQUAL(ti[2][2], 0);
		BOOST_CHECK_EQUAL(ti.getTraj(2,0), 2);
		BOOST_CHECK(!ti[3].exist(2));
		BOOST_CHECK(ti[4].exist(2));
		BOOST_CHECK_EQUAL(ti[4][2], 1);
		BOOST_CHECK_EQUAL(ti.getTraj(2,1), 4);

		//load a smaller frame
		boost::array<size_t, 2> fr0 = {{2, 1}};
		ti.add_Frame(fr0.begin(), fr0.end());
		BOOST_REQUIRE_EQUAL(ti.nbFrames(), 4);
		BOOST_REQUIRE_EQUAL(ti.size(), 5);
		BOOST_REQUIRE_EQUAL(ti.getInverse(3).size(), 2);
		BOOST_CHECK_EQUAL(ti[0].size(), 2);
		BOOST_REQUIRE_EQUAL(ti[1].size(), 4);
		BOOST_CHECK_EQUAL(ti[1][3], 1);
		BOOST_CHECK_EQUAL(ti.getTraj(3,1), 1);
		BOOST_CHECK(ti[2].exist(3));
		BOOST_CHECK_EQUAL(ti[2][3], 0);
		BOOST_CHECK_EQUAL(ti.getTraj(3,0), 2);
		BOOST_CHECK(!ti[3].exist(3));
		BOOST_CHECK(!ti[4].exist(3));

		//load a larger frame
		boost::array<size_t, 3> fr1 = {{2, 1, 5}};
		ti.add_Frame(fr1.begin(), fr1.end());
		BOOST_REQUIRE_EQUAL(ti.nbFrames(), 5);
		BOOST_REQUIRE_EQUAL(ti.size(), 6);
		BOOST_REQUIRE_EQUAL(ti.getInverse(4).size(), 3);
		BOOST_REQUIRE_EQUAL(ti[1].size(), 5);
		BOOST_CHECK_EQUAL(ti[1][4], 1);
		BOOST_CHECK_EQUAL(ti.getTraj(4,1), 1);
		BOOST_REQUIRE(ti[2].exist(4));
		BOOST_CHECK_EQUAL(ti[2][4], 0);
		BOOST_CHECK_EQUAL(ti.getTraj(4,0), 2);
		BOOST_REQUIRE(ti[5].exist(4));
		BOOST_CHECK_EQUAL(ti[5][4], 2);
		BOOST_CHECK_EQUAL(ti.getTraj(4,2), 5);

		//frame with no link
		ti.add_Frame(4, std::vector<double>(), std::vector<size_t>(), std::vector<size_t>());
		BOOST_REQUIRE_EQUAL(ti.nbFrames(), 6);
		BOOST_REQUIRE_EQUAL(ti.size(), 10);
		BOOST_REQUIRE_EQUAL(ti.getInverse(5).size(), 4);
		BOOST_REQUIRE_EQUAL(ti[1].size(), 5);
		BOOST_CHECK_EQUAL(ti[1][4], 1);
		BOOST_CHECK_EQUAL(ti.getTraj(4,1), 1);
		BOOST_REQUIRE(!ti[2].exist(5));
		BOOST_CHECK_EQUAL(ti[2][4], 0);
		BOOST_CHECK_EQUAL(ti.getTraj(4,0), 2);
		BOOST_REQUIRE(!ti[5].exist(5));
		BOOST_CHECK_EQUAL(ti[5][4], 2);
		BOOST_CHECK_EQUAL(ti.getTraj(4,2), 5);

	}
BOOST_AUTO_TEST_SUITE_END()
