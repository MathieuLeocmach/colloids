/**
    Copyright 2010 Mathieu Leocmach

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
#include <boost/program_options.hpp>
#include "src/center.hpp"
#include <cmath>

namespace po = boost::program_options;
using namespace std;
using namespace Colloids;

int main(int ac, char ** av)
{
	try
    {
		double maxlength = 2.0;
		po::options_description cmdline_options("Generic options");
		cmdline_options.add_options()
			("length,l", po::value<double>(&maxlength)->default_value(2.0), "Maximum relative bond length")
			("number,n", po::value<size_t>(), "Number of particles (optional)")
			("verbose,v", "Output debugging information")
			("help", "Show this help and exit")
			;
		po::variables_map vm;
		po::store(po::parse_command_line(ac, av, cmdline_options), vm);
		po::notify(vm);

		if(!!vm.count("help"))
		{
			std::cout<< "Read x y z r on standard input and output bonds as i j distance"<<std::endl;
			std::cout << cmdline_options << std::endl;
			return EXIT_SUCCESS;
		}
		std::vector<Center3D> centers;
		if(!!vm.count("number"))
			centers.reserve(vm["number"].as<size_t>());
		//read stdin
		while(std::cin.good())
		{
			Center3D c;
			std::cin >> c[0] >> c[1] >> c[2] >>c.r;
			if(std::cin.good())
				centers.push_back(c);
		}
		//std::cout<<centers.size()<<" particles"<<std::endl;
		//maximum radius
		const double rmax = std::max_element(centers.begin(), centers.end(), compare_radii<3>())->r;
		const double mdsq = maxlength*maxlength;

		//Spatial index
		typedef RStarTree<size_t, 3, 4, 32, double> RTree;
		RTree tree;
		for(size_t p=0; p<centers.size();++p)
			tree.Insert(p, get_bb(centers[p]));

		//construct bonds
		for(size_t p=0; p<centers.size()-1;++p)
		{
			std::list<size_t> ngbs;
			tree.Query(
					typename RTree::AcceptOverlapping(get_bb_margin(
							centers[p],
							maxlength*(centers[p].r+rmax))
							),
					Gatherer<3>(ngbs));
			ngbs.sort();
			for(std::list<size_t>::const_iterator q=std::lower_bound(ngbs.begin(), ngbs.end(), p+1); q!=ngbs.end(); ++q)
			{
				//Distance
				const double distsq = centers[p] - centers[*q];
				if(distsq < mdsq * pow(centers[p].r + centers[*q].r, 2))
					std::cout << p <<"\t" << *q << "\t" << sqrt(distsq) << "\n";
			}
		}
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

