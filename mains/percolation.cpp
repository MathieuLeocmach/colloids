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
**/

#include "../dynamicClusters.hpp"

using namespace std;
using namespace Colloids;

int main(int argc, char ** argv)
{
    if(argc<3)
    {
        cout << "Syntax : percolation [path]filename.dat radius range" << endl;
        cout << "\tOR"<<endl;
        cout << "\tpercolation [path]filename.traj range" << endl;
        cout << " range is in diameter unit" << endl;
        return EXIT_FAILURE;
    }

    cout << "Cluster percolation" << endl;
    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.rfind("."));
    const string ext = filename.substr(inputPath.size()+1);
    const string path = filename.substr(0, filename.find_last_of("/\\")+1);

    double radius,range;

	try
	{
		DynamicParticles parts(filename);
		//cout << parts.trajectories.size() << " particles"<<endl;
		radius = parts.radius;
		range = atof(argv[2]);

		//fetch the file name pattern directly in the file.traj
		string pattern, token;
		size_t offset, size;
		{
			ifstream trajfile(filename.c_str(), ios::in);
			getline(trajfile, pattern);
			getline(trajfile, pattern); //pattern is on the 2nd line
			getline(trajfile, token); //token is on the 3rd line
			trajfile >> offset >> size;
			trajfile.close();
		}
		FileSerie datSerie(path+pattern, token, size, offset),
			bondSerie = datSerie.changeExt(".bonds"),
			clustSerie = datSerie.changeExt(".cluster");


		//create the population of trajectories to split into clusters : all trajectories
		set<size_t> alltraj;
		for(size_t tr=0;tr<parts.trajectories.size();++tr)
			alltraj.insert(alltraj.end(),tr);

		//neighbour list at each time
		try
		{
			for(size_t t=0; t<parts.getNbTimeSteps();++t)
			{
				BondSet bonds = loadBonds(bondSerie%t);
				parts.positions[t].makeNgbList(bonds);
			}
		}
		catch(invalid_argument& e)
		{
			for(size_t t=0; t<parts.getNbTimeSteps();++t)
			{
				parts.positions[t].makeRTreeIndex();
				parts.positions[t].makeNgbList(1.3);
				BondSet bonds = parts.positions[t].getBonds();
				ofstream bondFile((inputPath+".bonds").c_str(), ios::out | ios::trunc);
				for(BondSet::const_iterator b=bonds.begin(); b!= bonds.end();++b)
					bondFile<<b->low()<<" "<<b->high()<<"\n";
			}
		}


		cout<<"finding the clusters ... ";
		DynamicClusters clusters(parts, alltraj);
		cout<<"done!"<<endl;

		clusters.save(clustSerie);

	}
	catch(const exception &e)
	{
		cerr<< e.what()<<endl;
		return EXIT_FAILURE;
	}
    return EXIT_SUCCESS;
}

