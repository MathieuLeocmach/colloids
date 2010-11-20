/**
    Copyright 2008,2009,2010 Mathieu Leocmach

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

#include "dynamicParticles.hpp"

#include <boost/progress.hpp>

using namespace std;
using namespace Colloids;

int main(int argc, char ** argv)
{
    if(argc<3)
    {
    	cout << "Average the positions over time"<<endl;
    	cout << "avpos [path]filename.traj AveragingInterval" << endl;
		return EXIT_FAILURE;
    }

    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    const string ext = filename.substr(filename.find_last_of(".")+1);
    const string path = filename.substr(0, filename.find_last_of("/\\")+1);
    const size_t tau = atoi(argv[2]);

    try
    {
    	ostringstream os;
		os<<"av"<<tau<<"/";
		system(("mkdir "+path+os.str()).c_str());

    	double radius, dt;
		string pattern, token;
		size_t offset, size;
		{
			ifstream trajfile(filename.c_str(), ios::in);
			if(!trajfile.good())
				throw invalid_argument((filename+" doesn't exist").c_str() );
			trajfile >> radius >> dt;
			trajfile.ignore(1); //escape the endl
			getline(trajfile, pattern); //pattern is on the 2nd line
			getline(trajfile, token); //token is on the 3rd line
			trajfile >> offset >> size;
			trajfile.close();
		}
		//File series
		FileSerie datSerie(path+os.str()+pattern, token, size, offset);
		cout<<"load"<<endl;
    	DynamicParticles parts(filename);
    	cout<<"remove drift"<<endl;
    	parts.removeDrift();
		cout << parts.trajectories.size() << " particles in "<<size<<" time steps"<<endl;

    	//average positions
    	cout<<"create VectorDynamicField"<<endl;
    	VectorDynamicField pos(parts.trajectories, tau, "positions");
    	for(size_t t=0; t<size; ++t)
			pos.push_back(VectorField(parts.positions[t].begin(), parts.positions[t].end()));
		//export
		for(size_t t=0; t<size; ++t)
		{
			VectorField f = pos[t];
			copy(
				f.values.begin(), f.values.end(),
				parts.positions[t].begin()
				);
			parts.positions[t].exportToFile(datSerie%t);
		}
	}
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


