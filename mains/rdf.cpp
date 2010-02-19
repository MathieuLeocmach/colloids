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

//Define the preprocessor variable "use_periodic" if you want periodic boundary conditions
#include "../periodic.hpp"

using namespace std;
using namespace Colloids;

int main(int argc, char ** argv)
{
    try
    {

		if(argc<5)
		{
			cout << "Syntax : [periodic_]rdf [path]filename radius NbOfBins range" << endl;
			cout << " range is in diameter unit" << endl;
			return EXIT_FAILURE;
		}

		cout << "Radial Distribution function" << endl;
		const string filename(argv[1]);
		const string inputPath = filename.substr(0,filename.find_last_of("."));
		const double radius = atof(argv[2]),
				nbDiameterCutOff = atof(argv[4]);
		const size_t Nbins = atoi(argv[3]);

		//construct the particle container out of the datafile
	#ifdef use_periodic
		if(argc<9)
		{
			cout << "Syntax : periodic_rdf [path]filename radius NbOfBins range Nb dx dy dz" << endl;
			cout << " range is in diameter unit" << endl;
			return EXIT_FAILURE;
		}
		const size_t Nb = atoi(argv[5]);
		BoundingBox b;
		for(size_t d=0;d<3;++d)
		{
			b.edges[d].first=0.0;
			b.edges[d].second = atof(argv[6+d]);
		}
		PeriodicParticles Centers(Nb,b,filename,radius);
		cout << "With periodic boundary conditions"<<endl;
	#else
		Particles Centers(filename,radius);
	#endif
		cout << Centers.size() << " particles ... spatial indexing ... ";
		Centers.makeRTreeIndex();
		cout << "g(r) ... ";

		//get g(r)
		vector<double> g = Centers.getRdf(Nbins,nbDiameterCutOff);
		cout << " done !" << endl;

		ofstream output((inputPath + ".rdf").c_str(), ios::out | ios::trunc);
		output<<"#r\tg"<<endl;
		const double scale = Nbins/nbDiameterCutOff;
		for(size_t r=0;r<g.size();++r)
            output<< r/scale <<"\t"<< g[r] << "\n";
    }
    catch(const std::exception &e)
    {
        cerr<<e.what() << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
