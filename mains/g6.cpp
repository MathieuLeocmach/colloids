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

int main(int argc, char ** argv)
{
	invalid_argument er(
		"Bond angle correlation function\n"
#ifdef use_periodic
		"Syntax : periodic_g6 [path]filename radius NbOfBins range Nb dx dy dz\n"
#else
		"Syntax : g6 [path]filename radius NbOfBins range\n"
#endif
		" range is in diameter unit");

    try
    {

		if(argc<5) throw er;

		const string filename(argv[1]);
		const string inputPath = filename.substr(0,filename.find_last_of("."));
		double radius,nbDiameterCutOff;
		sscanf(argv[2],"%lf",&radius);
		size_t Nbins;
		sscanf(argv[3],"%u",&Nbins);
		sscanf(argv[4],"%lf",&nbDiameterCutOff);

		//construct the particle container out of the datafile
#ifdef use_periodic
		if(argc<9) throw er;
		size_t Nb;
		sscanf(argv[5],"%u",&Nb);
		BoundingBox b;
		for(size_t d=0;d<3;++d)
		{
			b.edges[d].first=0.0;
			sscanf(argv[6+d],"%lf",&b.edges[d].second);
		}
		PeriodicParticles Centers(Nb,b,radius,filename);
		cout << "With periodic boundary conditions"<<endl;
#else
		IndexedParticles Centers(filename,radius);
#endif
		cout << Centers.size() << " particles ... ";

		//read q6m
		map<size_t,BooData> allBoo;
		Centers.load_q6m(inputPath+".q6m",allBoo);

		//create binner
		IndexedParticles::G6Binner binner(Centers,Nbins,nbDiameterCutOff,allBoo);

		//select particles and bin them
		set<size_t> inside = Centers.selectInside(2.0*radius*(nbDiameterCutOff+2));
		binner << inside;

		//get g6(r)
		deque< vector<double> > a;
		/*a.push_back(binner.g);
		a.push_back(binner.g6);*/


		binner.normalize(inside.size());
		cout << " done !" << endl;
		a.push_back(binner.g6);
		a.push_back(binner.g);
		saveTable(a.begin(),a.end(),inputPath + ".g6","r\tg6(r)\tg(r)",nbDiameterCutOff/Nbins);
		//saveRDF(binner.g6,inputPath + ".g6",((double)Nbins)/nbDiameterCutOff);
    }
    catch(const std::exception &e)
    {
        cerr<<e.what() << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

