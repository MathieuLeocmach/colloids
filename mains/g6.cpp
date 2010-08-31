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
#include "periodic.hpp"

using namespace std;
using namespace Colloids;

int main(int argc, char ** argv)
{
	invalid_argument er(
		"Bond angle correlation function\n"
#ifdef use_periodic
		"Syntax : periodic_g6 [path]filename radius NbOfBins range Nb dx dy dz [mode=0 [l=6]]\n"
#else
		"Syntax : g6 [path]filename radius NbOfBins range [mode=0 [l=6 [zmin zmax]]]\n"
#endif
		" range is in diameter unit\n"
		" mode\t0 raw boo\n"
		"     \t1 coarse grained boo\n"
		);

    try
    {

		if(argc<5) throw er;

		const string filename(argv[1]);
		const string inputPath = filename.substr(0,filename.find_last_of("."));
		const double radius = atof(argv[2]),
				nbDiameterCutOff = atof(argv[4]);
		const size_t Nbins = atoi(argv[3]);

		//construct the particle container out of the datafile
#ifdef use_periodic
		if(argc<9) throw er;
		const size_t Nb = atoi(argv[5]);
		BoundingBox b;
		for(size_t d=0;d<3;++d)
		{
			b.edges[d].first=0.0;
			b.edges[d].second = atof(argv[6+d]);
		}
		const bool mode = (argc<10)?0:atoi(argv[9]);
		const size_t l = (argc<11)?6:atoi(argv[10]);
		PeriodicParticles Centers(Nb,b,filename,radius);
		cout << "With periodic boundary conditions"<<endl;
#else
        const bool mode = (argc<6)?0:atoi(argv[5]);
        const size_t l = (argc<7)?0:atoi(argv[6]);
		Particles Centers(filename,radius);
#endif
		cout << Centers.size() << " particles ... ";
		Centers.makeRTreeIndex();

		//read q6m
		vector<BooData> rawBoo, allBoo;
		Centers.load_qlm(inputPath+".qlm", rawBoo);
		if(mode)
		{
		    Centers.makeNgbList(loadBonds(inputPath+".bonds"));
		    Centers.getCgBOOs(Centers.selectInside(4.0*radius*1.3), rawBoo, allBoo);
		}
		else
            allBoo=rawBoo;

		//create binner
		Particles::GlBinner binner(Centers, Nbins, nbDiameterCutOff, allBoo, l);

		//select particles and bin them
		vector<size_t> inside = Centers.selectInside(2.0*radius*(nbDiameterCutOff+(1+mode)*1.3/2));
		//Case where we ask to discard all points out of the interval [zmin, zmax]
		if(argc>8)
		{
		    const double zmin = atof(argv[6]),
                        zmax = atof(argv[7]);
            vector<size_t> in;
            in.reserve(inside.size());
            for(vector<size_t>::const_iterator p=inside.begin(); p!=inside.end(); ++p)
                if(Centers[*p][2]>=zmin && Centers[*p][2]<=zmax)
                    in.push_back(*p);
            in.swap(inside);
		}
		binner << inside;

		binner.normalize(inside.size());
		cout << " done !" << endl;
		ostringstream rdffilename;
		rdffilename << inputPath << (mode?".cg":".g") << l;
		ofstream rdfFile(rdffilename.str().c_str(), ios::out | ios::trunc);
		rdfFile << "#r\tg"<<l<<"(r)\tg(r)"<<endl;
		for(size_t r=0; r<Nbins; ++r)
			rdfFile<< r/(double)Nbins*nbDiameterCutOff <<"\t"<< binner.gl[r] <<"\t"<< binner.g[r] <<"\n";

    }
    catch(const std::exception &e)
    {
        cerr<<e.what() << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

