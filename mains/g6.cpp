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
		"Syntax : periodic_g6 [path]filename NbOfBins Nb dx dy dz [mode=0 [l=6]]\n"
#else
		"Syntax : g6 [path]filename radius NbOfBins range [mode=0 [l=6 [zmin zmax]]]\n"
		" range is in diameter unit\n"
#endif
		" mode\t0 raw boo\n"
		"     \t1 coarse grained boo\n"
		);

    try
    {

		if(argc<3) throw er;

		const string filename(argv[1]);
		const string inputPath = filename.substr(0,filename.find_last_of("."));

		vector<size_t> inside;

		//construct the particle container out of the datafile
#ifdef use_periodic
		if(argc<7) throw er;
		const size_t Nbins = atoi(argv[2]);
		const size_t Nb = atoi(argv[3]);
		double nbDiameterCutOff =atof(argv[4]);
		BoundingBox b;
		for(size_t d=0;d<3;++d)
		{
			b.edges[d].first=0.0;
			b.edges[d].second = atof(argv[4+d]);
			if(nbDiameterCutOff>b.edges[d].second)
                nbDiameterCutOff=b.edges[d].second;
		}
		nbDiameterCutOff /=4.0;
		const bool mode = (argc<8)?0:atoi(argv[7]);
		const size_t l = (argc<9)?6:atoi(argv[8]);
		PeriodicParticles Centers(Nb,b,filename,1.0);
		cout << "With periodic boundary conditions"<<endl;
#else
        const double radius = atof(argv[2]),
				nbDiameterCutOff = atof(argv[4]);
        const size_t Nbins = atoi(argv[3]);
        const bool mode = (argc<6)?0:atoi(argv[5]);
        const size_t l = (argc<7)?0:atoi(argv[6]);
		Particles Centers(filename,radius);
#endif
		cout << Centers.size() << " particles ... ";

		//read q6m
		vector<BooData> rawBoo, allBoo;
		Centers.load_qlm(inputPath+".qlm", rawBoo);
		if(mode)
		{
		    Centers.makeNgbList(loadBonds(inputPath+".bonds"));
		    //select all particles who's all neighbours' qlms are not null
		    vector<size_t> second_inside;
		    second_inside.reserve(Centers.size());
		    for(size_t p=0; p<Centers.getNgbList().size(); ++p)
		    {
		        bool hasNullNgb = false;
		        for(size_t n=0; n<Centers.getNgbList()[p].size(); ++n)
                    hasNullNgb = (hasNullNgb || rawBoo[Centers.getNgbList()[p][n]].isnull());
                if(!hasNullNgb)
                    second_inside.push_back(p);
		    }
		    cout<< second_inside.size()<< " have coarse grained qlms ... ";
		    Centers.getCgBOOs(second_inside, rawBoo, allBoo);
		}
		else
            allBoo=rawBoo;

#ifndef use_periodic
        vector<size_t> slab;
        slab.reserve(Centers.size());
        //Case where we ask to discard all points out of the interval [zmin, zmax]
        if(argc>8)
        {
            const double zmin = atof(argv[7]),
                        zmax = atof(argv[8]);
            //look for the particles that are in the slab and have non null qlm
            for(size_t p=0; p<Centers.size(); ++p)
                if(Centers[p][2]>zmin && Centers[p][2]<zmax && !allBoo[p].isnull())
                    slab.push_back(p);
        }
        else
            for(size_t p=0; p<Centers.size(); ++p)
                if(!allBoo[p].isnull())
                    slab.push_back(p);

        //reduce the centers list
        Particles cen;
        cen.radius = Centers.radius;
        cen.reserve(slab.size());
        for(size_t p=0; p<slab.size(); ++p)
            cen.push_back(Centers[slab[p]]);
        Centers.swap(cen);
        Centers.delNgbList();

        //reduce the qlm
        vector<BooData> boo(slab.size());
        for(size_t p=0; p<slab.size(); ++p)
            boo[p] = allBoo[slab[p]];
        allBoo.swap(boo);
#endif
        vector<size_t> nb(Nbins, 0);
        vector<double> gl(Nbins, 0.0);
        const double maxdistsq = pow(2.0*nbDiameterCutOff, 2);
        for(size_t p=0; p<Centers.size(); ++p)
            for(size_t q=p+1; q<Centers.size(); ++q)
            {
                Coord diff = Centers.getDiff(p, q);
                const double distsq = (diff*diff).sum();
                if(distsq<maxdistsq)
                {
                    const size_t r = sqrt(distsq) / (2.0*nbDiameterCutOff) * Nbins;
                    nb[r]++;
                    gl[r] += allBoo[p].innerProduct(allBoo[q], l);
                }
            }

        cout << " done !" << endl;
		ostringstream rdffilename;
		rdffilename << inputPath << (mode?".cg":".g") << l;
		ofstream rdfFile(rdffilename.str().c_str(), ios::out | ios::trunc);
		rdfFile << "#r\tg"<<l<<"(r)\tg(r)"<<endl;
		for(size_t r=0; r<Nbins; ++r)
			rdfFile<< r/(double)Nbins*nbDiameterCutOff <<"\t"<< gl[r] <<"\t"<< nb[r] <<"\n";

        //all particles are "inside"
        /*inside.reserve(Centers.size());
        for(size_t p=0; p<Centers.size(); ++p)
            inside[p] = p;
        cout << Centers.size() << " inside ... ";

        Centers.makeRTreeIndex();*/
		//create binner
		//Particles::GlBinner binner(Centers, Nbins, nbDiameterCutOff, allBoo, l);

		//select particles and bin them
		//vector<size_t> inside = Centers.selectInside(2.0*radius*(nbDiameterCutOff+(1+mode)*1.3/2));
		//Case where we ask to discard all points out of the interval [zmin, zmax]
		/*if(argc>8)
		{
		    const double zmin = atof(argv[7]),
                        zmax = atof(argv[8]);
            //cout<<"\tl="<<l<<"\tzmin="<<zmin<<"\tzmax="<<zmax;
            vector<size_t> in;
            in.reserve(inside.size());
            for(vector<size_t>::const_iterator p=inside.begin(); p!=inside.end(); ++p)
                if(Centers[*p][2]>=zmin && Centers[*p][2]<=zmax)
                    in.push_back(*p);
            in.swap(inside);
		}*/

		/*binner << inside;

		//binner.normalize(inside.size());
		cout << " done !" << endl;
		ostringstream rdffilename;
		rdffilename << inputPath << (mode?".cg":".g") << l;
		ofstream rdfFile(rdffilename.str().c_str(), ios::out | ios::trunc);
		rdfFile << "#r\tg"<<l<<"(r)\tg(r)"<<endl;
		for(size_t r=0; r<Nbins; ++r)
			rdfFile<< r/(double)Nbins*nbDiameterCutOff <<"\t"<< binner.gl[r] <<"\t"<< binner.g[r] <<"\n";*/

    }
    catch(const std::exception &e)
    {
        cerr<<e.what() << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

