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

//Define the preprocessor variable "periodic" if you want periodic boundary conditions
#include "../periodic.hpp"
//#include "../pv.hpp"
#include "../dynamicParticles.hpp"

using namespace std;
using namespace Colloids;

int main(int argc, char ** argv)
{
#ifdef use_periodic
	if(argc<7)
	{
		cerr<<"Syntax : periodic_boo [path]filename.grv radius Nb Dx Dy Dz" << endl;
		return EXIT_FAILURE;
	}
#else
    if(argc<3)
	{
		cerr<<"Syntax : boo [path]filename.dat radius" << endl;
		return EXIT_FAILURE;
	}
#endif

	try
    {
		const string filename(argv[1]);
		const string inputPath = filename.substr(0,filename.find_last_of("."));
		const string ext = filename.substr(filename.find_last_of(".")+1);
		const string head = filename.substr(0,filename.rfind("_t"));
		const string neck = filename.substr(head.size(), inputPath.size()-head.size());

		double radius;
		sscanf(argv[2],"%lf",&radius);
#ifdef use_periodic
		size_t Nb;
		sscanf(argv[3],"%u",&Nb);
		BoundingBox b;
		for(size_t d=0;d<3;++d)
		{
			b.edges[d].first=0.0;
			sscanf(argv[4+d],"%lf",&b.edges[d].second);
		}
		PeriodicParticles parts(Nb,b,filename,radius);
#else
		Particles parts(filename, radius);
#endif
		//select interesting particles and load (or make) bonds
		set<size_t> inside, secondInside;
		BondList bonds = loadBonds(inputPath+".bonds");
		if(bonds.empty())
		{
			parts.makeRTreeIndex();
			parts.makeNgbList(1.3);
			bonds = parts.getBonds();
			ofstream bondFile((inputPath+".bonds").c_str(), ios::out | ios::trunc);
			for(deque<pair<size_t, size_t> >::const_iterator b=bonds.begin(); b!= bonds.end();++b)
				bondFile<<b->first<<" "<<b->second<<"\n";
			inside = parts.selectInside(1.3*radius);
			secondInside = parts.selectInside(2.0*1.3*radius);
		}
		else
		{
			parts.makeNgbList(bonds);
			inside = parts.selectInside_noindex(1.3*radius);
			secondInside = parts.selectInside_noindex(2.0*1.3*radius);
			cout<<"\tinside is "<<inside.size()<<" particles"<<endl;
			cout<<"second inside is "<<secondInside.size()<<" particles"<<endl;
		}

		//calculate and export qlm
		vector<BooData> qlm, qlm_cg;
		parts.getBOOs(inside, qlm);
		parts.getCgBOOs(secondInside, qlm, qlm_cg);
		ofstream qlmFile((inputPath+".qlm").c_str(), ios::out | ios::trunc);
		copy(
			qlm_cg.begin(), qlm_cg.end(),
			ostream_iterator<BooData>(qlmFile,"\n")
			);

		//calculate and export invarients
		ofstream cloudFile((inputPath+".cloud").c_str(), ios::out | ios::trunc);
		cloudFile<<"#Q4\tQ6\tW4\tW6"<<endl;
		transform(
			qlm.begin(), qlm.end(),
			ostream_iterator<string>(cloudFile,"\n"),
			cloud_exporter()
			);
		cloudFile.close();

		ofstream cloud_cgFile((head+"_space"+neck+".cloud").c_str(), ios::out | ios::trunc);
		cloud_cgFile<<"#Q4\tQ6\tW4\tW6"<<endl;
		transform(
			qlm_cg.begin(), qlm_cg.end(),
			ostream_iterator<string>(cloud_cgFile,"\n"),
			cloud_exporter()
			);
		cloud_cgFile.close();
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
