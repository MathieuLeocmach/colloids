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
#include <boost/progress.hpp>

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

		const double radius = atof(argv[2]);
#ifdef use_periodic
		const size_t Nb = atoi(argv[3]);
		BoundingBox b;
		for(size_t d=0;d<3;++d)
		{
			b.edges[d].first=0.0;
			b.edges[d].second = atof(argv[4+d]);
		}
		PeriodicParticles parts(Nb,b,filename,radius);
#else
		Particles parts(filename, radius);
#endif
		//select interesting particles and load (or make) bonds
		vector<size_t> inside, secondInside;
		try
		{
            BondSet bonds = loadBonds(inputPath+".bonds");
            cout<<bonds.size()<<" bonds"<<endl;
            parts.makeNgbList(bonds);
			inside = parts.selectInside_noindex(1.3*radius);
			secondInside = parts.selectInside_noindex(2.0*1.3*radius);
		}
		catch(invalid_argument &e)
		{
		    cout<<"bond network ";
            boost::progress_timer ti;
			parts.makeRTreeIndex();
			parts.makeNgbList(1.3);
			BondSet bonds = parts.getBonds();
			cout<<bonds.size()<<" bonds"<<endl;
			ofstream bondFile((inputPath+".bonds").c_str(), ios::out | ios::trunc);
			for(BondSet::const_iterator b=bonds.begin(); b!= bonds.end();++b)
				bondFile<<b->low()<<" "<<b->high()<<"\n";
			inside = parts.selectInside(1.3*radius);
			secondInside = parts.selectInside(2.0*1.3*radius);
		}
		const bool empty = inside.empty(), empty2 = secondInside.empty();
		if(empty)
			for(size_t p=0; p<parts.size(); ++p)
				inside.insert(inside.end(), p);
		if(empty2)
			for(size_t p=0; p<parts.size(); ++p)
				secondInside.insert(secondInside.end(), p);

		//calculate and export qlm
		vector<BooData> qlm, qlm_cg, qlm_1551, qlm_2331, qlm_flip;
        {
            cout<<"boo without surface bonds for "<<inside.size()<<" particles ";
            boost::progress_timer ti;
            parts.getBOOs(qlm);
            if(!empty)
            	parts.removeOutside(inside, qlm);
        }
		{
            cout<<"coarse grained for "<<secondInside.size()<<" particles ";
            boost::progress_timer ti;
            parts.getCgBOOs(secondInside, qlm, qlm_cg);
		}
		BondSet bonds1551, bonds2331, bondsBoth;
		{
			boost::progress_timer ti;
			bonds1551 = parts.get1551pairs();
			cout<<"(1551) pairs: "<<bonds1551.size()<<" in "<<ti.elapsed()<<"s"<<endl;
			parts.getFlipBOOs(qlm, qlm_1551, bonds1551);
		}
		{
			boost::progress_timer ti;
			bonds2331 = parts.get2331pairs();
			cout<<"(2331) pairs: "<<bonds2331.size()<<" in "<<ti.elapsed()<<"s"<<endl;
			parts.getFlipBOOs(qlm, qlm_2331, bonds2331);
		}
		{
			boost::progress_timer ti;
			set_union(
				bonds1551.begin(), bonds1551.end(),
				bonds2331.begin(), bonds2331.end(),
				inserter(bondsBoth, bondsBoth.end())
				);
			cout<<"both pairs: "<<bondsBoth.size()<<" in "<<ti.elapsed()<<"s"<<endl;
			parts.getFlipBOOs(qlm, qlm_flip, bondsBoth);
		}

		//calculate and export invarients
		ofstream cloudFile((inputPath+".cloud").c_str(), ios::out | ios::trunc);
		cloudFile<<"#q4\tq6\tq8\tq10\tw4\tw6\tw8\tw10"<<endl;
		transform(
			qlm.begin(), qlm.end(),
			ostream_iterator<string>(cloudFile,"\n"),
			cloud_exporter()
			);
		cloudFile.close();

		ofstream cloud_cgFile((head+"_space"+neck+".cloud").c_str(), ios::out | ios::trunc);
		cloud_cgFile<<"#Q4\tQ6\tQ8\tQ10\tW4\tW6\tW8\tW10"<<endl;
		transform(
			qlm_cg.begin(), qlm_cg.end(),
			ostream_iterator<string>(cloud_cgFile,"\n"),
			cloud_exporter()
			);
		cloud_cgFile.close();

		ofstream cloud_1551File((head+"_1551"+neck+".cloud").c_str(), ios::out | ios::trunc);
		cloud_1551File<<"#Q4\tQ6\tW4\tW6"<<endl;
		transform(
			qlm_1551.begin(), qlm_1551.end(),
			ostream_iterator<string>(cloud_1551File,"\n"),
			cloud_exporter()
			);
		cloud_1551File.close();

		ofstream cloud_2331File((head+"_2331"+neck+".cloud").c_str(), ios::out | ios::trunc);
		cloud_2331File<<"#Q4\tQ6\tW4\tW6"<<endl;
		transform(
			qlm_2331.begin(), qlm_2331.end(),
			ostream_iterator<string>(cloud_2331File,"\n"),
			cloud_exporter()
			);
		cloud_2331File.close();

		ofstream cloud_flipFile((head+"_flip"+neck+".cloud").c_str(), ios::out | ios::trunc);
		cloud_flipFile<<"#Q4\tQ6\tW4\tW6"<<endl;
		transform(
			qlm_flip.begin(), qlm_flip.end(),
			ostream_iterator<string>(cloud_flipFile,"\n"),
			cloud_exporter()
			);
		cloud_flipFile.close();
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

