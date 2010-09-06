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
#include "periodic.hpp"
//#include "pv.hpp"
#include "dynamicParticles.hpp"
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
		BondSet bonds;
		try
		{
            bonds = loadBonds(inputPath+".bonds");
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
			bonds = parts.getBonds();
			ofstream bondFile((inputPath+".bonds").c_str(), ios::out | ios::trunc);
			for(BondSet::const_iterator b=bonds.begin(); b!= bonds.end();++b)
				bondFile<<b->low()<<" "<<b->high()<<"\n";
			inside = parts.selectInside(1.3*radius);
			secondInside = parts.selectInside(2.0*1.3*radius);
		}
		cout<<bonds.size()<<" bonds"<<endl;
		const bool empty = inside.empty(), empty2 = secondInside.empty();
		if(empty)
			for(size_t p=0; p<parts.size(); ++p)
				inside.insert(inside.end(), p);
		if(empty2)
			for(size_t p=0; p<parts.size(); ++p)
				secondInside.insert(secondInside.end(), p);

		//calculate and export qlm
		vector<BooData> qlm, qlm_cg, qlm_first, qlm_second;
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
		{
			boost::progress_timer ti;
			cout<<"flip first shell "<<bonds.size()<<" bonds ";
			parts.getFlipBOOs(qlm, qlm_first, bonds);
		}
		{
			boost::progress_timer ti;
			BondSet bonds2nd = parts.getSecondShell();
			cout<<"flip second shell: "<<bonds2nd.size()<<" bonds in "<<ti.elapsed()<<"s"<<endl;
			parts.getFlipBOOs(qlm, qlm_second, bonds2nd);
		}

		/*{
			boost::progress_timer ti;
			set_union(
				bonds1551.begin(), bonds1551.end(),
				bonds2331.begin(), bonds2331.end(),
				inserter(bondsBoth, bondsBoth.end())
				);
			cout<<"both pairs: "<<bondsBoth.size()<<" in "<<ti.elapsed()<<"s"<<endl;
			parts.getFlipBOOs(qlm, qlm_flip, bondsBoth);
		}*/
		//BondSet bonds2331 = parts.get2331pairs();
		vector<size_t> nb_bonds6(parts.size(), 0), nb_bonds6_ref = nb_bonds6, nb_bonds6_rot = nb_bonds6;
		vector<double> sim(bonds.size(),-2.0), sim_ref=sim, sim_rot=sim;
		{
			boost::progress_timer ti;
			cout<<"boo product on first shell bonds ";
			int c=-1;
			for(BondSet::const_iterator b=bonds.begin(); b!=bonds.end(); ++b)
			{
				c++;
				if(qlm[b->low()][0]==0.0 || qlm[b->high()][0]==0.0
					|| !binary_search(inside.begin(), inside.end(), b->low())
					|| !binary_search(inside.begin(), inside.end(), b->high())
					)
					continue;
				BooData //rotated = qlm[b->high()].rotate_by_Pi(parts.getDiff(b->low(), b->high())),
					//reflected = qlm[b->high()].reflect(parts.getDiff(b->low(), b->high()), 6);
					reflected = qlm[b->high()].rotate_by_Pi(parts.getDiff(b->low(), b->high()));

				sim[c] = qlm[b->low()].normedProduct(qlm[b->high()], 6);
				//sim_rot[c] = qlm[b->low()].normedProduct(rotated, 6);
				sim_ref[c] = qlm[b->low()].normedProduct(reflected, 6);

				if(sim[c]>0.5)
				{
					nb_bonds6[b->low()]++;
					nb_bonds6[b->high()]++;
				}
				/*if(sim_rot[c]>0.5)
				{
					nb_bonds6_rot[b->low()]++;
					nb_bonds6_rot[b->high()]++;
				}*/
				if(sim_ref[c]>0.5)
				{
					nb_bonds6_ref[b->low()]++;
					nb_bonds6_ref[b->high()]++;
				}
				/*cout<<"\nleft";
				for(size_t m=0; m<=6;++m)
					cout<<"\t"<<qlm[b->low()](6,m);
				cout<<endl;
				cout<<"right";
				for(size_t m=0; m<=6;++m)
					cout<<"\t"<<qlm[b->high()](6,m);
				cout<<endl;
				cout<<"flipped";
				for(size_t m=0; m<=6;++m)
					cout<<"\t"<<reflected(6,m);
				cout<<endl;
				complex<double> sum(0.0,0.0);
				for(int m=-6; m<=6; ++m)
					sum += qlm[b->low()](6,m)*conj(reflected(6,m));
				cout<<"prod\t"<<sum<<endl;
				cout<<"norm left\t"<<qlm[b->low()].getSumNorm(6)<<endl;
				cout<<"norm right\t"<<reflected.getSumNorm(6)<<endl;
				cout<<"normed prod\t"<<qlm[b->low()].normedProduct(reflected, 6)<<endl;
				cout<<"self normed prod\t"<<qlm[b->low()].normedProduct(qlm[b->low()], 6)<<endl;*/
			}
		}

		vector<double> q6(parts.size(), 0.0), w6(parts.size(), 0.0);
		for(size_t p=0; p<parts.size(); ++p)
			qlm[p].getInvarients(6, q6[p], w6[p]);
		vector<double> Q6(parts.size(), 0.0), W6(parts.size(), 0.0);
		for(size_t p=0; p<parts.size(); ++p)
			qlm_cg[p].getInvarients(6, Q6[p], W6[p]);

		ofstream simvtkFile((head+"_sim"+neck+".vtk").c_str(), ios::out | ios::trunc);
		parts.toVTKstream(simvtkFile, "similarities");
		toVTKstream(simvtkFile, bonds);
		simvtkFile<<"POINT_DATA "<<parts.size()<<"\n";
		simvtkFile<<"SCALARS nb_bonds_raw double\n"
				"LOOKUP_TABLE default\n";
		copy(
			nb_bonds6.begin(), nb_bonds6.end(),
			ostream_iterator<double>(simvtkFile,"\n")
			);
		/*simvtkFile<<"SCALARS nb_bonds_rot double\n"
				"LOOKUP_TABLE default\n";
		copy(
			nb_bonds6_rot.begin(), nb_bonds6_rot.end(),
			ostream_iterator<size_t>(simvtkFile, "\n")
			);*/
		simvtkFile<<"SCALARS nb_bonds_ref double\n"
				"LOOKUP_TABLE default\n";
		copy(
			nb_bonds6_ref.begin(), nb_bonds6_ref.end(),
			ostream_iterator<size_t>(simvtkFile, "\n")
			);
		simvtkFile<<"SCALARS q6 double\n"
				"LOOKUP_TABLE default\n";
		copy(
			q6.begin(), q6.end(),
			ostream_iterator<double>(simvtkFile, "\n")
			);
		simvtkFile<<"SCALARS w6 double\n"
				"LOOKUP_TABLE default\n";
		copy(
			w6.begin(), w6.end(),
			ostream_iterator<double>(simvtkFile, "\n")
			);
		simvtkFile<<"SCALARS Q6 double\n"
				"LOOKUP_TABLE default\n";
		copy(
			Q6.begin(), Q6.end(),
			ostream_iterator<double>(simvtkFile, "\n")
			);
		simvtkFile<<"SCALARS W6 double\n"
				"LOOKUP_TABLE default\n";
		copy(
			W6.begin(), W6.end(),
			ostream_iterator<double>(simvtkFile, "\n")
			);

		simvtkFile<<"CELL_DATA "<<bonds.size()<<"\n";
		simvtkFile<<"SCALARS sim6_raw double\n"
				"LOOKUP_TABLE default\n";
		copy(
			sim.begin(), sim.end(),
			ostream_iterator<double>(simvtkFile, "\n")
			);
		/*simvtkFile<<"SCALARS sim6_rot double\n"
				"LOOKUP_TABLE default\n";
		copy(
			sim_rot.begin(), sim_rot.end(),
			ostream_iterator<double>(simvtkFile, "\n")
			);*/
		simvtkFile<<"SCALARS sim6_ref double\n"
				"LOOKUP_TABLE default\n";
		copy(
			sim_ref.begin(), sim_ref.end(),
			ostream_iterator<double>(simvtkFile, "\n")
			);
		simvtkFile.close();
		/*{
			boost::progress_timer ti;
			cout<<"boo product on 1551 bonds ";
			BondSet bonds1551 = parts.get1551pairs();
			ofstream simFile((head+"_1551"+neck+".sim").c_str(), ios::out | ios::trunc);
			ofstream simRotFile((head+"_rot_1551"+neck+".sim").c_str(), ios::out | ios::trunc);
			ofstream simRefFile((head+"_ref_1551"+neck+".sim").c_str(), ios::out | ios::trunc);
			for(BondSet::const_iterator b=bonds1551.begin(); b!=bonds1551.end(); ++b)
			{
				if(qlm[b->low()][0]==0.0 || qlm[b->high()][0]==0.0
					|| !binary_search(inside.begin(), inside.end(), b->low())
					|| !binary_search(inside.begin(), inside.end(), b->high())
					)
					continue;
				BooData rotated = qlm[b->high()].rotate_by_Pi(parts.getDiff(b->low(), b->high())),
					reflected = qlm[b->high()].reflect(parts.getDiff(b->low(), b->high()));
				for(int l=0; l<=10; l+=2)
				{
					simFile<<qlm[b->low()].normedProduct(qlm[b->high()], l)<<"\t";
					simRotFile<<qlm[b->low()].normedProduct(rotated, l)<<"\t";
					simRefFile<<qlm[b->low()].normedProduct(reflected, l)<<"\t";
				}
				simFile<<"\n";
				simRotFile<<"\n";
				simRefFile<<"\n";
			}
			simFile.close();
			simRotFile.close();
			simRefFile.close();
		}
		{
			boost::progress_timer ti;
			cout<<"boo product on 2331 bonds ";
			BondSet bonds2331 = parts.get2331pairs();
			ofstream simFile((head+"_2331"+neck+".sim").c_str(), ios::out | ios::trunc);
			ofstream simRotFile((head+"_rot_2331"+neck+".sim").c_str(), ios::out | ios::trunc);
			ofstream simRefFile((head+"_ref_2331"+neck+".sim").c_str(), ios::out | ios::trunc);
			for(BondSet::const_iterator b=bonds2331.begin(); b!=bonds2331.end(); ++b)
			{
				if(qlm[b->low()][0]==0.0 || qlm[b->high()][0]==0.0
					|| !binary_search(inside.begin(), inside.end(), b->low())
					|| !binary_search(inside.begin(), inside.end(), b->high())
					)
					continue;
				BooData rotated = qlm[b->high()].rotate_by_Pi(parts.getDiff(b->low(), b->high())),
					reflected = qlm[b->high()].reflect(parts.getDiff(b->low(), b->high()));
				for(int l=0; l<=10; l+=2)
				{
					simFile<<qlm[b->low()].normedProduct(qlm[b->high()], l)<<"\t";
					simRotFile<<qlm[b->low()].normedProduct(rotated, l)<<"\t";
					simRefFile<<qlm[b->low()].normedProduct(reflected, l)<<"\t";
				}
				simFile<<"\n";
				simRotFile<<"\n";
				simRefFile<<"\n";
			}
			simFile.close();
			simRotFile.close();
			simRefFile.close();
		}*/

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

		ofstream cloud_1File((head+"_flip"+neck+".cloud").c_str(), ios::out | ios::trunc);
		cloud_1File<<"#Q4\tQ6\tW4\tW6"<<endl;
		transform(
			qlm_first.begin(), qlm_first.end(),
			ostream_iterator<string>(cloud_1File,"\n"),
			cloud_exporter()
			);
		cloud_1File.close();

		ofstream cloud_2File((head+"_flip2"+neck+".cloud").c_str(), ios::out | ios::trunc);
		cloud_2File<<"#Q4\tQ6\tW4\tW6"<<endl;
		transform(
			qlm_second.begin(), qlm_second.end(),
			ostream_iterator<string>(cloud_2File,"\n"),
			cloud_exporter()
			);
		cloud_2File.close();

		/*ofstream cloud_flipFile((head+"_flip"+neck+".cloud").c_str(), ios::out | ios::trunc);
		cloud_flipFile<<"#Q4\tQ6\tW4\tW6"<<endl;
		transform(
			qlm_flip.begin(), qlm_flip.end(),
			ostream_iterator<string>(cloud_flipFile,"\n"),
			cloud_exporter()
			);
		cloud_flipFile.close();*/
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

