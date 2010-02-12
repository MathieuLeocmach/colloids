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

#include "../dynamicParticles.hpp"

using namespace std;
using namespace Colloids;

int errorMessage()
{
    cout << "traj2vtk [path]filename.traj (AveragingInterval)" << endl;
    return EXIT_FAILURE;
}

int main(int argc, char ** argv)
{
    if(argc<2) return errorMessage();

    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    const string ext = filename.substr(filename.find_last_of(".")+1);
    const string path = filename.substr(0, filename.find_last_of("/\\")+1);

    //const double shell = 1.3;

    try
    {
		DynamicParticles parts(filename);
		const size_t stop = parts.getNbTimeSteps()-1;
		cout << parts.trajectories.size() << " particles in "<<stop+1<<" time steps ... ";

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
		FileSerie datSerie(path+pattern, token, size, offset);

		cout<<"remove drift ... ";
		parts.removeDrift();

		size_t tau;
		if(argc<3)
		{
			cout<<"isf to find relaxation time ... ";
			//Dynamics calculation and export
			vector<double> msd;
			vector< vector<double> > isf;
			parts.makeDynamics(msd, isf);
			ofstream msd_f((inputPath + ".msd").c_str());
			msd_f << "#t\tMSD"<<endl;
			for(size_t t=0; t<msd.size(); ++t)
				msd_f << t*parts.dt<<"\t"<< msd[t]<<endl;

			ofstream isf_f((inputPath + ".isf").c_str());
			isf_f << "#t\tx\ty\tz\tav"<<endl;
			for(size_t t=0; t<msd.size(); ++t)
			{
				isf_f << t*parts.dt;
				for(size_t d=0; d<4; ++d)
					isf_f<<"\t"<< isf[t][d];
				isf_f<<endl;
			}

			//Find relaxation time
			if(isf.back().back() > exp(-1.0))
				tau = parts.getNbTimeSteps();
			else
				tau = parts.getNbTimeSteps() - (upper_bound(isf.back().rbegin(),isf.back().rend(),exp(-1.0))-isf.back().rbegin());
		}
		else
			sscanf(argv[2],"%u",&tau);

		cout<<"relaxation time is "<<tau<<" steps, ie "<<tau*parts.dt<<"s"<<endl;
		//const size_t halfInterval = tau/2;

		//name pattern of the .cloud files
		cout<<"Load and average bond orientational order ... ";
		FileSerie booSerie = datSerie.changeExt(".cloud"),
				 SbooSerie = booSerie.addPostfix("_space");

		boost::multi_array<double, 2> qw;
		vector<ScalarDynamicField> scalars(9, ScalarDynamicField(parts.trajectories, tau, 0));
		scalars.reserve(9);
		for(size_t i=0;i<8;++i)
			scalars[i].name = string(i/4?"cg":"")+string((i/2)%2?"W":"Q")+string(i%2?"6":"4");
		for(size_t t=0; t<parts.getNbTimeSteps(); ++t)
		{
			//read raw boo from file
			parts.positions[t].loadBoo(booSerie%t, qw);
			//bin into the average
			cout<<"t="<<t<<" bin raw ... ";
			for(size_t i=0; i<qw.begin()->size(); ++i)
				scalars[i].push_back(ScalarField(qw.begin(), qw.end(), "", i));
			//same for coarse-grained boo
			cout<<" bin cg ... ";
			parts.positions[t].loadBoo(SbooSerie%t, qw);
			//bin into the average
			for(size_t i=0; i<qw.begin()->size(); ++i)
				scalars[i+4].push_back(ScalarField(qw.begin(), qw.end(), "", i));
		}

		cout<<"Neighbour lists ... ";
		FileSerie bondSerie = datSerie.changeExt(".bonds");
		for(size_t t=0; t<parts.getNbTimeSteps(); ++t)
			parts.positions[t].makeNgbList(loadBonds(bondSerie%t));

		cout<<"Lost neighbours ... ";
		boost::ptr_vector< vector<double> > lngb(parts.getNbTimeSteps());
		for(size_t t=0; t<parts.getNbTimeSteps(); ++t)
		{
			lngb.push_back(new vector<double>(parts.trajectories.size()));
			lngb.back() = parts.getNbLostNgbs(t, tau/2);
		}
		scalars.back().assign(lngb);
		scalars.back().name= "lostNgb";


		cout<<"Velocities ... ";
		boost::ptr_vector< vector<Coord> > vel(parts.getNbTimeSteps());
		for(size_t t=0; t<parts.getNbTimeSteps(); ++t)
		{
			vel.push_back(new vector<Coord>(parts.trajectories.size(), Coord(0.0,3)));
			vel.back() = parts.velocities(t, tau/2);
		}
		vector<VectorDynamicField> vectors(1, VectorDynamicField(parts.trajectories, vel, "V"));

		cout<<"export VTK ... ";
		FileSerie vtkSerie = datSerie.addPostfix("_dynamic", ".vtk");
		parts.exportToVTK(vtkSerie, scalars, vectors);
		cout<<endl;
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

