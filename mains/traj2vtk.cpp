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

#include <boost/progress.hpp>

using namespace std;
using namespace Colloids;

void export_lostNgb(const TrajIndex &trajectories, const size_t &tau, FileSerie &bondSerie, FileSerie &lngbSerie)
{
    const size_t size = trajectories.inverse.size();
    typedef vector<size_t>	Ngbs;
    typedef boost::ptr_vector<Ngbs> NgbFrame;
    typedef boost::ptr_vector<NgbFrame> DynNgbs;
    typedef list<size_t> ListNgb;
    typedef vector<ListNgb> ListNgbFrame;

    //what are the trajectories neighbouring the position p at time t
    DynNgbs dyn(size);
    for(size_t t=0; t<size;++t)
	{
		dyn.push_back(new NgbFrame(trajectories.inverse[t].size()));
		for(size_t p=0; p<trajectories.inverse[t].size();++p)
			dyn[t].push_back(new Ngbs());
	}


	for(size_t t=0; t<size;++t)
	{
		//Easily filled, disordered but easy to sort container. Bad memory perf. But it's just one frame.
		ListNgbFrame easy(dyn[t].size());
		//load bonds from file and bin their ends to neighbour list
		size_t a, b;
		ifstream f((bondSerie%t).c_str());
		while(f.good())
		{
			f>>a>>b;
			easy[a].push_back(trajectories.inverse[t][b]);
			easy[b].push_back(trajectories.inverse[t][a]);
		}
		f.close();

		#pragma omp parallel for shared(easy, dyn) schedule(dynamic)
		for(size_t p=0;p<easy.size();++p)
		{
			//sort each neighbour list
			easy[p].sort();
			//ensure uniqueness (could be asserted, but anyway...)
			easy[p].unique();
			//fill in the memory-efficient container
			dyn[t][p].assign(easy[p].begin(), easy[p].end());
		}
    }

    for(size_t t0=0; t0<size; ++t0)
	{
	    vector<size_t> lngb(trajectories.inverse[t0].size(),0);
	    #pragma omp parallel for shared(trajectories, lngb, t0) schedule(dynamic)
	    for(size_t p=0;p<trajectories.inverse[t0].size(); ++p)
	    {
	        const Traj &tr = trajectories[trajectories.inverse[t0][p]];
	        const size_t start = max(tr.start_time+tau/2, (size_t)t0)-tau/2,
					stop = min(tr.last_time(), (size_t)t0+tau/2);
            if(start!=stop)
            {
                ListNgb lost;
                set_difference(
                    dyn[stop][tr[stop]].begin(), dyn[stop][tr[stop]].end(),
                    dyn[start][tr[start]].begin(), dyn[start][tr[start]].end(),
                    back_inserter(lost)
                    );
                lngb[p] = lost.size();
            }
	    }
		ofstream f((lngbSerie%t0).c_str(), ios::out | ios::trunc);
		copy(
			lngb.begin(), lngb.end(),
			ostream_iterator<size_t>(f,"\n")
			);
		f.close();
	}
}

void export_timeBoo(const TrajIndex& trajectories, const size_t& tau, FileSerie &timeBooSerie, const string &prefix="")
{
	const size_t size = trajectories.inverse.size();
	FileSerie booSerie = timeBooSerie.changeExt(".cloud");
	boost::multi_array<double, 2> qw;
	vector<ScalarDynamicField> scalars(4, ScalarDynamicField(trajectories, tau, 0));
	for(size_t i=0;i<4;++i)
		scalars[i].name = prefix+string((i/2)%2?"W":"Q")+string(i%2?"6":"4");
	for(size_t t=0; t<size; ++t)
	{
		//read boo from file
		boost::array<size_t, 2> shape = {{trajectories.inverse[t].size(), 4}};
		qw.resize(shape);
		ifstream cloud((booSerie%t).c_str(), ios::in);
		string trash;
		//trashing the header
		getline(cloud, trash);
		copy(
			istream_iterator<double>(cloud), istream_iterator<double>(),
			qw.origin()
			);
		//bin into the average
		for(size_t i=0; i<qw.begin()->size(); ++i)
			scalars[i].push_back(ScalarField(qw.begin(), qw.end(), "", i));
	}
	//export
	for(size_t t=0; t<size; ++t)
	{
		vector<ScalarField> s(4, ScalarField("", trajectories.inverse[t].size()));
		for(int i=0; i<4; ++i)
			s[i] = scalars[i][t];
		ofstream f((timeBooSerie%t).c_str(), ios::out | ios::trunc);
		f<<"#Q4\tQ6\tW4\tW6"<<endl;
		for(size_t p=0; p<s.front().size(); ++p)
		{
			for(size_t i=0;i<4;++i)
				f << s[i][p]<<"\t";
			f<<"\n";
		}
		f.close();
	}
}

void export_volume(const TrajIndex& trajectories, const size_t& tau, FileSerie &volSerie, FileSerie &timeVolSerie)
{
	const size_t size = trajectories.inverse.size();
	ScalarDynamicField volume(trajectories, tau, 0, "Volume");
	for(size_t t=0; t<size; ++t)
	{
		vector<double> vol(trajectories.inverse[t].size());
		ifstream f((volSerie%t).c_str(), ios::in);
		copy(
			istream_iterator<double>(f), istream_iterator<double>(),
			vol.begin()
			);
		volume.push_back(ScalarField(vol.begin(), vol.end()));
	}
	//export
	for(size_t t=0; t<size; ++t)
	{
		ofstream f((timeVolSerie%t).c_str(), ios::out | ios::trunc);
		f<<volume[t];
		f.close();
	}
}


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
		cout << parts.trajectories.size() << " particles in "<<stop+1<<" time steps"<<endl;

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
			velSerie = datSerie.changeExt(".vel"),
			bondSerie = datSerie.changeExt(".bonds"),
			timeBooSerie = datSerie.changeExt(".boo"),
			timecgBooSerie = datSerie.addPostfix("_space",".boo"),
			lngbSerie = datSerie.changeExt(".lngb"),
			volSerie = datSerie.changeExt(".vol"),
			timeVolSerie = volSerie.addPostfix("_time"),
			vtkSerie = datSerie.addPostfix("_dynamic", ".vtk");

		cout<<"remove drift ... ";
		parts.removeDrift();

		size_t tau;
		if(argc<3)
		{
			cout<<"isf to find relaxation time ... ";
			vector< vector<double> > isf(4,vector<double>(size));
			ifstream in((inputPath + ".isf").c_str());
			if(in.good())
			{
				cout<<"read from file ... "<<endl;
				string s;
				getline(in, s);
				for(size_t t=0; t<isf.front().size(); ++t)
				{
					in >> s;
					for(size_t d=0; d<4; ++d)
						in >> isf[d][t];
				}
				in.close();
			}
			else
			{
				cout<<"calculate ... "<<endl;
				//Dynamics calculation and export
				vector<double> msd;

				parts.makeDynamics(msd, isf);
				ofstream msd_f((inputPath + ".msd").c_str());
				msd_f << "#t\tMSD"<<endl;
				for(size_t t=0; t<msd.size(); ++t)
					msd_f << t*parts.dt<<"\t"<< msd[t]<<"\n";

				msd_f.close();
				ofstream isf_f((inputPath + ".isf").c_str(), ios::out | ios::trunc);
				isf_f << "#t\tx\ty\tz\tav"<<endl;
				for(size_t t=0; t<isf.front().size(); ++t)
				{
					isf_f << t*parts.dt;
					for(size_t d=0; d<4; ++d)
						isf_f<<"\t"<< isf[d][t];
					isf_f<<"\n";
				}
				isf_f.close();
			}
			//Find relaxation time
			if(isf.back().back() > exp(-1.0))
				tau = parts.getNbTimeSteps()-1;
			else
				tau = parts.getNbTimeSteps()-1 - (upper_bound(isf.back().rbegin(),isf.back().rend(),exp(-1.0))-isf.back().rbegin());
		}
		else
			tau = atoi(argv[2]);

		cout<<"relaxation time is "<<tau<<" steps, ie "<<tau*parts.dt<<"s"<<endl;



		cout<<"Velocities ... ";
		if(ifstream((velSerie%0).c_str()).good() && ifstream((velSerie%(size-1)).c_str()).good())
			cout<<"have already been calculated"<<endl;
		else
		{
			cout<<"calculate"<<endl;
			#pragma omp parallel for schedule(runtime) shared(parts, tau, velSerie)
			for(ssize_t t=0; t<parts.getNbTimeSteps(); ++t)
			{
				vector<Coord> vel = parts.velocities(t, (tau+1)/2);
				ofstream v_f((velSerie%t).c_str(), ios::out | ios::trunc);
				v_f << VectorField(vel.begin(), vel.end(), "V");
				v_f.close();
			}
		}

		cout<<"Lost Neighbours ..."<<endl;
		if(ifstream((lngbSerie%0).c_str()).good() && ifstream((lngbSerie%(size-1)).c_str()).good())
			cout<<"have already been calculated"<<endl;
		else
		{
			cout<<"calculate"<<endl;
			export_lostNgb(parts.trajectories, tau, bondSerie, lngbSerie);
		}


		cout<<"Time averaged bond orientational order ... ";
		if(ifstream((timeBooSerie%0).c_str()).good() && ifstream((timeBooSerie%(size-1)).c_str()).good())
			cout<<"have already been calculated"<<endl;
		else
		{
			cout<<"calculate"<<endl;
			export_timeBoo(parts.trajectories, tau, timeBooSerie);
		}

		cout<<"Time averaged coarse grained bond orientational order ... ";
		if(ifstream((timecgBooSerie%0).c_str()).good() && ifstream((timecgBooSerie%(size-1)).c_str()).good())
			cout<<"have already been calculated"<<endl;
		else
		{
			cout<<"calculate"<<endl;
			export_timeBoo(parts.trajectories, tau, timecgBooSerie, "cg");
		}
		cout<<"Voronoi cell volume ... ";
		bool haveVolume = ifstream((timeVolSerie%0).c_str()).good() && ifstream((timeVolSerie%(size-1)).c_str()).good();
		if(haveVolume)
			cout<<"have already been calculated"<<endl;
		else
		{
			haveVolume = ifstream((volSerie%0).c_str()).good() && ifstream((volSerie%(size-1)).c_str()).good();
			if(haveVolume)
			{
				cout<<"calculate"<<endl;
				export_volume(parts.trajectories, tau, volSerie, timeVolSerie);
			}
			else
				cout<<"are not present"<<endl;
		}

		cout<<"export VTK"<<endl;
		for(size_t t=0; t<size; ++t)
		{
			//load bonds
			BondSet bonds = loadBonds(bondSerie%t);
			//load boo
			boost::multi_array<double, 2> qw, cg_qw;
			parts.positions[t].loadBoo(timeBooSerie%t, qw);
			parts.positions[t].loadBoo(timecgBooSerie%t, cg_qw);

			//load lost neighbours
			vector<double> lngb(parts.positions[t].size());
			{
				ifstream in((lngbSerie%t).c_str());
				copy(
					istream_iterator<size_t>(in), istream_iterator<size_t>(),
					lngb.begin()
					);
				in.close();
			}

			//load volumes
			vector<double> volumes(parts.positions[t].size());
			if(haveVolume)
			{
				ifstream in((timeVolSerie%t).c_str());
				string trash;
				getline(in, trash);
				getline(in, trash);
				copy(
					istream_iterator<double>(in), istream_iterator<double>(),
					volumes.begin()
					);
				in.close();
			}

			//load velocities
			vector<Coord> vel(parts.positions[t].size(), Coord(3));
			{
				ifstream in((velSerie%t).c_str());
				string trash;
				getline(in, trash);
				for(vector<Coord>::iterator p=vel.begin(); p!=vel.end(); ++p)
						in>>(*p);
				in.close();
			}

			vector<ScalarField> scalars(haveVolume?10:9, ScalarField(lngb.begin(), lngb.end(), "lostNgb"));
			for(int i=0;i<8;++i)
			{
				boost::multi_array<double, 2> &data = (i/4)?cg_qw:qw;
				scalars[i] = ScalarField(
					data.begin(), data.end(),
					string((i/4)?"cg":"")+string((i/2)%2?"W":"Q")+string(i%2?"6":"4"),
					i%4
					);
			}
			if(haveVolume)
				scalars.back() = ScalarField(volumes.begin(), volumes.end(), "Volume");

			vector<VectorField> vectors(1, VectorField(vel.begin(), vel.end(), "V"));

			parts.positions[t].exportToVTK(vtkSerie%t, bonds, scalars, vectors);
		}
		cout<<endl;
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

