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
	    vector<double> lngb(trajectories.inverse[t0].size(),0);
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
                lngb[p] = lost.size();// * (1+tau*2)/(double)(stop-start);
            }
	    }
		ofstream f((lngbSerie%t0).c_str(), ios::out | ios::trunc);
		copy(
			lngb.begin(), lngb.end(),
			ostream_iterator<double>(f,"\n")
			);
		f.close();
	}
}

void export_timeBoo(const TrajIndex& trajectories, const size_t& tau, FileSerie &timeBooSerie, const string &prefix="")
{
	const size_t size = trajectories.inverse.size();
	FileSerie booSerie = timeBooSerie.changeExt(".cloud");
	boost::multi_array<double, 2> qw;
	vector<ScalarDynamicField> scalars(8, ScalarDynamicField(trajectories, tau));
	for(size_t i=0;i<scalars.size();++i)
	{
	    stringstream o(prefix);
	    o<<((i/4)?"W":"Q")<<((i%4)*2+4);
		scalars[i].name = o.str();
	}
	for(size_t t=0; t<size; ++t)
	{
		//read qw from file
		boost::array<size_t, 2> shape = {{trajectories.inverse[t].size(), scalars.size()}};
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
		vector<ScalarField> s(scalars.size(), ScalarField("", trajectories.inverse[t].size()));
		for(int i=0; i<s.size(); ++i)
			s[i] = scalars[i][t];
		ofstream f((timeBooSerie%t).c_str(), ios::out | ios::trunc);
		f<<"#Q4\tQ6\tQ8\tQ10\tW4\tW6\tW8\tW10"<<endl;
		for(size_t p=0; p<s.front().size(); ++p)
		{
			for(size_t i=0;i<s.size();++i)
				f << s[i][p]<<"\t";
			f<<"\n";
		}
		f.close();
	}
}

void export_phi(const TrajIndex& trajectories, const double &radius, const size_t& tau, FileSerie &volSerie, FileSerie &phiSerie)
{
	const size_t size = trajectories.inverse.size();
	const double unitVolume = 4/3*pow(radius, 3.0)*M_PI;
	ScalarDynamicField phi(trajectories, tau, "phi");
	for(size_t t=0; t<size; ++t)
	{
		vector<double> vol(trajectories.inverse[t].size());
		ifstream f((volSerie%t).c_str(), ios::in);
		copy(
			istream_iterator<double>(f), istream_iterator<double>(),
			vol.begin()
			);
		for(size_t p=0; p<vol.size(); ++p)
			if(1.0 + vol[p]*vol[p] > 1.0)
				vol[p] = unitVolume / vol[p];
		phi.push_back(ScalarField(vol.begin(), vol.end()));
	}
	//export
	for(size_t t=0; t<size; ++t)
	{
		ofstream f((phiSerie%t).c_str(), ios::out | ios::trunc);
		f<<phi[t];
		f.close();
	}
}

size_t loadTau(const string &filename, const size_t &size)
{
	cout<<"load isf to find relaxation time ... ";
	vector< vector<double> > isf(4,vector<double>(size));
	ifstream in(filename.c_str());
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

	//Find relaxation time
	if(isf.back().back() > exp(-1.0))
		return size-1;
	else
		return size-1 - (upper_bound(isf.back().rbegin(),isf.back().rend(),exp(-1.0))-isf.back().rbegin());
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
    	//construct the trajectory index anyway
    	TrajIndex trajectories;
    	double radius, dt;
		string pattern, token;
		size_t offset, size, tau;
		{
			ifstream trajfile(filename.c_str(), ios::in);
			if(!trajfile.good())
				throw invalid_argument((filename+" doesn't exist").c_str() );
			trajfile >> radius >> dt;
			trajfile.ignore(1); //escape the endl
			getline(trajfile, pattern); //pattern is on the 2nd line
			getline(trajfile, token); //token is on the 3rd line
			trajfile >> offset >> size;
			trajfile >> trajectories;
			trajfile.close();
			trajectories.makeInverse(trajectories.getFrameSizes(size));
		}
		cout << trajectories.size() << " particles in "<<size<<" time steps"<<endl;
		//File series
		FileSerie datSerie(path+pattern, token, size, offset),
			velSerie = datSerie.changeExt(".vel"),
			bondSerie = datSerie.changeExt(".bonds"),
			timeBooSerie = datSerie.changeExt(".boo"),
			timecgBooSerie = datSerie.addPostfix("_space", ".boo"),
			timesfBooSerie = datSerie.addPostfix("_surf", ".boo"),
			lngbSerie = datSerie.changeExt(".lngb"),
			volSerie = datSerie.addPostfix("_space", ".vol"),
			phiSerie = datSerie.addPostfix("_space", ".phi"),
			vtkSerie = datSerie.addPostfix("_dynamic", ".vtk");

		const bool hasVel = ifstream((velSerie%0).c_str()).good() && ifstream((velSerie%(size-1)).c_str()).good(),
			hasISF = ifstream((inputPath + ".isf").c_str()).good();
		if((argc<3 && !hasISF) || !hasVel)
		{
			cout<<"load positions"<<endl;
			DynamicParticles parts(filename);
			cout<<"remove drift ... ";
			parts.removeDrift();

			if(argc<3 && !hasISF)
			{
				cout<<"calculate ISF"<<endl;
				//Dynamics calculation and export
				vector< vector<double> > isf(4,vector<double>(size));
				vector<double> msd, ngp;
				parts.makeDynamics(msd, isf, ngp);

				ofstream msd_f((inputPath + ".msd").c_str());
				msd_f << "#t\tMSD"<<endl;
				for(size_t t=0; t<msd.size(); ++t)
					msd_f << t*parts.dt<<"\t"<< msd[t]<<"\n";
				msd_f.close();

				ofstream ngp_f((inputPath + ".ngp").c_str());
				ngp_f << "#t\tNGP"<<endl;
				for(size_t t=0; t<ngp.size(); ++t)
					ngp_f << t*parts.dt<<"\t"<< ngp[t]<<"\n";
				ngp_f.close();

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

			if(!hasVel)
			{
				if(argc<3)
					tau = loadTau(inputPath + ".isf", size);
				else
					tau = atoi(argv[2]);
				cout<<"calculate velocities"<<endl;
				#pragma omp parallel for schedule(runtime) shared(parts, tau, velSerie)
				for(ssize_t t=0; t<parts.getNbTimeSteps(); ++t)
				{
					vector<Coord> vel = parts.velocities(t, (tau+1)/2);
					ofstream v_f((velSerie%t).c_str(), ios::out | ios::trunc);
					v_f << VectorField(vel.begin(), vel.end(), "V");
					v_f.close();
				}
			}
		} //no more positions in memory

		if(hasVel)
		{
			if(argc<3)
				tau = loadTau(inputPath + ".isf", size);
			else
				tau = atoi(argv[2]);
		}
		cout<<"relaxation time is "<<tau<<" steps, ie "<<tau*dt<<"s"<<endl;

		cout<<"Lost Neighbours ..."<<endl;
		if(ifstream((lngbSerie%0).c_str()).good() && ifstream((lngbSerie%(size-1)).c_str()).good())
			cout<<"have already been calculated"<<endl;
		else
		{
			cout<<"calculate"<<endl;
			export_lostNgb(trajectories, tau, bondSerie, lngbSerie);
		}

		cout<<"Time averaged bond orientational order ... ";
		if(ifstream((timeBooSerie%0).c_str()).good() && ifstream((timeBooSerie%(size-1)).c_str()).good())
			cout<<"have already been calculated"<<endl;
		else
		{
			cout<<"calculate"<<endl;
			export_timeBoo(trajectories, tau, timeBooSerie);
		}

		cout<<"Time averaged coarse grained bond orientational order ... ";
		if(ifstream((timecgBooSerie%0).c_str()).good() && ifstream((timecgBooSerie%(size-1)).c_str()).good())
			cout<<"have already been calculated"<<endl;
		else
		{
			cout<<"calculate"<<endl;
			export_timeBoo(trajectories, tau, timecgBooSerie, "cg");
		}

		/*cout<<"Time averaged surface bond orientational order ... ";
		if(ifstream((timesfBooSerie%0).c_str()).good() && ifstream((timesfBooSerie%(size-1)).c_str()).good())
			cout<<"have already been calculated"<<endl;
		else
		{
			cout<<"calculate"<<endl;
			export_timeBoo(trajectories, tau, timesfBooSerie, "sf");
		}*/

		cout<<"Voronoi cell volume fraction ... ";
		bool haveVolume = ifstream((phiSerie%0).c_str()).good() && ifstream((phiSerie%(size-1)).c_str()).good();
		if(haveVolume)
			cout<<"have already been calculated"<<endl;
		else
		{
			haveVolume = ifstream((volSerie%0).c_str()).good() && ifstream((volSerie%(size-1)).c_str()).good();
			if(haveVolume)
			{
				cout<<"calculate"<<endl;
				export_phi(trajectories, radius, tau, volSerie, phiSerie);
			}
			else
				cout<<"voronoi cell volume files are not present"<<endl;
		}

		cout<<"export VTK"<<endl;
		boost::progress_display show_progress(size);
		for(size_t t=0; t<size; ++t)
		{
			//load positions
			Particles parts(datSerie%t);
			//load bonds
			BondSet bonds = loadBonds(bondSerie%t);
			//load boo
			boost::multi_array<double, 2> cg_qw;
			parts.loadBoo(timecgBooSerie%t, cg_qw);

			//load lost neighbours
			vector<double> lngb(parts.size());
			{
				ifstream in((lngbSerie%t).c_str());
				copy(
					istream_iterator<double>(in), istream_iterator<double>(),
					lngb.begin()
					);
				in.close();
			}

			//load volumes
			vector<double> phi(parts.size());
			if(haveVolume)
			{
				ifstream in((phiSerie%t).c_str());
				string trash;
				getline(in, trash);
				getline(in, trash);
				copy(
					istream_iterator<double>(in), istream_iterator<double>(),
					phi.begin()
					);
				in.close();
			}

			//load velocities
			vector<Coord> vel(parts.size(), Coord(3));
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
			    stringstream o;
                o<<((i/4)?"W":"Q")<<((i%4)*2+4);
                scalars[i].name = o.str();
				scalars[i] = ScalarField(
					cg_qw.begin(), cg_qw.end(),
					o.str(),
					i%8
					);
			}
			if(haveVolume)
				scalars.back() = ScalarField(phi.begin(), phi.end(), "phi");

			vector<VectorField> vectors(1, VectorField(vel.begin(), vel.end(), "velocity"));

			parts.exportToVTK(vtkSerie%t, bonds, scalars, vectors);
			++show_progress;
		}
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

