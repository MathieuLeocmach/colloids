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

TrajIndex export_post_lostNgb(const TrajIndex &trajectories, const size_t &tau, FileSerie &bondSerie, FileSerie &lngbSerie)
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
	{
	    boost::progress_timer ti;
	    #pragma omp parallel for schedule(static)
	for(size_t t=0; t<size;++t)
	{
		//Easily filled, disordered but easy to sort container. Bad memory perf. But it's just one frame.
		ListNgbFrame easy(dyn[t].size());
		//load bonds from file and bin their ends to neighbour list
		string bondfile;
		#pragma omp critical
		{
		    bondfile = bondSerie%t;
		}
		ifstream f(bondfile.c_str());
		size_t a,b;
		while(f.good())
		{
			f>>a>>b;
			easy[a].push_back(trajectories.inverse[t][b]);
			easy[b].push_back(trajectories.inverse[t][a]);
		}
		f.close();


		//#pragma omp parallel for shared(easy, dyn) schedule(dynamic)
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
    cout<<"input & assignement\t";
	}

	//look at the neighbourhood difference between t0 and t0+tau
    for(size_t t0=0; t0<size-tau; ++t0)
	{
	    vector<int> lngb(trajectories.inverse[t0].size(),-1);
	    #pragma omp parallel for shared(trajectories, lngb, t0) schedule(dynamic)
	    for(size_t p=0;p<trajectories.inverse[t0].size(); ++p)
	    {
	        const Traj &tr = trajectories[trajectories.inverse[t0][p]];
	        if(tr.last_time()>=t0+tau)
	        {
	        	ListNgb lost;
                set_difference(
                    dyn[t0+tau][tr[t0+tau]].begin(), dyn[t0+tau][tr[t0+tau]].end(),
                    dyn[t0][tr[t0]].begin(), dyn[t0][tr[t0]].end(),
                    back_inserter(lost)
                    );
                lngb[p] = lost.size();
	        }
	    }
		ofstream f((lngbSerie%t0).c_str(), ios::out | ios::trunc);
		copy(
			lngb.begin(), lngb.end(),
			ostream_iterator<int>(f,"\n")
			);
		f.close();
	}

	//bonus : split the trajectories in neighbourhoods (same neighbours all along)
	TrajIndex split;
	{
	    boost::progress_timer ti;
	for(TrajIndex::const_iterator tr=trajectories.begin(); tr!=trajectories.end(); ++tr)
    {
		split.push_back(Traj(tr->start_time, tr->steps.front()));
		Traj *s = &split.back();
		Ngbs *ref_ngb = &dyn[s->start_time][s->steps.front()];
		for(size_t t=tr->start_time+1; t<=tr->last_time(); ++t)
		{
			if(ref_ngb->size()==dyn[t][(*tr)[t]].size() &&
				equal(
					ref_ngb->begin(), ref_ngb->end(),
					dyn[t][(*tr)[t]].begin()
					))
			{
				//in case of unchanged neighbourhood, the trajectory is prolonged
				s->push_back((*tr)[t]);
			}
			else
			{
				//if the neighbourhood has changed, a new trajectory begins
				split.push_back(Traj(t, (*tr)[t]));
				s = &split.back();
				ref_ngb = &dyn[s->start_time][s->steps.front()];
			}
		}
	}
	cout<<"trajectory splitting\t";
	}
	return split;
}

int main(int argc, char ** argv)
{
    if(argc<3)
    {
    	cerr<<"lostngb [path]filename.traj AveragingInterval" << endl;
		return EXIT_FAILURE;
    }
    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    const string ext = filename.substr(filename.find_last_of(".")+1);
    const string path = filename.substr(0, filename.find_last_of("/\\")+1);

    try
    {
    	//construct the trajectory index
    	TrajIndex trajectories;
    	double radius, dt;
		string pattern, token;
		size_t offset, size, tau = atoi(argv[2]);
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
		cout << "relaxation time is "<<tau<<" steps, ie "<<tau*dt<<"s"<<endl;
		//File series
		FileSerie datSerie(path+pattern, token, size, offset),
			bondSerie = datSerie.changeExt(".bonds"),
			lngbSerie = datSerie.addPostfix("_post", ".lngb");

		cout<<"Lost Neighbours between t0 and t0+tau..."<<endl;
		TrajIndex split = export_post_lostNgb(trajectories, tau, bondSerie, lngbSerie);

		ofstream output((inputPath+"_split.traj").c_str(), ios::out | ios::trunc);
		if(output.good())
		{
			//header
			output << radius << "\t" << dt << endl;

			//link to the positions files serie
			output << pattern << endl;
			output << token << endl;
			output << offset << "\t" << size << endl;

			//trajectory data
			output << split;

			output.close();
		}
		else cerr << " cannot open the file";
	}
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
