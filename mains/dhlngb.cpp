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

#include "dynamicParticles.hpp"

#include <boost/progress.hpp>

using namespace std;
using namespace Colloids;

int main(int argc, char ** argv)
{
    if(argc<3)
    {
        cerr<<"Compute the lost neighbours between an initial time step and subsequent time steps"<<endl;
    	cerr<<"dhlngb [path]filename.traj initial" << endl;
		return EXIT_FAILURE;
    }
    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    const string ext = filename.substr(filename.find_last_of(".")+1);
    const string path = filename.substr(0, filename.find_last_of("/\\")+1);
    const size_t initial = atoi(argv[2]);

    try
    {
    	//construct the trajectory index
    	TrajIndex trajectories;
    	double radius, dt;
		string pattern, token;
		size_t offset, size;
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
		if(initial+1 >= size)
            throw invalid_argument("Initial time step after the end !");

        //File series
		FileSerie datSerie(path+pattern, token, size, offset),
			bondSerie = datSerie.changeExt(".bonds"),
			lngbSerie = datSerie.addPostfix("_"+string(argv[2])+"dynhet", ".lngb");

        typedef vector<size_t>	Ngbs;
        typedef boost::ptr_vector<Ngbs> NgbFrame;
        typedef boost::ptr_vector<NgbFrame> DynNgbs;
        typedef list<size_t> ListNgb;
        typedef vector<ListNgb> ListNgbFrame;

        //what are the trajectories neighbouring the position p at time t
        DynNgbs dyn(size-initial);
        for(size_t t=initial; t<size;++t)
        {
            dyn.push_back(new NgbFrame(trajectories.inverse[t].size()));
            for(size_t p=0; p<trajectories.inverse[t].size();++p)
                dyn[t].push_back(new Ngbs());
        }
        {
            boost::progress_timer ti;
            #pragma omp parallel for schedule(static)
            for(size_t t=initial; t<size;++t)
            {
                //Easily filled, disordered but easy to sort container. Bad memory perf. But it's just one frame.
                ListNgbFrame easy(dyn[t-size].size());
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
                    dyn[t-initial][p].assign(easy[p].begin(), easy[p].end());
                }
            }
            cout<<"input & assignement\t";
        }

        //calculate and export the neighbourhood difference between initial and t>initial
        for(size_t t=initial+1; t<size; ++t)
        {
            //trajectories starting after initial time step get a -1
            vector<int> lngb(trajectories.inverse[initial].size(),-1);
            #pragma omp parallel for shared(trajectories, lngb, t) schedule(dynamic)
            for(size_t p=0;p<trajectories.inverse[initial].size(); ++p)
            {
                const Traj &tr = trajectories[trajectories.inverse[initial][p]];
                if(tr.last_time()>=t)
                {
                    ListNgb lost;
                    //how many neighbours have been lost by trajectory tr (pth particle of initial time step) at time t ?
                    set_difference(
                        dyn[t-initial][tr[t-initial]].begin(), dyn[t-initial][tr[t-initial]].end(),
                        dyn[0][p].begin(), dyn[0][p].end(),
                        back_inserter(lost)
                        );
                    lngb[p] = lost.size();
                }
            }
            //export
            ofstream f((lngbSerie%t).c_str(), ios::out | ios::trunc);
            copy(
                lngb.begin(), lngb.end(),
                ostream_iterator<int>(f,"\n")
                );
            f.close();
        }
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
