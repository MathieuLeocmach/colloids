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
#include <voro++.cc>
#include <boost/progress.hpp>

using namespace std;
using namespace Colloids;

//Class extending voro++ container, in order to use protected data
class VoroContainer : public container
{
    public:
        VoroContainer(Particles &parts, bool periodic) :
        container(
            parts.bb.edges[0].first,parts.bb.edges[0].second,
            parts.bb.edges[1].first,parts.bb.edges[1].second,
            parts.bb.edges[2].first,parts.bb.edges[2].second,
            (parts.bb.edges[0].second - parts.bb.edges[0].first)*0.9/parts.radius + 1,
            (parts.bb.edges[1].second - parts.bb.edges[1].first)*0.9/parts.radius + 1,
            (parts.bb.edges[2].second - parts.bb.edges[2].first)*0.9/parts.radius + 1,
            periodic, periodic, periodic, 8
            )
        {
            for(int p=0; p<parts.size(); ++p)
                this->put(p, parts[p][0], parts[p][1], parts[p][2]);
        };

        void get_cgVolumes(vector<double> &cgVolumes)
        {
            vector<double> volumes(cgVolumes.size(),0.0);
            vector< vector<int> > neighbours(cgVolumes.size());
            {
                //compute voronoi diagram and export it in a memory buffer (only way get data from voro++)
                stringstream buffer;
                this->print_all_custom("%i %v %s %n", buffer);
                //read back the buffer
                int p, s;
                while(buffer.good())
                {
                    buffer >> p;
                    buffer >> volumes[p];
                    buffer >> s;
                    neighbours[p].resize(s);
                    for(int i=0; i<s; ++i)
                        buffer >> neighbours[p][i];
                }
            }
            //sort ngbs
            for(vector< vector<int> >::iterator p = neighbours.begin(); p!= neighbours.end(); ++p)
                sort(p->begin(), p->end());
            //coarse grain the volumes
            vector<size_t> div(cgVolumes.size(),0);
            for(size_t p=0; p<cgVolumes.size(); ++p)
            {
                //the center cell volume is significant only if not near a wall
                if(neighbours[p].front()<0) continue;
                cgVolumes[p] += volumes[p];
                div[p]++;
                for(vector<int>::const_iterator n = neighbours[p].begin(); n!=neighbours[p].end(); ++n)
                {
                    cgVolumes[*n] += volumes[p];
                    div[*n]++;
                }
            }
            deque<size_t> secondChance;
            for(size_t p=0; p<cgVolumes.size(); ++p)
                if(div[p]>0)
                    cgVolumes[p] /= div[p];
                else
                    secondChance.push_back(p);

            //second chance for the particles whose neighbours where all near wall
            while(!secondChance.empty())
            {
                const size_t p = secondChance.front();
                secondChance.pop_front();
                int nb = 0;
                for(vector<int>::const_iterator n = lower_bound(neighbours[p].begin(), neighbours[p].end(), 0); n!=neighbours[p].end(); ++n)
                    if(1.0+cgVolumes[*n]*cgVolumes[*n] > 1.0)
                    {
                        cgVolumes[p] += cgVolumes[*n];
                        nb++;
                    }
                if(nb>0)
                    cgVolumes[p]/=(double)nb;
                else
                    secondChance.push_back(p);
            }
        }
};


int main(int argc, char ** argv)
{
    #ifdef use_periodic
	if(argc<6)
	{
		cerr<<"Syntax : cgVoro [path]filename.grv Nb Dx Dy Dz" << endl;
		return EXIT_FAILURE;
	}
    #else
    if(argc<2)
	{
		cerr<<"Syntax : cgVoro [path]filename[.dat|.traj]" << endl;
		return EXIT_FAILURE;
	}
    #endif

	try
    {
		const string filename(argv[1]);
		const string ext = filename.substr(filename.find_last_of("."));
		string inputPath = filename.substr(0,filename.find_last_of("."));
		const string path = filename.substr(0, filename.find_last_of("/\\")+1);
		//no traj for periodic version
		#ifndef use_periodic
		if(ext==".traj")
		{
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
            }
            FileSerie datSerie(path+pattern, token, size, offset),
                volSerie = datSerie.addPostfix("_space", ".vol");

            boost::progress_display *showProgress;
            #pragma omp parallel shared(radius, showProgress) firstprivate(datSerie, volSerie)
            {
                #pragma omp single
                showProgress = new boost::progress_display(size);
                #pragma omp for
                for(size_t t=0; t<size;++t)
                {
                    Particles parts(datSerie%t,radius);
                    valarray<double> maxi = parts.front(), mini = parts.front();
                    for(Particles::const_iterator p= parts.begin(); p!=parts.end(); ++p)
                        for(int d=0; d<3;++d)
                        {
                            maxi[d] = max(maxi[d], (*p)[d]);
                            mini[d] = min(mini[d], (*p)[d]);
                        }
                    for(int d=0; d<3;++d)
                    {
                        parts.bb.edges[d].first = mini[d]-1;
                        parts.bb.edges[d].second = maxi[d]+1;
                    }
                    VoroContainer con(parts, false);
                    vector<double> cgVolumes(parts.size(),0.0);
                    con.get_cgVolumes(cgVolumes);

                    //export
                    ofstream out((volSerie%t).c_str(), ios::out | ios::trunc);
                    copy(
                        cgVolumes.begin(), cgVolumes.end(),
                        ostream_iterator<double>(out, "\n")
                        );
                    ++(*showProgress);
                }
            }
		}
		else
		#endif
		{
            #ifdef use_periodic
            const size_t Nb = atoi(argv[3]);
            BoundingBox b;
            for(size_t d=0;d<3;++d)
            {
                b.edges[d].first=0.0;
                b.edges[d].second = atof(argv[4+d]);
            }
            PeriodicParticles parts(Nb,b,filename,5.0);
            VoroContainer con(parts, true);
            #else
            Particles parts(filename,5.0);
            valarray<double> maxi = parts.front(), mini = parts.front();
            for(Particles::const_iterator p= parts.begin(); p!=parts.end(); ++p)
                for(int d=0; d<3;++d)
                {
                    maxi[d] = max(maxi[d], (*p)[d]);
                    mini[d] = min(mini[d], (*p)[d]);
                }
            for(int d=0; d<3;++d)
            {
                parts.bb.edges[d].first = mini[d]-1;
                parts.bb.edges[d].second = maxi[d]+1;
            }
            VoroContainer con(parts, false);
            #endif
            vector<double> cgVolumes(parts.size(),0.0);
            con.get_cgVolumes(cgVolumes);

            //export
            ofstream csv((inputPath+".csv").c_str(), ios::out | ios::trunc);
            csv<<"x,y,z,v\n";
            for(size_t p=0; p<parts.size(); ++p)
            {
                for(int d=0;d<3;++d)
                    csv<<parts[p][d]<<",";
                csv<<cgVolumes[p]<<"\n";
            }
            const string outName = inputPath.insert(inputPath.rfind("_t"), "_space");
            ofstream out((outName+".vol").c_str(), ios::out | ios::trunc);
            copy(
                cgVolumes.begin(), cgVolumes.end(),
                ostream_iterator<double>(out, "\n")
                );
		}

    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
