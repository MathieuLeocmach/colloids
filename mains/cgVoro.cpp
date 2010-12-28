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
#include <voro++.cc>
#include <boost/progress.hpp>

using namespace std;
using namespace Colloids;

//Class extending voro++ container, in order to use protected data
#ifdef use_periodic
class VoroContainer : public container_poly
#else
class VoroContainer : public container
#endif
{
    public:
		vector< vector<int> > neighbours;
		#ifdef use_periodic
        VoroContainer(Particles &parts, vector<double> &radii, bool periodic) :
            container_poly(
        #else
        VoroContainer(Particles &parts, bool periodic) :
            container(
        #endif
                parts.bb.edges[0].first,parts.bb.edges[0].second,
                parts.bb.edges[1].first,parts.bb.edges[1].second,
                parts.bb.edges[2].first,parts.bb.edges[2].second,
                (parts.bb.edges[0].second - parts.bb.edges[0].first)*0.9/parts.radius + 1,
                (parts.bb.edges[1].second - parts.bb.edges[1].first)*0.9/parts.radius + 1,
                (parts.bb.edges[2].second - parts.bb.edges[2].first)*0.9/parts.radius + 1,
                periodic, periodic, periodic, 8
                )
        {
            for(size_t p=0; p<parts.size(); ++p)
                #ifdef use_periodic
                this->put(p, parts[p][0], parts[p][1], parts[p][2], radii[p]);
                #else
                this->put(p, parts[p][0], parts[p][1], parts[p][2]);
                #endif
        };

        void get_volumes(vector<double> &volumes)
        {
            neighbours.resize(volumes.size());
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
        };

        vector<double> get_cgVolumes(vector<double> &cgVolumes)
        {
            vector<double> volumes(cgVolumes.size());
            this->get_volumes(volumes);
            //coarse grain the volumes
            vector<size_t> div(cgVolumes.size(),0);
            for(size_t p=0; p<cgVolumes.size(); ++p)
            {
                //cout<<"p="<<p<<endl;
                //the center cell volume is significant only if not near a wall
                if(neighbours[p].empty() || neighbours[p].front()<0) continue;
                //cout<<"\tnot wall"<<endl;
                cgVolumes[p] += volumes[p];
                div[p]++;
                //cout<<"\tadded"<<endl;
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
                    if(!neighbours[p].empty())
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
            return volumes;
        }

        /**	\brief get the index of the particles neighbouring any wall and the index of the particles neighbouring them */
        void getOutsides(list<size_t> &outside, list<size_t> &secondOutside)
        {
        	outside.clear();
        	secondOutside.clear();
        	for(size_t p=0; p<neighbours.size(); ++p)
        		if(neighbours[p].front()<0)
        		{
        			outside.push_back(p);
        			secondOutside.push_back(p);
        			copy(
						lower_bound(neighbours[p].begin(), neighbours[p].end(), 0),
						neighbours[p].end(),
						back_inserter(secondOutside)
						);
        		}
        	secondOutside.sort();
        	secondOutside.unique();
        }

        /** print the bonds between centers (not walls) to a stream */
        void print_bonds(ostream &out)
        {
        	for(size_t p=0; p<neighbours.size(); ++p)
				for(
					vector<int>::const_iterator q = lower_bound(neighbours[p].begin(), neighbours[p].end(), p);
					q!=neighbours[p].end();	++q)
					out<<p<<" "<< *q <<"\n";
        }
};


int main(int argc, char ** argv)
{
    #ifdef use_periodic
	if(argc<7)
	{
		cerr<<"Syntax : cgVoro [path]filename.grv radiiFile Nb Dx Dy Dz [token size [offset]]" << endl;
		return EXIT_FAILURE;
	}
    #else
    if(argc<2)
	{
		cerr<<"Syntax : cgVoro [path]filename[.dat [token size [offset]] |.traj]" << endl;
		return EXIT_FAILURE;
	}
    #endif

	try
    {
		const string filename(argv[1]);
		const string ext = filename.substr(filename.find_last_of("."));
		string inputPath = filename.substr(0,filename.find_last_of("."));
		const string path = filename.substr(0, filename.find_last_of("/\\")+1);

		#ifdef use_periodic
		const size_t Nb = atoi(argv[3]);
		vector<double> radii(Nb);
		//load radii from file
		{
		    ifstream radiiFile(argv[2]);
		    if(!radiiFile.good())
                throw invalid_argument("No such file as "+string(argv[2]));
            copy(istream_iterator<double>(radiiFile), istream_iterator<double>(), radii.begin());
		}
		BoundingBox b;
		for(size_t d=0;d<3;++d)
		{
			b.edges[d].first=0.0;
			b.edges[d].second = atof(argv[4+d]);
		}
		if(argc>8)
		#else
		if(argc>2 || ext==".traj")
		#endif
		{
			double radius = 5.0;
			size_t offset, size;
			#ifdef use_periodic
			size = atol(argv[8]);
			offset = (argc<10)?0:atol(argv[9]);
            FileSerie datSerie(filename, string(argv[7]), size, offset);
            #else
            string pattern, token;
            if(ext==".traj")
			{
				double dt;
				ifstream trajfile(filename.c_str(), ios::in);
                if(!trajfile.good())
                    throw invalid_argument((filename+" doesn't exist").c_str() );
                trajfile >> radius >> dt;
                trajfile.ignore(1); //escape the endl
                getline(trajfile, pattern); //pattern is on the 2nd line
                getline(trajfile, token); //token is on the 3rd line
                trajfile >> offset >> size;
                pattern.insert(0, path);
			}
			else
			{
				pattern = filename;
 				token = string(argv[2]);
				size = atol(argv[3]);
				offset = (argc<5)?0:atol(argv[4]);
			}
			FileSerie datSerie(pattern, token, size, offset),
				outsideSerie = datSerie.changeExt(".outside"),
				secondOutsideSerie = datSerie.changeExt(".outside2");
			#endif
			FileSerie volSerie = datSerie.changeExt(".vol"),
                cgVolSerie = datSerie.addPostfix("_space", ".vol"),
				bondSerie = datSerie.addPostfix("_voro", ".bonds");

            boost::progress_display *showProgress;
            #pragma omp parallel shared(radius, showProgress) firstprivate(datSerie, volSerie)
            {
                #pragma omp single
                showProgress = new boost::progress_display(size);
                #pragma omp for
                for(size_t t=0; t<size;++t)
                {
                    string datfile;
                    #pragma omp critical
                    {
                        datfile = datSerie%t;
                    }
                	#ifdef use_periodic
					PeriodicParticles parts(Nb,b,datfile,radius);
					VoroContainer con(parts, radii, true);
					#else
                    Particles parts(datfile,radius);
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
                    vector<double> cgVolumes(parts.size(),0.0),
                        volumes = con.get_cgVolumes(cgVolumes);

                    //export
                    ofstream out((volSerie%t).c_str(), ios::out | ios::trunc);
                    copy(
                        volumes.begin(), volumes.end(),
                        ostream_iterator<double>(out, "\n")
                        );
					out.close();
					ofstream outcg((cgVolSerie%t).c_str(), ios::out | ios::trunc);
                    copy(
                        cgVolumes.begin(), cgVolumes.end(),
                        ostream_iterator<double>(outcg, "\n")
                        );
					outcg.close();

					ofstream bondfile((bondSerie%t).c_str(), ios::out | ios::trunc);
					con.print_bonds(bondfile);

					#ifndef use_periodic
					list<size_t> outside, secondOutside;
					con.getOutsides(outside, secondOutside);
					ofstream outsidefile((outsideSerie%t).c_str(), ios::out | ios::trunc);
					copy(outside.begin(), outside.end(), ostream_iterator<size_t>(outsidefile,"\n"));
					outsidefile.close();
					ofstream secondOutsidefile((secondOutsideSerie%t).c_str(), ios::out | ios::trunc);
					copy(secondOutside.begin(), secondOutside.end(), ostream_iterator<size_t>(secondOutsidefile,"\n"));
					secondOutsidefile.close();
					#endif

                    ++(*showProgress);
                }
            }
		}
		else
		{
            #ifdef use_periodic
            PeriodicParticles parts(Nb,b,filename,5.0);
            VoroContainer con(parts, radii, true);
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
            csv.close();

			ofstream bondfile((inputPath+".bonds").c_str(), ios::out | ios::trunc);
            con.print_bonds(bondfile);

            #ifndef use_periodic
            list<size_t> outside, secondOutside;
            con.getOutsides(outside, secondOutside);
            ofstream outsidefile((inputPath+".outside").c_str(), ios::out | ios::trunc);
            copy(outside.begin(), outside.end(), ostream_iterator<size_t>(outsidefile,"\n"));
            outsidefile.close();
            ofstream secondOutsidefile((inputPath+".outside2").c_str(), ios::out | ios::trunc);
            copy(secondOutside.begin(), secondOutside.end(), ostream_iterator<size_t>(secondOutsidefile,"\n"));
            secondOutsidefile.close();

            const string outName = inputPath.insert(inputPath.rfind("_t"), "_space");
            #else
            const string outName = inputPath+"_space";
            #endif
            ofstream out((outName+".vol").c_str(), ios::out | ios::trunc);
            copy(
                cgVolumes.begin(), cgVolumes.end(),
                ostream_iterator<double>(out, "\n")
                );
			out.close();
		}

    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
