/**
    Copyright 2010 Mathieu Leocmach

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

using namespace std;

int main(int argc, char ** argv)
{
	try
    {
		if(argc<3) throw invalid_argument("Syntax : dat2vtk radius coordinateFile [neigboursFile [dataFile]]");

		const string filename(argv[1]);
		const string inputPath = filename.substr(0,filename.find_last_of("."));
		const string ext = filename.substr(filename.find_last_of(".")+1);
		double radius;
		sscanf(argv[2],"%lf",&radius);

		IndexedParticles parts(filename,radius);
		cout<<"calculate neighbour list"<<endl;
		/*const double shell =1.3;
		vector< set<size_t> >ngbList;
		parts.getNgbList(shell,ngbList);*/
		deque<pair<size_t, size_t> > bonds = parts.getBonds(1.3*2.0*radius);

		map<size_t,valarray<double> > qw,Sqw;
		cout<<"get raw Boo from files"<<endl;
		parts.getBooFromFile(inputPath+".cloud",qw);
		cout<<"get coarse grained Boo from files"<<endl;
		parts.getBooFromFile(inputPath+"_space.cloud",Sqw);

		vector<string> datanames(8);
		for(size_t i=0;i<datanames.size();++i)
			datanames[i] = string(i/4?"cg":"")+string((i/2)%2?"W":"Q")+string(i%2?"6":"4");

		ofstream output((inputPath+".vtk").c_str(), ios::out | ios::trunc);
		if(!output)
			throw invalid_argument("Cannot write on "+filename);

		output<<"# vtk DataFile Version 3.0\n"
				<<"particles\n"
				"ASCII\n"
				"DATASET POLYDATA\n"
				"POINTS "<<parts.size()<<" double\n";
		for(IndexedParticles::const_iterator p=parts.begin();p!=parts.end();++p)
		{
			for(size_t d=0;d<3;++d)
				output<<(*p)[d]<<" ";
			output<<"\n";
		}

		output << "LINES "<<bonds.size()<<" "<<bonds.size()*3<<endl;
		for(deque< pair<size_t,size_t> >::const_iterator b= bonds.begin();b!=bonds.end();++b)
			output<<"2 "<< (*b).first<<" "<<(*b).second<<endl;

		output<<"POINT_DATA "<<parts.size()<<endl;
		for(size_t s=0;s<datanames.size();++s)
		{
			output<<"SCALARS "<< datanames[s] <<" double\n"
					"LOOKUP_TABLE default\n";

			size_t p=0;
			for(map<size_t,valarray<double> >::const_iterator l = ((s/4)?Sqw:qw).begin();l!=((s/4)?Sqw:qw).end();++l)
			{
				while(p<(*l).first)
				{
					output<<0<<endl;
					p++;
				}
				p++;
				output<<(*l).second[s%4]<<endl;
			}
			while(p<parts.size())
			{
				output<<0<<endl;
				p++;
			}
		}
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
