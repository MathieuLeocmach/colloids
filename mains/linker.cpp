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


#include "../dynamicParticles.hpp"
//#include "../time_tracker.hpp"
//#include "../files_series.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    if(argc<7)
    {
        cout << "Syntax : linker [path]filename token radius time_step t_offset t_span " << endl;
        return EXIT_FAILURE;
    }

    const string filename(argv[1]), token(argv[2]);
    double radius,time_step;
    sscanf(argv[3],"%lf",&radius);
    sscanf(argv[4],"%lf",&time_step);
    size_t t_offset,t_span;
    sscanf(argv[5],"%u",&t_offset);
    sscanf(argv[6],"%u",&t_span);

    try
    {
        DynamicParticles parts(radius,time_step,filename,token,t_offset,t_span);
        cout << "total of " << parts.trajectories.size() << " trajectories" << endl;

        const string head = filename.substr(0,filename.rfind(token));
        const string noExt = filename.substr(0,filename.find_last_of("."));
        const string nofolder = filename.substr(filename.find_last_of("/\\")+1);
        const string folder = filename.substr(0,filename.find_last_of("/\\")+1);
        parts.save(head+".traj",nofolder,token,t_offset,t_span);

		//marks a few trajectories to verify the time tracking
		vector<map <size_t,unsigned char> > labels(t_span);
		set<size_t> spanning = parts.getSpanning(0,t_span-1);
		set<size_t>::iterator tr = spanning.begin();
		for(size_t i=0;i<5 && tr!=spanning.end();++i)
		{
			for(size_t t=0;t<t_span;++t)
				labels[t].insert(make_pair(parts.trajectories[*tr][t],(unsigned char)15));
			tr++;
		}
        parts.exportToPV(noExt+"_linked.pv",labels);
    }
    catch(const std::exception &e)
    {
        cerr <<e.what()<<endl;;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

