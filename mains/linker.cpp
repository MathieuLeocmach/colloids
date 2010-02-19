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

using namespace std;
using namespace Colloids;

int main(int argc, char ** argv)
{
    if(argc<6)
    {
        cout << "Syntax : linker [path]filename token radius time_step t_span t_offset(0)" << endl;
        return EXIT_FAILURE;
    }

    const string filename(argv[1]), token(argv[2]);
    const double radius = atof(argv[3]),
			time_step= atof(argv[4]);
    const size_t t_span = atoi(argv[5]),
		t_offset = (argc<7)?0:atoi(argv[6]);

    try
    {
    	FileSerie files(filename, token, t_span, t_offset);
        DynamicParticles parts(files, radius, time_step);
        cout << "total of " << parts.trajectories.size() << " trajectories" << endl;

        const string head = filename.substr(0,filename.rfind(token));
        const string nofolder = filename.substr(filename.find_last_of("/\\")+1);
        parts.save(head+".traj",nofolder,token,t_offset,t_span);

    }
    catch(const std::exception &e)
    {
        cerr <<e.what()<<endl;;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

