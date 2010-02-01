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

