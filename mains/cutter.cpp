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

//ask for the definition of the class IndexedParticles
#include "../particles.hpp"
#include "../files_series.hpp"
#include <boost/progress.hpp>

using namespace std;
using namespace Colloids;

int main(int argc, char ** argv)
{
    if(argc<4)
    {
        cout << "Syntax : cutter [path]filename radius minSep" << endl;
        cout << "OR : cutter [path]filename token radius t_offset t_span minSep" << endl;
        cout << " minSep is in diameter unit" << endl;
        return EXIT_FAILURE;
    }

    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    const string ext = filename.substr(filename.find_last_of(".")+1);
    double radius,minSep;

    try
    {
        if(argc>4)
        {
            cout<<"file serie"<<endl;
            const string token(argv[2]);
            sscanf(argv[3],"%lf",&radius);
            size_t t_offset,t_span;
            sscanf(argv[4],"%u",&t_offset);
            sscanf(argv[5],"%u",&t_span);
            sscanf(argv[6],"%lf",&minSep);
            const double sep = 2.0*radius*minSep;

            vector<string> tokens(1,token);
            TokenTree tt(tokens,filename);
            vector<size_t> v(1,t_offset);
            boost::progress_display show_progress(t_span);
            while(v[0]<t_offset+t_span)
            {
                Particles(tt(v)).cut(sep).exportToFile(tt(v));
                v[0]++;
                ++show_progress;
            }
        }
        else
        {
            if(argc<4)
            {
                cout << "Syntax : cutter [path]filename radius minSep" << endl;
                cout << "OR : cutter [path]filename token radius time_step t_offset t_span minSep" << endl;
                cout << " minSep is in diameter unit" << endl;
                return EXIT_FAILURE;
            }
            sscanf(argv[2],"%lf",&radius);
            sscanf(argv[3],"%lf",&minSep);

            const double sep = 2.0*radius*minSep;
            Particles(filename).cut(sep).exportToFile(filename);
        }
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
