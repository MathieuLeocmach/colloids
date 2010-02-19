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
        cout << "OR : cutter [path]filename token radius t_span t_offset minSep" << endl;
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
            radius = atof(argv[3]);
            const size_t t_span = atoi(argv[4]),
					t_offset = atoi(argv[5]);
            minsep = atof(argv[6]);
            const double sep = 2.0*radius*minSep;

            FileSerie datSerie(filename, token, t_span, t_offset);
            boost::progress_display show_progress(t_span);
            for(size_t t=0; t<t_span; ++t)
            {
                Particles(datSerie%t).cut(sep).exportToFile(datSerie%t);
                ++show_progress;
            }
        }
        else
        {
            if(argc<4)
            {
                cout << "Syntax : cutter [path]filename radius minSep" << endl;
                cout << "OR : cutter [path]filename token radius t_span t_offset minSep" << endl;
                cout << " minSep is in diameter unit" << endl;
                return EXIT_FAILURE;
            }
            radius = atof(argv[2]);
            minsep = atof(argv[3]);

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
