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

#include "../indexedParticles.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    if(argc<2)
    {
        cout << "compute number density function of x,y and z"<<endl;
        cout << "Syntax : gradient [path]filename NbBins" << endl;
        return EXIT_FAILURE;
    }

    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    size_t NbBins;
    sscanf(argv[2],"%u",&NbBins);


    Particles parts(filename);

    vector< std::vector<size_t> > composition(3,vector<size_t>(NbBins,0));

    for(size_t d=0;d<3;++d)
    {
        const double length = (parts.bb.edges[d].second-parts.bb.edges[d].first)/NbBins;
        for(Particles::iterator p = parts.begin();p!=parts.end();++p)
            composition[d][(size_t)((*p)[d]/length)]++;
    }

    saveTable(composition.begin(),composition.end(),inputPath+".gradient","px\tx\ty\tz");
    return EXIT_SUCCESS;
}
