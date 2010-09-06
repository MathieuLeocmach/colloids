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
#include "periodic.hpp"

using namespace std;
using namespace Colloids;

int main(int argc, char ** argv)
{
	try
    {
		if(argc<2) throw invalid_argument("Syntax : coordinateFile [maxBondLength]");

		const string filename(argv[1]);
		const string inputPath = filename.substr(0,filename.find_last_of("."));
		const string ext = filename.substr(filename.find_last_of(".")+1);
		double maxBondLength = 0.0;
		BondSet bonds;
		Particles parts(filename,1);
		parts.makeRTreeIndex();
		if(argc>2)
			maxBondLength = atof(argv[2]);
		else
		{
			vector<double> g = parts.getRdf(200,15.0);
			//set the max bond length as the first minima of g(r)
			//the loop is here only to get rid of possible multiple centers at small r
			vector<double>::iterator first_peak = g.begin();
			size_t first_min;
			do
			{
				first_peak = max_element(g.begin(),g.end());
				first_min = distance(g.begin(), min_element(first_peak,g.end()));
			}
			while(g[first_min]==0.0);
		}
		parts.makeNgbList(maxBondLength);
		bonds = parts.getBonds();
		ofstream output((inputPath + ".bonds").c_str(), ios::out | ios::trunc);
		for(BondSet::const_iterator b=bonds.begin(); b!= bonds.end();++b)
			output<<b->low()<<" "<<b->high()<<"\n";
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

