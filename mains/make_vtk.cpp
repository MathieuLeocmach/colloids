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
    if(argc<3)
    {
        cout << "Syntax : make_vtk positions.dat radius" << endl;
        return EXIT_FAILURE;
    }

    try
    {
		const string filename(argv[1]);
		const string inputPath = filename.substr(0,filename.find_last_of("."));
		double radius;
		sscanf(argv[2],"%lf",&radius);

		IndexedParticles parts(filename,radius);

		std::vector< scalarField > scalars(4);
		scalars[0].first = "Q4";
		scalars[1].first = "Q6";
		scalars[2].first = "cgQ4";
		scalars[3].first = "cgQ6";

		parts.getBooFromFile(inputPath+".cloud",scalars[0].second,scalars[1].second);
		parts.getBooFromFile(inputPath+"_space.cloud",scalars[2].second,scalars[3].second);
		parts.exportToVTK(inputPath+".vtk",scalars,string("Bond Orientational Order"));
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
