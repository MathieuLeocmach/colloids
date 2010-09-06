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
		if(argc<3) throw invalid_argument("Syntax : dat2vtk coordinateFile");

		const string filename(argv[1]);
		const string inputPath = filename.substr(0,filename.find_last_of("."));
		const string ext = filename.substr(filename.find_last_of(".")+1);
		const size_t tokenPos = inputPath.rfind("_t");

		Particles parts(filename);
		BondSet bonds = loadBonds(inputPath+".bonds");

		boost::multi_array<double, 2> qw, Sqw;
		parts.loadBoo(inputPath+".cloud", qw);

		parts.loadBoo(
			inputPath.substr(0, tokenPos)+"_space"+inputPath.substr(tokenPos)+".cloud",
			Sqw);
		vector<ScalarField> scalars;
		scalars.reserve(8);
		for(size_t i=0;i<4;++i)
			scalars.push_back(ScalarField(qw.begin(), qw.end(), "", i));
		for(size_t i=0;i<4;++i)
			scalars.push_back(ScalarField(Sqw.begin(), Sqw.end(), "", i));
		for(size_t i=0;i<8;++i)
			scalars[i].name = string(i/4?"cg":"")+string((i/2)%2?"W":"Q")+string(i%2?"6":"4");

		parts.exportToVTK(inputPath+".vtk", bonds, scalars, vector<VectorField>());
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
