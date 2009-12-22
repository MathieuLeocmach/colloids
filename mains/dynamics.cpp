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

int main(int argc, char ** argv)
{
    if(argc<3)
    {
        cout << "compute both Mean Square displacement and Self fIntermediate scattering function"<<endl;
        cout << "Syntax : dynamic [path]filename mode" << endl;
        cout<<"\tmode=0\t No drift removal (default)"<<endl;
        cout<<"\tmode=1\t Drift removed"<<endl;
        cout<<"\tmode=2\t 0 then 1"<<endl;
        return EXIT_FAILURE;
    }

    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    size_t mode =0;
    if(argc>2)
        sscanf(argv[2],"%u",&mode);

    try
    {
        DynamicParticles parts(filename);
        if(mode==0 || mode==2)
        {
            cout <<"No drift removal"<<endl;
            parts.exportDynamics(inputPath+"_drift");
        }
        if(mode==1 || mode==2)
        {
            cout <<"Removing drift ... ";
            parts.removeDrift();
            cout<<"ok"<<endl;
            parts.exportDynamics(inputPath);
        }
    }
    catch(const std::exception &e)
    {
        cerr<<e.what() << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}




