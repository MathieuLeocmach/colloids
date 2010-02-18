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
#include <boost/format.hpp>
#include "../dynamicParticles.hpp"


using namespace std;
using namespace Colloids;

int main(int argc, char ** argv)
{
    if(argc<4)
    {
        cout << "compute Self Intermediate scattering function for sub-time intervals"<<endl;
        cout << "Syntax : ageing [path]filename start1 stop1 av1 [start2 stop2 av2 [...]]" << endl;
        return EXIT_FAILURE;
    }

    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    const size_t nbSub = (argc-2)/3;
    size_t start, stop, av;

    try
    {
        DynamicParticles parts(filename);
        parts.removeDrift();
        vector<double> ISF;
        boost::format name (inputPath+"_%1%from_%2%to_%3%av.isf");
        for(size_t i=0;i<nbSub;++i)
        {
        	sscanf(argv[3*i+2],"%u",&start);
        	sscanf(argv[3*i+3],"%u",&stop);
        	sscanf(argv[3*i+4],"%u",&av);
        	if(start+1>parts.getNbTimeSteps() || stop+av+1>parts.getNbTimeSteps())
				throw invalid_argument
				(
					(boost::format("[%1%,%2%] not included in [0,%3%]") % start % (stop+av) % (parts.getNbTimeSteps()-1)).str()
				);
			cout<<"["<<start<<","<<stop<<"] <"<<av<<">" << endl;
        	ISF = parts.getSelfISF(start,stop,av);

        	//export to file
        	ofstream output((name % start %stop %av).str().c_str(), std::ios::out | std::ios::trunc);
        	output <<"#t\tISF"<<endl;
        	for(size_t t=0;t<ISF.size();++t)
				output <<t*parts.dt<<"\t"<<ISF[t]<<"\n";
        	output.close();
        }
	}
    catch(const std::exception &e)
    {
        cerr<<e.what() << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
