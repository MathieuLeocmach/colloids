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

//Define the preprocessor variable "use_periodic" if you want periodic boundary conditions
#include "../periodic.hpp"
#include "../files_series.hpp"

using namespace std;

int main(int argc, char ** argv)
{
	try
    {
		if(argc<5)
		{
			cout << "Radial Distribution function" << endl;
			cout << "Syntax : totalRdf [path]filename _t NbOfBins range" << endl;
			cout << " range is in pixel unit." << endl;
			return EXIT_FAILURE;
		}

		const string filename(argv[1]), token(argv[2]);
		const string inputPath = filename.substr(0,filename.rfind(token));
		double range;
		size_t Nbins;
		sscanf(argv[3],"%u",&Nbins);
		sscanf(argv[4],"%lf",&range);

		vector<double> gTot(Nbins,0.0);
		vector<string> tokens(1,token);
		TokenTree tt(tokens,filename);
		size_t t=0;

		try
		{
			while(true)
			{
				IndexedParticles parts((tt%t++).str(),1.0);
				vector<double> g = parts.getRdf(Nbins,range);
				transform(g.begin(),g.end(),gTot.begin(),gTot.begin(),plus<double>());
			}
		}
		catch(invalid_argument &e){};
		transform(gTot.begin(),gTot.end(),gTot.begin(),bind2nd(divides<double>(),(double)t));
		saveRDF(gTot,inputPath + "_total.rdf",Nbins/range);
		vector<double>::iterator firstPeak = max_element(gTot.begin(),gTot.end());
		vector<double>::iterator firstMin = min_element(firstPeak,gTot.end());
		cout<<(firstPeak-gTot.begin())*range/Nbins<<"\t"<<2.0*(firstPeak-gTot.begin())/(double)(firstMin-gTot.begin())<<endl;
	}
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
