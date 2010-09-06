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
*/

#include "../lifTracker.hpp"

using namespace std;
using namespace Colloids;

#ifndef INSTAL_PATH
#define INSTAL_PATH "c:/bin/"
#endif

int main(int argc, char* argv[])
{
	if(argc<2)
	{
		cout << "Syntax : aquireWisdom [path]filename" << endl;
		return EXIT_FAILURE;
	}
	string inputFile(argv[1]);
	try
	{
	    LifReader reader(inputFile);
	    const size_t serie = reader.chooseSerieNumber();
	    LifTracker track(reader.getSerie(serie),0, FFTW_EXHAUSTIVE);
		cout << "tracker ok"<<endl;
		Tracker::saveWisdom(INSTAL_PATH "wisdom.fftw");
		cout<<"done !"<<endl;

	}
    catch(exception& e)
    {
        cout << e.what() << "\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
