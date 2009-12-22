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


#include <fstream>
#include <iostream>
#include <valarray>
using namespace std;

int main(int argc, char ** argv)
{
	if(argc<7)
    {
        cout << "Syntax : histogram [path]filename nbcolumns usefullColumn min max Nbins" << endl;
        return EXIT_FAILURE;
    }
    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    size_t nbcolumns,usefullColumn,Nbins;
    double minV,maxV;
    sscanf(argv[2],"%u",&nbcolumns);
    sscanf(argv[3],"%u",&usefullColumn);
    sscanf(argv[4],"%lf",&minV);
    sscanf(argv[5],"%lf",&maxV);
    sscanf(argv[6],"%u",&Nbins);

    //open file
    ifstream file(argv[1], ios::in);
    if(!file)
    {
        cerr<<"no such file as "<<filename<<endl;
        return EXIT_FAILURE;
    }

    //get the name of the usefull column
    string trash, name;
    size_t col=0;
    while(col<usefullColumn)
    {
    	file>>trash;
    	col++;
    }
    file>>name;
    col++;
    while(col<nbcolumns)
    {
    	file>>trash;
    	col++;
    }

    //bins to put the histogram in
    valarray<double> bins(0.0,Nbins+1);
    //binning
    double val;
    size_t bin;
    while(!file.eof())
    {
    	col=0;
    	while(col<usefullColumn)
		{
			file>>trash;
			col++;
		}
		file>>val;
		bin = (size_t)(Nbins*(val-minV)/(maxV-minV));
		if(bin<Nbins)
			bins[bin]++;
		else
			bins[Nbins]++;
		col++;
		while(col<nbcolumns)
		{
			file>>trash;
			col++;
		}
    }
    file.close();
    //normalizing
    bins /= bins.sum();

    //output
    ofstream output((inputPath+"_"+name+".histogram").c_str(), ios::out | ios::trunc);
    output<<"lowerbounds\t"<<name<<endl;
    const double step = (maxV-minV)/Nbins;
    for(size_t i=0;i<bins.size();++i)
		output<< minV+step*i <<"\t"<<bins[i]<<endl;
	output.close();
	return EXIT_SUCCESS;
}
