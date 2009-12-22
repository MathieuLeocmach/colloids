#include "../files_series.hpp"
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

//#include "../indexedParticles.hpp"
#include "../saveTable.hpp"
#include <valarray>

using namespace std;

int main(int argc, char ** argv)
{
    if(argc<6)
    {
        cout << "Syntax : sumrdf [path]filename token nbcolumns t_offset t_span" << endl;
        return 1;
    }

    const string filename(argv[1]);//,postfix(argv[6]);
    const string ext = filename.substr(filename.find_last_of(".")+1);

    vector<string> tokens(1,argv[2]);
    TokenTree tt(tokens,filename);

    size_t nbcolumns,t_offset,t_span;
    sscanf(argv[3],"%u",&nbcolumns);
    sscanf(argv[4],"%u",&t_offset);
    sscanf(argv[5],"%u",&t_span);

    //reading the dimension of the file at t=t_offset
    vector<size_t> v(1,t_offset);
    ifstream firstFile(tt(v).c_str(), ios::in);
    if(!firstFile)
    {
        cerr<<"no such file as "<<tt(v)<<endl;
        cerr << tt << endl;
        return EXIT_FAILURE;
    }
    string line,header;
    size_t linecount=0;
    getline(firstFile, header);
    while(getline(firstFile, line))
        linecount++;

    firstFile.close();

    valarray< valarray<double> > g(valarray<double>(0.0,nbcolumns),linecount);
    valarray<double>aline(0.0,nbcolumns);

    while(v[0]<t_span)
    {
        ifstream input(tt(v).c_str(), ios::in);
        if(!input)
        {
            cerr<<"no such file as "<<tt(v)<<endl;
            return EXIT_FAILURE;
        }
        getline(input, line);
        for(size_t l=0;l<g.size();++l)
        {
            for(size_t c=0;c<g[l].size();++c)
                input >> aline[c];
            g[l]+=aline;
        }
        input.close();
        v[0]++;
    }
    g/=valarray<double>((double)t_span,nbcolumns);

    string outputname = tt.getNoIndexNameNoExt()+"_total"+ext;
    cout << "export to " << outputname << endl;
     ofstream output(outputname.c_str(), ios::out | ios::trunc);
     if(output)
     {
       //header line
       output << header << endl;

       for(size_t l=0;l<g.size();++l)
       {
            output << g[l][0];
            for(size_t c=1;c<g[l].size();++c)
                output<< "\t"<<g[l][c];
            output << endl;
       }
       output.close();
     }
     else cout << " cannot open the output file";

    return EXIT_SUCCESS;
}

