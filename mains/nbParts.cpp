#include <string>
#include <iostream>
#include <fstream>
#include <valarray>
#include "../files_series.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    if(argc<5)
    {
        cout << "Syntax : nbParts [path]filename token t_offset t_span" << endl;
        return EXIT_FAILURE;
    }

    const string filename(argv[1]);

    vector<string> tokens(1,argv[2]);
    TokenTree tt(tokens,filename);
    cout << tt << endl;

    size_t t_offset,t_span;
    sscanf(argv[3],"%u",&t_offset);
    sscanf(argv[4],"%u",&t_span);

    vector<size_t> v(1,t_offset);

    string trash;
    valarray<double> nbs(0.0,t_span);
    ofstream output((tt.getNoIndexNameNoExt()+"_nb.txt").c_str(), ios::out | ios::trunc);
    output << "t\tNb"<<endl;
    while(v[0]<t_span)
    {
        ifstream input(tt(v).c_str(), ios::in);
        if(!input)
        {
            cerr<<"no such file as "<<tt(v)<<endl;
            return EXIT_FAILURE;
        }
        input >> trash;
        input >> nbs[v[0]-t_offset];
        output << v[0] <<"\t"<< nbs[v[0]-t_offset] << endl;
        input.close();
        v[0]++;
    }
    output.close();
    const double average = nbs.sum()/t_span;
    const double stdev = sqrt(((nbs-average)*(nbs-average)).sum()/t_span);

    cout << "average number "<<average<<" +- "<<stdev<<endl;
    return EXIT_SUCCESS;
}
