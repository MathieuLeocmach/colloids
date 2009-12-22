#include "../particles.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    if(argc<3)
    {
        cout << "Syntax : sim2xp [path]filename radius" << endl;
        cout << " convert data with periodic boundary conditions (simulation) to hard wall boundary condition (experiment)" << endl;
        return 1;
    }
    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    float radius;
    sscanf(argv[2],"%f",&radius);

    particles Centers(filename);

    valarray<double> trans((double)radius,3);
    Centers+=trans;
    for(size_t i=0;i<3;++i)
    {
        Centers.bb.edges[i].first  -= (double)radius;
        Centers.bb.edges[i].second += (double)radius;
    }
    Centers.exportToFile(inputPath + "_xp.dat");
    return EXIT_SUCCESS;
}
