//Define the preprocessor variable "use_periodic" if you want periodic boundary conditions
#include "../periodic.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    if(argc<6)
    {
        cout << "Syntax : [periodic_]ngb_count [path]filename radius start stop nbsteps" << endl;
        cout << " start and stop are in diameter unit" << endl;
        return EXIT_FAILURE;
    }

    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    double radius,start,stop,step;
    sscanf(argv[2],"%lf",&radius);
    sscanf(argv[3],"%lf",&start);
    sscanf(argv[4],"%lf",&stop);
    size_t nbsteps;
    sscanf(argv[5],"%u",&nbsteps);
    step = (stop-start)/(nbsteps+1);

    //construct the particle container out of the datafile
#ifdef use_periodic
    periodicParticles Centers(filename,radius);
    cout << "With periodic boundary conditions"<<endl;
#else
    indexedParticles Centers(filename,radius);
#endif
    cout << Centers.size() << " particles ... ";

    cout << "export to " << inputPath + "_nbNgb.txt" << endl;
    ofstream output((inputPath + "_nbNgb.txt").c_str(), ios::out | ios::trunc);
    if(output)
    {
        output << "cutoff\tN" << endl;

        for(size_t i=0;i<nbsteps+1;++i)
            output<<start+step*i<<"\t"<<Centers.getMeanNbNeighbours((start+step*i)*2.0*radius)<<endl;

        output.close();
    }
    else cout << " cannot open the file";
    return EXIT_SUCCESS;
}

