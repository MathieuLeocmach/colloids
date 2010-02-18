#include "../dynamicParticles.hpp"

using namespace std;
using namespace Colloids;

int main(int argc, char ** argv)
{
    if(argc<2)
    {
        cout << "Syntax : MSD [path]filename.traj start1 stop1 av1 [start2 stop2 av2 [...]]" << endl;
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
        vector<double> MSD;
        boost::format name (inputPath+"_%1%from_%2%to_%3%av.msd");
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
        	MSD = parts.getMSD(start,stop,av);

        	//export to file
        	ofstream output((name % start %stop %av).str().c_str(), std::ios::out | std::ios::trunc);
        	output <<"#t\tMSD"<<endl;
        	for(size_t t=0;t<MSD.size();++t)
				output <<t*parts.dt<<"\t"<<MSD[t]<<"\n";
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


