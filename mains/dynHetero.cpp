#include "../dynamicParticles.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    if(argc<2)
    {
        cout << "Syntax : dynHetero [path]filename" << endl;
        return EXIT_FAILURE;
    }

    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));

    try
    {
        dynamicParticles parts(filename);

        const size_t stop = parts.positions.size()-1;
        set<size_t> spanning = parts.getSpanning(0,stop);
        valarray<double> SD(0.0,spanning.size());
        valarray<double> diff(0.0,3);
        size_t i=0;
        for(set<size_t>::iterator tr=spanning.begin();tr!=spanning.end();++tr)
        {
            diff=parts.getDiff(*tr,0,*tr,stop);
            SD[i++]=(diff*diff).sum();
        }

        cout << "export square displacements to " << inputPath<<"_SD.txt" << endl;

        ofstream output((inputPath+"_SD.txt").c_str(), ios::out | ios::trunc);
        if(output)
        {
            output<<"t\tSD"<<endl;
            for(size_t i=0;i<SD.size();++i)
                output << i << "\t" << SD[i]<<endl;
            output.close();
        }
        const double meanD = SD.apply(sqrt).sum()/SD.size();
        const double varD = sqrt(SD.sum()/SD.size()-meanD*meanD);
        const double uptail = pow(meanD+varD,2.0);
        const double downtail = pow(meanD-varD,2.0);
        i=0;
        for(set<size_t>::iterator tr=spanning.begin();tr!=spanning.end();++tr)
        {
            if(SD[i]>uptail)
                parts.setTrajLabel(*tr,14);
            else
            {
                if (SD[i]<downtail)
                    parts.setTrajLabel(*tr,1);
                else
                    parts.setTrajLabel(*tr,7);
            }
            i++;
        }


        animation a= parts.getAnimation(0,1);
        a.exportToPV(inputPath+"_dynH.pv");

        return EXIT_SUCCESS;
    }
    catch(const trajerror &tre)
    {
        cerr<< tre.what()<<endl;
        return EXIT_FAILURE;
    }
}


