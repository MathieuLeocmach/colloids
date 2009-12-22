#include "../dynamicParticles.hpp"

int main(int argc, char ** argv)
{
    if(argc<2)
    {
        cout << "Syntax : drift [path]filename" << endl;
        return EXIT_FAILURE;
    }

    const string filename(argv[1]);
    const string noExt = filename.substr(0,filename.find_last_of("."));


    try
    {
        dynamicParticles parts(filename);
        cout << "total of " << parts.trajectories.size() << " trajectories" << endl;

        const size_t last_frame = parts.positions.size()-1;

        cout << "export to " << noExt+"_drift.txt" << endl;

        ofstream output((noExt+"_drift.txt").c_str(), ios::out | ios::trunc);
        valarray<double> drift(0.0,3);
        output << "t\tx\ty\tz" << endl;

        for(size_t t=0;t<last_frame;++t)
        {
            drift = parts.getDrift(parts.getSpanning(t,t+1),t,t+1);
            output << t;
            for(size_t i=0;i<3;++i)
                output << "\t" << drift[i];
            output << endl;
        }
        output.close();
    }
    catch(trajerror tre)
    {
        cout << tre;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

