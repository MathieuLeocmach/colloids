#include "../dynamicParticles.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    if(argc<2)
    {
        cout << "Syntax : ISF [path]filename mode=0" << endl;
        cout <<"\t modes are:"<<endl;
        cout <<"\t0: complete ISF"<<endl;
        cout <<"\t1: self-ISF"<<endl;
        cout <<"\t2: both"<<endl;
        return EXIT_FAILURE;
    }

    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));

    size_t mode =0;
    if(argc>2)
        sscanf(argv[2],"%u",&mode);

    dynamicParticles parts(filename);

    cout << "get Intermediate scattering function between t=0 and t=" << parts.getNbTimeSteps()-1 << " ...";
    if(mode==0 || mode==2)
    {
        valarray<double> q(0.0,3);
        q[0]=M_PI/parts.radius;
        vector<double> ISF = parts.getISF(q,0,parts.getNbTimeSteps()-1);
        cout << " done !"<<endl;

        cout << "export to " << inputPath<<"_ISF.txt" << endl;

        ofstream output((inputPath+"_ISF.txt").c_str(), ios::out | ios::trunc);
        if(output)
        {
            output<<"t\tISF"<<endl;
            for(size_t t=0;t<ISF.size();++t)
                output << t*parts.dt << "\t" << ISF[t]<<endl;
            output.close();
        }
    }
    if(mode==1 || mode==2)
    {
        cout << "self part "<< endl;
        vector< vector<double> > ISF(3,vector<double>(parts.getNbTimeSteps()-1,0));
        try
        {
            for(size_t d=0;d<3;++d)
            {
                cout<<d<<" ... "<<endl;
                valarray<double> q(0.0,3);
                q[d]=M_PI/parts.radius;
                ISF[d] = parts.getSelfISF(q,0,parts.getNbTimeSteps()-1);
            }
            cout << " done !"<<endl;
        }
        catch(const exception &e)
        {
            cerr << e.what() << endl;
            return EXIT_FAILURE;
        }
        saveTable(ISF.begin(),ISF.end(),inputPath+"_SelfISF.txt","t\tx\ty\tz",parts.dt);
    }

    return EXIT_SUCCESS;
}



