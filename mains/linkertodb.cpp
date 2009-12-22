#include "../dynamicParticles.hpp"
#include "../sqlite3x/sqlite3x.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    if(argc<9)
    {
        cout << "Syntax : linkertodb dbname measurement [path]filename token radius time_step t_offset t_span" << endl;
        return EXIT_FAILURE;
    }

    const string dbname(argv[1]),filename(argv[3]), token(argv[4]);
    double radius,time_step;
    sscanf(argv[5],"%lf",&radius);
    sscanf(argv[6],"%lf",&time_step);
    size_t measurement,t_offset,t_span;
    sscanf(argv[2],"%u",&measurement);
    sscanf(argv[7],"%u",&t_offset);
    sscanf(argv[8],"%u",&t_span);

    try
    {
        cout<<"linking :"<<endl;
        dynamicParticles parts(radius,time_step,filename,token,t_offset,t_span,0.3);
        cout << "total of " << parts.trajectories.size() << " trajectories" << endl;
        for(size_t t=0;t<parts.getNbTimeSteps();++t)
        {
            try{parts.positions[t].tree.getOverallBox();}
            catch(exception &e)
            {
                cerr<<"At time "<<t<<" before drift removal: "<<e.what()<<endl;
                exit(1);
            }
        }
        parts.removeDrift();
        cout<<"drift removed"<<endl;
        //exit(1);

        const string noExt = filename.substr(0,filename.find_last_of("."));
        const string nofolder = filename.substr(filename.find_last_of("/\\")+1);
        const string folder = filename.substr(0,filename.find_last_of("/\\")+1);
        parts.save(noExt+".traj",nofolder,token,t_offset,t_span);

        animation a = parts.getAnimation(t_offset,t_span);
        cout << "animation is " << a.size() << " frames" << endl;
        a.exportToPV(noExt+"_linked.pv");

        cout<<"saving to database"<<endl;
        parts.exportToDb(dbname,measurement);
        cout <<"done !"<<endl;
    }
    /*catch(const sqlite3x::database_error &e)
    {
        cout << e.what()<<endl;
        return EXIT_FAILURE;
    }*/
    catch(const exception &e)
    {
        cout << e.what()<<endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

