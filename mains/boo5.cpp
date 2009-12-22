#include "../dynamicParticles.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    if(argc<3)
    {
        cout << "Syntax : boo5 [path]filename.dat radius" << endl;
        return EXIT_FAILURE;
    }
    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    const string ext = filename.substr(filename.find_last_of(".")+1);
    double radius;
    sscanf(argv[2],"%lf",&radius);

    try
    {
        indexedParticles parts(filename,radius);
        cout << parts.size() << " particles"<<endl;

        //five fold bond orientational order
        map<size_t,boo_five> boof,Sboof;
        map<size_t,boo_five>::iterator booit;
        set<size_t> inside = parts.getRealInside(1.4*2.0*radius);
        for(set<size_t>::iterator p=inside.begin();p!=inside.end();++p)
        {
            set<size_t> EuNgb = parts.getEuclidianNeighbours(parts[*p],1.4*2.0*radius);
            if(EuNgb.size()>1)
            {
                booit = boof.insert(boof.end(),make_pair(*p,boo_five()));
                for(set<size_t>::iterator n=EuNgb.begin();n!=EuNgb.end();++n)
                    if( *p != *n)
                        booit->second += boo_five(parts.getDiff(*p,*n));
                booit->second /= (double)(EuNgb.size()-1);
            }
        }

        //coarse grained
        set<size_t> second_inside = parts.getRealInside(2.0*1.4*2.0*radius);
        map<size_t,boo_five>::const_iterator bofit;
        for(set<size_t>::iterator p=second_inside.begin();p!=second_inside.end();++p)
        {
            set<size_t> EuNgb = parts.getEuclidianNeighbours(parts[*p],1.4*2.0*radius);
            booit = Sboof.insert(Sboof.end(),make_pair(*p,boo_five()));
            for(set<size_t>::iterator n=EuNgb.begin();n!=EuNgb.end();++n)
            {
                bofit=boof.find(*p);
                if(bofit!=boof.end());
                    booit->second += (*bofit).second;
            }
        }

        vector<scalarField> q5(2);
        q5[0].first = new string("q5");
        q5[1].first = new string("Sq5");
        q5[0].second = new map<size_t,double>();
        q5[1].second = new map<size_t,double>();
        for(bofit=boof.begin();bofit!=boof.end();++bofit)
            q5[0].second->insert(q5[0].second->end(),make_pair(bofit->first,bofit->second.getQ5()));
        for(bofit=Sboof.begin();bofit!=Sboof.end();++bofit)
            q5[1].second->insert(q5[1].second->end(),make_pair(bofit->first,bofit->second.getQ5()));

        parts.exportToVTK(inputPath+"_five.vtk",q5);
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

