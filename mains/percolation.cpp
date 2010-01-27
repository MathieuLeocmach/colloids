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

#include "../dynamicClusters.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    if(argc<3)
    {
        cout << "Syntax : percolation [path]filename.dat radius range" << endl;
        cout << "\tOR"<<endl;
        cout << "\tpercolation [path]filename.traj range" << endl;
        cout << " range is in diameter unit" << endl;
        return EXIT_FAILURE;
    }

    cout << "Cluster percolation" << endl;
    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    const string ext = filename.substr(filename.find_last_of(".")+1);

    double radius,range;

    if(ext.compare("traj")==0)
    {
        try
        {
            DynamicParticles parts(filename);
            //cout << parts.trajectories.size() << " particles"<<endl;
            radius = parts.radius;
            sscanf(argv[2],"%lf",&range);

            //create the population of trajectories to split into clusters : all trajectories
            set<size_t> alltraj;
            for(size_t tr=0;tr<parts.trajectories.size();++tr)
                alltraj.insert(alltraj.end(),tr);
            cout<<"finding the clusters ... ";
            boost::ptr_vector< vector< set<size_t> > > ngbList;
            ngbList.swap(parts.getNgbList(range));
            DynamicClusters clusters(parts,alltraj,range);
            cout<<"done!"<<endl;

            vector<scalarDynamicField> scalars(1);
            scalars.front() = clusters.getLabels();
            vector<vectorDynamicField> vectors;
            parts.exportToVTK(clusters.getLabels(),vectors,1,ngbList);

            ofstream output((inputPath+"_clusters_largest.txt").c_str(), ios::out | ios::trunc);
            output<<"t\tmaxX\tmaxY\tmaxZ"<<endl;
            valarray<double> largest(0.0,3),dims=largest;

            for(size_t d=0;d<3;++d)
                dims[d] = parts.positions.front().bb.edges[d].second - parts.positions.front().bb.edges[d].first - 2.0*parts.radius;

            for(size_t t=0;t<parts.positions.size();++t)
            {
                largest=clusters.getLargestDelta(t);
                output << t*parts.dt;
                for(size_t d=0;d<3;++d)
                    output<<"\t"<<largest[d] /dims[d];
                output<<endl;
            }
            output.close();

        }
        /*catch(const IdTrajError &e)
        {
            cerr<< e.what()<<endl;
            return EXIT_FAILURE;
        }*/
        catch(const exception &e)
        {
            cerr<< e.what()<<endl;
            return EXIT_FAILURE;
        }
    }
    else
    {
        if(argc<4)
        {
            cout << "Syntax : percolation [path]filename.dat radius range" << endl;
            cout << "\tOR "<<endl;
            cout << "\tpercolation [path]filename.traj range" << endl;
            cout << " range is in diameter unit" << endl;
            return EXIT_FAILURE;
        }
        sscanf(argv[2],"%lf",&radius);
        sscanf(argv[3],"%lf",&range);

        IndexedParticles Centers(filename,radius);
        cout << Centers.size() << " particles ... ";

        //neighbour list
        vector< set<size_t> > ngbList;
        Centers.getNgbList(1.3,ngbList);
        //cluster data container
        deque< set<size_t> > clusters;
        //get the clusters
        segregateAll(clusters, Centers);

        //labelling
        string name = "clusters";
        vector<scalarField> scalars(1);
        scalars.front().first = *name;
        scalars.front().second = new map<size_t,double>();

        Centers.labels.assign(Centers.size(),(unsigned char)0);
        double label = 1;
        for(deque< set<size_t> >::const_iterator k=clusters.begin();k!=clusters.end();++k)
            if((*k).size()>1)
            {
                transform(
                    k->begin(),k->end(),
                    inserter(*scalars.front().second,scalars.front().second->end()),
                    std::bind2nd(make_pair,label)
                    );
                label++;
            }


        Centers.exportToVTK(inputPath + "_clusters",scalars);
    }
    return EXIT_SUCCESS;
}

