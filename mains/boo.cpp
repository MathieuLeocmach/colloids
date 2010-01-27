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

//Define the preprocessor variable "periodic" if you want periodic boundary conditions
#include "../periodic.hpp"
//#include "../pv.hpp"
#include "../dynamicParticles.hpp"

using namespace std;
using namespace Colloids;

int errorMessage()
{
#ifdef use_periodic
    cout << "Syntax : periodic_boo [path]filename.grv radius Nb Dx Dy Dz" << endl;
#else
    cout << "Syntax : boo [path]filename.dat radius" << endl;
    cout << "\tOR\t  boo [path]filename.traj (AveragingInterval)" << endl;
#endif
    return EXIT_FAILURE;
}

int main(int argc, char ** argv)
{
    if(argc<2) return errorMessage();

    const string filename(argv[1]);
    const string inputPath = filename.substr(0,filename.find_last_of("."));
    const string ext = filename.substr(filename.find_last_of(".")+1);

    //const double shell = 1.3;

    try
    {

        if(ext.compare("traj")==0)
        {

            DynamicParticles parts(filename);
            const size_t stop = parts.getNbTimeSteps()-1;
            cout << parts.trajectories.size() << " particles in "<<stop+1<<" time steps"<<endl;

            parts.removeDrift();
            cout<<"drift removed"<<endl;

			size_t tau;
            if(argc<3)
            {
            	//Dynamics calculation and export
				vector< vector<double> > MSD(1);
				vector< vector<double> > ISF;
				parts.makeDynamics(MSD.front(),ISF);
				saveTable(MSD.begin(),MSD.end(),inputPath + ".msd","t\tMSD",parts.dt);
				saveTable(ISF.begin(),ISF.end(),inputPath + ".isf","t\tx\ty\tz\tav",parts.dt);

				//Find relaxation time
				if(ISF.back().back() > exp(-1.0))
					tau = parts.getNbTimeSteps();
				else
					tau = parts.getNbTimeSteps() - (upper_bound(ISF.back().rbegin(),ISF.back().rend(),exp(-1.0))-ISF.back().rbegin());
				cout<<"relaxation time is "<<tau<<" steps, ie "<<tau*parts.dt<<"s"<<endl;
            }
            else
				sscanf(argv[2],"%u",&tau);
            //const size_t halfInterval = tau/2;

            vector< map<size_t,valarray<double> > > qw,Sqw,Tqw,STqw;
            cout<<"get raw Boo from files"<<endl;
            set<size_t> inside = parts.getBooFromFile("",qw);
            cout<<"get coarse grained Boo from files"<<endl;
            set<size_t> secondInside = parts.getBooFromFile("_space",Sqw);

            cout<<"Neighbour lists ";
            const double range = 1.3;
			DynamicNeigbourList ngbList = parts.getNgbList(range);
			cout<<"and lost neighbours"<<endl;
			vector<scalarDynamicField> scalars(9);
            scalars[8] = parts.getNbLostNgbs(inside,ngbList,tau);

            cout<<"Time averaging over Tau="<<parts.dt*tau<<" s ... ";
            parts.makeSlidingTimeAverage(inside,tau,qw,Tqw);
            parts.makeSlidingTimeAverage(secondInside,tau,Sqw,STqw);
            vectorDynamicField vel = parts.averageVelocities(inside,1,tau);
            cout<<"done!"<<endl;

            ofstream cloudT((inputPath + "_time.cloud").c_str(), ios::out | ios::trunc);
            cloudT << "#tr\tt\tQ4\tQ6\tW4\tW6"<<endl;
            for(vector< map<size_t,valarray<double> > >::const_iterator t=Tqw.begin();t!=Tqw.end();++t)
                for(map<size_t,valarray<double> >::const_iterator tr=t->begin();tr!=t->end();++tr)
                {
                    cloudT<<t-Tqw.begin()<<"\t"<<(*tr).first;
                    for(size_t i=0;i<tr->second.size();++i)
						cloudT<<"\t"<<tr->second[i];
					cloudT<<endl;
                }
            cloudT.close();

            ofstream cloudST((inputPath + "_space_time.cloud").c_str(), ios::out | ios::trunc);
            cloudST << "#tr\tt\tQ4\tQ6\tW4\tW6"<<endl;
            size_t ti;
            for(vector< map<size_t,valarray<double> > >::const_iterator t=STqw.begin();t!=STqw.end();++t)
            {
            	ti = t-STqw.begin();
                for(map<size_t,valarray<double> >::const_iterator tr=t->begin();tr!=t->end();++tr)
                {
                    cloudST<<tr->first<<"\t"<<ti;
                    for(size_t i=0;i<tr->second.size();++i)
						cloudST<<"\t"<<tr->second[i];
					cloudST<<"\t"<< ((ti<vel.second->size())?
						sqrt(pow(vel.second->at(ti)[tr->first],2.0).sum()):
						sqrt(pow(vel.second->back()[tr->first],2.0).sum())
					);
					cloudST<<"\t"<< ((ti<scalars[8].second->size())?
						(size_t)scalars[8].second->at(ti)[tr->first]:
						(size_t)scalars[8].second->back()[tr->first]
					);
					cloudST<<endl;
                }
            }
            cloudST.close();

            //make 4 scalarDynamicfields out of boo data
            //vector<scalarDynamicField> scalars(8);
            //string[8] fields_names = {"Q4","Q6","W4","W6"};
            for(size_t i=0;i<8;++i)
            {
				scalars[i].first = string(i<4?"":"cg")+string((i%4<2)?"Q":"W")+string((i%2)?"6":"4");
				/*fields_names[i];
				scalars[i+4].first = "cg"+fields_names[i];*/
				scalars[i].second = new vector< map<size_t,double> >(Tqw.size());
            }
            for(size_t t=0;t<Tqw.size();++t)
            {
                for(map<size_t,valarray<double> >::const_iterator tr=Tqw[t].begin();tr!=Tqw[t].end();++tr)
                	for(size_t i=0;i<tr->second.size();++i)
						scalars[i].second->at(t).insert(scalars[i].second->at(t).end(),make_pair(tr->first,tr->second[i]));
				for(map<size_t,valarray<double> >::const_iterator tr=STqw[t].begin();tr!=STqw[t].end();++tr)
                	for(size_t i=0;i<tr->second.size();++i)
						scalars[i+4].second->at(t).insert(scalars[i+4].second->at(t).end(),make_pair(tr->first,tr->second[i]));
            }

			//make 1 vector field out of velocity data
			vector<vectorDynamicField> vectors(1,vel);

		    parts.exportToVTK(scalars,vectors);

		    cout<<"VTK exported"<<endl;


        }
        else
        {
            if(argc<3) return errorMessage();
            double radius;
            sscanf(argv[2],"%lf",&radius);

            //construct the particle container out of the datafile
#ifdef use_periodic
            if(argc<7) return errorMessage();
            size_t Nb;
            sscanf(argv[3],"%u",&Nb);
            BoundingBox b;
            for(size_t d=0;d<3;++d)
            {
                b.edges[d].first=0.0;
                sscanf(argv[4+d],"%lf",&b.edges[d].second);
            }
            PeriodicParticles parts(Nb,b,radius,filename);
#else
            IndexedParticles parts(filename,radius);
#endif
            //cout << parts.size() << " particles"<<endl;

            const double shell =1.3;
            vector< set<size_t> >ngbList;
			parts.getNgbList(shell,ngbList);
			ofstream ngb_file((inputPath + ".ngb").c_str(), ios::out | ios::trunc);
			for(vector< set<size_t> >::const_iterator it =ngbList.begin();it!=ngbList.end();++it)
				ngb_file<<it->size()-1<<endl;
			ngb_file.close();

            //cout<<"get bond orientational order data for each usefull particles ... ";
            set<size_t> inside = parts.getRealInside(shell*2.0*radius);
            map<size_t,BooData>allBoo;
            for(set<size_t>::const_iterator p=inside.begin();p!=inside.end();++p)
                allBoo.insert(allBoo.end(),make_pair(*p,parts.getBOO(*p,ngbList[*p])));

            //Disabled. Take disk space for no use
            //parts.exportQlm(allBoo,inputPath+".qlm");

            //cout<<"Coarse graining ... ";
            //because we take the first shell, margin*2
            set<size_t> secondInside = parts.getRealInside(2.0*shell*2.0*radius);
            map<size_t,BooData> SallBoo;
            for(set<size_t>::const_iterator p=secondInside.begin();p!=secondInside.end();++p)
                SallBoo.insert(SallBoo.end(),make_pair(*p,parts.getAvBOO(allBoo,*p,ngbList[*p])));

            parts.exportQ6m(SallBoo,inputPath+".q6m");

            /*vector< set<size_t> > remarkableSets(3);
            vector< string > remarkableSetsNames(3);
            remarkableSetsNames[0] = "all";
            remarkableSetsNames[1] = "q6>0.5";
            remarkableSetsNames[2] = "Sq6>0.3";
            const double q6threshold=0.5, Sq6threshold=0.3;*/
            valarray<double> vboo(0.0,4);

            //cout<<"Export data ... ";
            ofstream cloud((inputPath + ".cloud").c_str(), ios::out | ios::trunc);
            if(!cloud)
                throw invalid_argument("No such file as "+inputPath + ".cloud");
            cloud << "#p\tQ4\tQ6\tW4\tW6"<<endl;
            for(map<size_t,BooData>::const_iterator p=allBoo.begin();p!=allBoo.end();++p)
            {
                cloud<<p->first;
                (*p).second.getInvarients(4,vboo[0],vboo[2]);
                (*p).second.getInvarients(6,vboo[1],vboo[3]);
                for(size_t i=0;i<vboo.size();++i)
					cloud<<"\t"<<vboo[i];
                cloud<<endl;
                //remarkableSets[0].insert(remarkableSets[0].end(),(*p).first);
                //if(val_q6 >q6threshold)
                //    remarkableSets[1].insert(remarkableSets[1].end(),(*p).first);
            }
            cloud.close();
            ofstream Scloud((inputPath + "_space.cloud").c_str(), ios::out | ios::trunc);
            Scloud << "#p\tQ4\tQ6\tW4\tW6"<<endl;
            for(map<size_t,BooData>::const_iterator p=SallBoo.begin();p!=SallBoo.end();++p)
            {
                Scloud<<p->first;
                (*p).second.getInvarients(4,vboo[0],vboo[2]);
                (*p).second.getInvarients(6,vboo[1],vboo[3]);
                for(size_t i=0;i<vboo.size();++i)
					Scloud<<"\t"<<vboo[i];
                Scloud<<endl;
                //if(val_q6 >Sq6threshold)
                    //remarkableSets[2].insert(remarkableSets[2].end(),(*p).first);
            }
            Scloud.close();
            //cout<<"done"<<endl;
            //parts.rdf_angD(remarkableSets,remarkableSetsNames,inputPath);
        }
    }
    catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
