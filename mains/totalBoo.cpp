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

#include "../dynamicParticles.hpp"
#include <boost/timer.hpp>
#include <boost/progress.hpp>

using namespace std;

int main(int argc, char ** argv)
{
	if(argc<2)
	{
		cerr <<"SYNTAX : totalBoo [path]file.traj"<<endl;
		return EXIT_FAILURE;
	}
	try
	{
		const string filename(argv[1]);
		const string inputPath = filename.substr(0,filename.find_last_of("."));
		const string ext = filename.substr(filename.find_last_of(".")+1);

		DynamicParticles parts(filename);
		const size_t stop = parts.getNbTimeSteps()-1;
		cout << parts.trajectories.size() << " particles in "<<stop+1<<" time steps"<<endl;
		parts.removeDrift();
		cout<<"drift removed"<<endl;

		vector<size_t> v(1,parts.trajectories.t_offset);
		const string positionFile = (parts.trajectories.tt % parts.trajectories.t_offset).str();
		const string base = positionFile.substr(0,positionFile.find_last_of("."));
		TokenTree base_tt(parts.trajectories.tt.tokens,base);

		//Dynamics calculation and export
		vector< vector<double> > MSD(1);
		vector< vector<double> > ISF;
		parts.makeDynamics(MSD.front(),ISF);
		saveTable(MSD.begin(),MSD.end(),inputPath + ".msd","t\tMSD",parts.dt);
		saveTable(ISF.begin(),ISF.end(),inputPath + ".isf","t\tx\ty\tz\tav",parts.dt);

		//Find relaxation time
		size_t tau;
		if(ISF.back().back() > exp(-1.0))
			tau = parts.getNbTimeSteps();
		else
			tau = parts.getNbTimeSteps() - (upper_bound(ISF.back().rbegin(),ISF.back().rend(),exp(-1.0))-ISF.back().rbegin());
		cout<<"relaxation time is "<<tau<<" steps, ie "<<tau*parts.dt<<"s"<<endl;

		cout<<"Treating each time step"<<endl;
		boost::progress_display show_progress(parts.getNbTimeSteps());

		//declaration of the data containers
		const double range = 1.3;
		DynamicNeigbourList ngbList = parts.getNgbList(range);
		vector< map<size_t,valarray<double> > > qw(parts.getNbTimeSteps()),Sqw(parts.getNbTimeSteps()),Tqw,STqw;
		//identify trajectories spanning the whole time steps
		set<size_t> insideTr = parts.getSpanning(0,stop), secondInsideTr=insideTr;

		clock_t interT =0;


		for(size_t t=0;t<parts.getNbTimeSteps();++t)
		{
			set<size_t> inside = parts.positions[t].getRealInside(range*2.0*parts.radius),
				secondInside = parts.positions[t].getRealInside(2.0*range*2.0*parts.radius);

			//boo containers
			map<size_t,BooData>boo;
			map<size_t,BooData>::iterator b;
			valarray<double> vboo(0.0,4);
			ofstream cloud(((base_tt% (t+parts.trajectories.t_offset)).str() + ".cloud").c_str(), ios::out | ios::trunc);
			cloud << "#p\tQ4\tQ6\tW4\tW6"<<endl;
            ofstream Scloud(((base_tt% (t+parts.trajectories.t_offset)).str() + "_space.cloud").c_str(), ios::out | ios::trunc);
            Scloud << "#p\tQ4\tQ6\tW4\tW6"<<endl;

			//raw boo
			for(set<size_t>::const_iterator p=inside.begin();p!=inside.end();++p)
			{
				cloud<<*p<<"\t";
				b = boo.insert(boo.end(), make_pair(*p,parts.positions[t].getBOO(*p,ngbList[t][*p])) );
				b->second.getInvarients(4,vboo[0],vboo[2]);
				b->second.getInvarients(6,vboo[1],vboo[3]);
				qw[t].insert(qw[t].end(),make_pair(*p,vboo));
				for(size_t i=0;i<vboo.size();++i)
					cloud<<"\t"<<vboo[i];
                cloud<<"\n";
			}
			cloud.close();
			//disabled to save disk space
			//parts.positions[t].exportQlm(boo,(base_tt% (t+parts.trajectories.t_offset)).str()+".qlm");

            //coarse grained boo
            BooData Sb;
			for(set<size_t>::const_iterator p=secondInside.begin();p!=secondInside.end();++p)
			{
				Scloud<<*p<<"\t";
				//b = Sboo.insert(Sboo.end(), make_pair(*p,parts.positions[t].getAvBOO(boo,*p,ngbList[t])) );
				Sb = parts.positions[t].getAvBOO(boo,*p,ngbList[t][*p]);
				Sb.getInvarients(4,vboo[0],vboo[2]);
				Sb.getInvarients(6,vboo[1],vboo[3]);
				Sqw[t].insert(Sqw[t].end(),make_pair(*p,vboo));
				for(size_t i=0;i<vboo.size();++i)
					Scloud<<"\t"<<vboo[i];
                Scloud<<"\n";
			}
			cloud.close();
			parts.positions[t].exportQ6m(boo,(base_tt% (t+parts.trajectories.t_offset)).str()+".q6m");

			//What are the trajectories that have boo in all time steps ?
			//Which one is the quickest between that and the STindex ? => to be tested
			clock_t startInter = clock();
			/*if(t==0)
			{
				transform(inside.begin(),inside.end(),inserter(insideTr,insideTr.end()),TrajIndex::Inverser(t,parts.trajectories));
				transform(secondInside.begin(),secondInside.end(),inserter(secondInsideTr,secondInsideTr.end()),TrajIndex::Inverser(t,parts.trajectories));
			}
			else
			{*/
				set<size_t> transf,inter;
				transform(
					inside.begin(),
					inside.end(),
					inserter(transf,transf.end()),
					TrajIndex::Inverser(t,parts.trajectories)
					);
				set_intersection(insideTr.begin(),insideTr.end(),transf.begin(),transf.end(),inserter(inter,inter.end()));
				insideTr.swap(inter);

				transf.clear();
				inter.clear();
				transform(secondInside.begin(),secondInside.end(),inserter(transf,transf.end()),TrajIndex::Inverser(t,parts.trajectories));
				set_intersection(secondInsideTr.begin(),secondInsideTr.end(),transf.begin(),transf.end(),inserter(inter,inter.end()));
				secondInsideTr.swap(inter);
			//}
			interT += clock() - startInter;
			++show_progress;
		}
		/*cout<<interT/(double)CLOCKS_PER_SEC<<"s to get insides by insersection"<<endl;
		boost::timer STindexT;
		set<size_t> a = parts.getSpanningInside(0,stop,1.3*2.0*parts.radius);
		set<size_t> b = parts.getSpanningInside(0,stop,2.0*1.3*2.0*parts.radius);
		cout<<STindexT.elapsed()<<"s to get insides by STindex"<<endl;
		if(equal(a.begin(),a.end(),insideTr.begin()))
			cout <<"insides are equal"<<endl;
		if(equal(b.begin(),b.end(),secondInsideTr.begin()))
			cout <<"second insides are equal"<<endl;
		else
		{
			cout<<insideTr.size()<<"!="<<a.size()<<"?"<<endl;
			cout<<secondInsideTr.size()<<"!="<<b.size()<<"?"<<endl;
		}*/
		cout<<"Preparing data for export"<<endl;
		vector<scalarDynamicField> scalars(9);
		scalars[8] = parts.getNbLostNgbs(insideTr,ngbList,tau);
		parts.makeSlidingTimeAverage(insideTr,tau,qw,Tqw);
		parts.makeSlidingTimeAverage(secondInsideTr,tau,Sqw,STqw);
		vectorDynamicField vel = parts.averageVelocities(secondInsideTr,1,tau);

		ofstream cloudT((inputPath + "_time.cloud").c_str(), ios::out | ios::trunc);
		cloudT << "#tr\tt\tQ4\tQ6\tW4\tW6"<<endl;
		for(vector< map<size_t,valarray<double> > >::const_iterator t=Tqw.begin();t!=Tqw.end();++t)
			for(map<size_t,valarray<double> >::const_iterator tr=t->begin();tr!=t->end();++tr)
			{
				cloudT<<t-Tqw.begin()<<"\t"<<(*tr).first;
				for(size_t i=0;i<tr->second.size();++i)
					cloudT<<"\t"<<tr->second[i];
				cloudT<<"\n";
			}
		cloudT.close();

		ofstream cloudST((inputPath + "_space_time.cloud").c_str(), ios::out | ios::trunc);
		cloudST << "#tr\tt\tQ4\tQ6\tW4\tW6\tV\tLostNgb"<<endl;
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
				cloudST<<"\n";
			}
		}
		cloudST.close();

		//make 4 scalarDynamicfields out of boo data
		//vector<scalarDynamicField> scalars(9);
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

		parts.exportToVTK(scalars,vectors,1,ngbList);

		cout<<"VTK exported"<<endl;

	}
	catch(const exception &e)
    {
        cerr<< e.what()<<endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
