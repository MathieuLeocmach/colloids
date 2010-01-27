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

#include <ctime>
#include <numeric>
#include <boost/progress.hpp>
#include <boost/bind.hpp>

#include "dynamicParticles.hpp"
#include "files_series.hpp"
#include "saveTable.hpp"

using namespace std;
using namespace Colloids;
//using namespace tvmet;
/**
    \brief Makes the smallest bounding box enclosing all positions of the trajectory
*/
BoundingBox DynamicParticles::boundsTrajectory(const size_t &tr) const
{
	size_t t = trajectories[tr].start_time;
	const size_t last = trajectories[tr].last_time();
    BoundingBox bb = Particles::bounds((*this)(tr,t));

    while(t<last)
        bb.stretch(Particles::bounds((*this)(tr,++t)));

    return bb;
}


/** @brief get the maximum bounding box  */
BoundingBox DynamicParticles::getMaxBox() const
{
	boost::ptr_vector<Particles>::const_iterator it=positions.begin();
	BoundingBox bb = (*it).bb;
	it++;
	while(it!=positions.end())
	{
		bb.stretch((*it).bb);
		it++;
	}
	return bb;
}

/** @brief get the size of the most populated frame  */
size_t DynamicParticles::getMaxSimultaneousParticles() const
{
	size_t N=0;
	for(boost::ptr_vector<Particles>::const_iterator it=positions.begin();it!=positions.end();++it)
		N = max(N,it->size());
	return N;
}

/** \brief Force the spatio-temporal indexation of the trajectories */
/*void DynamicParticles::makeSTindex(const bool reindexAllFrames)
{
    if(reindexAllFrames)
    {
    	boost::progress_display show_progress(positions.size());
        for(size_t t=0;t<positions.size();++t)
        {
            positions[t].makeIndex();
            ++show_progress;
        }
    }

    trajectories.makeInverse();

    //reserve memory
    STindex.clear();
    STindex.reserve(getNbTimeSteps());// = boost::ptr_vector< boost::ptr_vector<RTree> >(positions.size(),boost::ptr_vector<RTree>(positions.size()));
    for(size_t i=0;i<getNbTimeSteps();++i)
    {
    	STindex.push_back(new boost::ptr_vector<RTree>(getNbTimeSteps()));
        for(size_t j=0;j<getNbTimeSteps();++j)
            STindex.back().push_back(new RTree());
    }

    //cout << "STindex initialized ... ";

    for(size_t tr=0;tr<trajectories.size();++tr)
        STindex[trajectories[tr].start_time][trajectories[tr].last_time()].Insert(tr,boundsTrajectory(tr));

}*/

/** @brief constructor of an empty object (no trajectory, no position)
  *
  * \param rad radius of the particles
  * \param time_step The duration of one time step in your time unit
    \param t_size number of time steps
  */
 DynamicParticles::DynamicParticles(const double &rad,const double &time_step,const size_t &t_size)
{
    radius = rad;
    dt = time_step;
    positions.resize(t_size,new Particles(0,rad));
}

/** \brief constructor from the first time step */
DynamicParticles::DynamicParticles(Particles *parts,const double &time_step)
{
    positions.resize(1,parts);
    for(size_t p=0;p<parts->size();++p)
        trajectories.push_back(Traj(0,p));
    dt=time_step;
    radius = parts->radius;
    return;
}

/** \brief constructor from file */
DynamicParticles::DynamicParticles(const string &filename) : trajectories(filename)
{
    radius = trajectories.radius;
    dt = trajectories.dt;
    //cout << dt <<"\t"<<radius << endl;

    //reading the positions
    cout << "Reading position data ... "<<endl;
    //cout << trajectories.tt << endl;

    vector<size_t> v(1,trajectories.t_offset);
    boost::progress_display show_progress(trajectories.t_size-trajectories.t_offset);
    while(v[0]<trajectories.t_size)
        try
        {
            positions.push_back(new Particles(trajectories.tt(v),radius));
            v[0]++;
            ++show_progress;
        }
        catch(const exception &e)
        {
            cerr<<trajectories.tt<<endl;
            throw;
        }
	return;
}

/**
    \brief constructor from position file serie. Links the positions at each time step into trajectories
    \param rad The radius
    \param time_step The duration of one time step in your time unit
    \param filename The name of the trajectory file
    \param base_name A path to a file of the positions files serie, like "myDynamicParticles/series05_t000.dat"
    \param token A substring of base_name just before the time step number, like "_t"
    \param t_offset First time step
    \param t_size number of time steps
*/
DynamicParticles::DynamicParticles(
    const double &rad, const double &time_step,
    const string &base_name, const string &token,
    const size_t &t_offset, const size_t &t_size
    )
{
    radius = rad;
    dt = time_step;

    vector<string> tokens(1,token);
    TokenTree tt(tokens,base_name);

    //displacement vector, dummy if no .displ file
    vector<Coord> displ(t_size, Coord(0));
    double trash;
    ifstream displFile((base_name.substr(0,base_name.rfind(token))+".displ").c_str(), ios::in);
    if(!!displFile)
    {
		//get the relative displacements form file
		cout<<"Using "<< base_name.substr(0,base_name.rfind(token)) << ".displ as an hint for global displacements between frames"<<endl;
		for(size_t t=0;t<t_offset;++t)
			displFile >> trash >>trash;
		for(size_t t=0;t+t_offset<t_size;++t)
			displFile >> displ[t][0] >> displ[t][1]; //only 2D displacements
		displFile.close();

		//compute objective displacement vector
		partial_sum(displ.begin(),displ.end(),displ.begin());
		//largest absolute negative displacement
		Coord maxd(0.0,3);
		for(vector<Coord>::const_iterator t=displ.begin();t!=displ.end();++t)
		{
			maxd[0] = max(maxd[0],(*t)[0]);
			maxd[1] = max(maxd[1],(*t)[1]);
		}
		for(vector<Coord>::iterator t=displ.begin();t!=displ.end();++t)
            *t -= maxd;
    }
    else
		cout<< "Impossible to find "<<base_name.substr(0,base_name.rfind(token)) << ".displ\n"
		"Supposes there is no global displacement between frames"<<endl;

    //push frames by frame
    try
    {
        vector<size_t> v(1,t_offset);

        size_t nbTraj = 0;
        double Error=0;
        double maxError=0;
        double sumError =0;
        boost::progress_display show_progress(t_size);

        while(v[0]<t_offset+t_size)
        {
			Particles p(tt(v), radius);
			if(abs(displ[v[0]-t_offset]).max()>0.0)
				p += displ[v[0]-t_offset];

            push_back(p);

			try{positions.back().getOverallBox();}
			catch(const exception &e)
			{
				cerr<<"At time "<<v[0]<<" after insersion: ";
				throw;
			}
            v[0]++;

            if(v[0]>t_offset+1)
            {
                Error = (trajectories.size() - nbTraj)/(double)trajectories.size();
                sumError+=Error;
                if(maxError<Error) maxError=Error;
            }
            nbTraj = trajectories.size();
            ++show_progress;
        }
        cout<<endl;
        cout<<"Trajectory creation rate : mean="<<100.0*sumError/(positions.size()-1)<<"%\tmax="<<100.0*maxError<<"%"<<endl;
        return;
    }
    catch(const exception &e)
    {
        cerr<<tt<<endl;
        throw;
    }
}

/**
    \brief export the trajectroy data and a relative link to the position files into an ASCII file
    \param filename The name of the trajectory file
    \param base_name A relative link to a file of the positions files serie, like "myDynamicParticles/series05_t000.dat"
    \param token A substring of base_name just before the time step number, like "_t"
    \param t_offset First time step
    \param t_size number of time steps
*/
void DynamicParticles::save(const string &filename,const string &base_name,const string &token,const size_t &t_offset, const size_t &t_size) const
{
    cout << "export to " << filename << endl;

    ofstream output(filename.c_str(), ios::out | ios::trunc);
    if(output)
    {
        //header
        output << radius << "\t" << dt << endl;

        //link to the positions files serie
        output << base_name << endl;
        output << token << endl;
        output << t_offset << "\t" << t_size << endl;

        //trajectory data
        for(deque<Traj>::const_iterator tr=trajectories.begin();tr!=trajectories.end();++tr)
            output << *tr << endl;

        output.close();
    }
    else cerr << " cannot open the file";
}

/**
    \brief export the trajectory data and the position files
    \param filename The name of the trajectory file
    \param base_name A relative path to a file of the positions files serie, like "myDynamicParticles/series05_t000.dat"
    \param token A substring of base_name just before the time step number, like "_t"
*/
void DynamicParticles::saveAll(const string &filename,const string &base_name,const string &token) const
{
    //save trajectories data
    save(filename,base_name,token,0,getNbTimeSteps());
    //extract the path from the filename
    const size_t endfolder = filename.find_last_of("/\\");
    const string folder = filename.substr(0,endfolder+1);

    cout << "writting positions data"<<endl;
    vector<string> tokens(1,token);
    TokenTree tt(tokens,folder+base_name);
    cout << tt << endl;

    vector<size_t> v(1,0);
    while(v[0]<getNbTimeSteps())
    {
        positions[v[0]].exportToFile(tt(v));
        v[0]++;
    }
    cout << "Positions data written" << endl;
}

/** \brief Links the particles of a new frame to existing trajectories, creating trajectories if necessary */
void DynamicParticles::push_back(Particles &parts)
{
    const size_t maxtime = getNbTimeSteps()-1;
    const double range = 2.0*radius;

    //if the Particles object is not indexed, use a R*Tree spatial index
    if(!parts.hasIndex())
        parts.makeRTreeIndex();

    //get possible followers of each trajectory sorted by euclidian distance squared
    vector< multimap<double,size_t> > followersByDist(trajectories.size());
    for(size_t tr = 0;tr<trajectories.size();++tr)
        if(trajectories[tr].exist(maxtime))
            followersByDist[tr] = parts.getEuclidianNeighboursBySqDist((*this)(tr,maxtime),range);

    //link the new positions to existing trajectories, or create new trajectories starting at maxtime+1 if needed
    trajectories.addTimeStep(maxtime,parts.size(),followersByDist);

    //add the frame to the positions
    positions.push_back(&parts);
}

/** @brief export Structured Cell Data Format
	\param filenameFLD The main file to export to
	\param postfix The postfix to append at the end of the labl files
	\param labels Each frame label
	\param t0 Time offset
	\param nsteps Number of steps
	\param stepSize Every N time steps

	The coordinate files are used
*/
void DynamicParticles::exportToFLD(const std::string &postfix,const std::vector<std::map<size_t,double> > &labels,const size_t &stepSize,const double &threshold) const
{
	if(labels.empty())
        throw invalid_argument("DynamicParticles::exportToFLD : empty label");
    if(stepSize*labels.size()>getNbTimeSteps())
    {
        cerr<<"stepSize="<<stepSize<<"\tnsteps="<<labels.size()<<"\tsize="<<getNbTimeSteps()<<endl;
        throw invalid_argument("DynamicParticles::exportToFLD : stepSize*nsteps>size()");
    }
    const size_t halfInterval = stepSize/2;


	//determination of the name patterns
	const string path = trajectories.tt.getPath();
	const string filenameFLD = trajectories.tt.getNoIndexNameNoExt()+postfix+".fld";

	const string positionFilesPattern = trajectories.tt.getPatternNoPath();
	TokenTree positionSerie(trajectories.tt.tokens,positionFilesPattern);
	const string labelFilesPattern = positionFilesPattern.substr(0,positionFilesPattern.find_last_of("."))+postfix+".label";
	TokenTree labelSerie(trajectories.tt.tokens,labelFilesPattern);

	//opening the main file
	ofstream fld(filenameFLD.c_str(), ios::out | ios::trunc );
    if(!fld)
		throw invalid_argument("No such file as "+filenameFLD);

	//header
	const size_t dim = getMaxSimultaneousParticles();
	fld << "# AVS field file \n"
		"# This is a sample for field 1D scalar 3-space irregular\n"
		"\n"
		"nstep = "<<labels.size()<<"\n"
		"ndim = 1\n" //only one set of data
		"dim1 = "<<dim<<"\n"
		"nspace = 3\n" //3D
		"veclen = 2\n" //two one data fields (size and q6)
		"data = double\n"
		"field = irregular\n"; //no grid

	//bounding box (optional)
	const BoundingBox bb = getMaxBox();
	fld << "max_ext =";
	for(int d=0;d<3;++d)
		fld<<" "<<bb.edges[d].second;
	fld<<"\n";
	fld << "min_ext =";
	for(int d=0;d<3;++d)
		fld<<" "<<bb.edges[d].first;

	fld<<endl;


	//data export
	//boost::format varFormat("variable %1%  file=./%2%  filetype=ascii  skip=1  offset=%3%  stride=2\n");
	//boost::format coordFormat("coord %1%  file=./%2%  filetype=ascii  skip=2  offset=%3%  stride=3\n");
	vector<size_t> v(1,halfInterval);
	//boost::progress_display show_progress(nsteps);
	//cout<<"t from 0 to "<<labels.size()<<endl;
	//cout<<"v[0] from "<<halfInterval<<" to "<<labels.size()*stepSize + halfInterval<<endl;
	for(size_t t=0;t<labels.size() && v[0]<getNbTimeSteps();++t)
    {
    	//cout<<t<<"\t";
    	//v[0] = t*stepSize + halfInterval;
    	//open the label file
    	ofstream labelFile((path+labelSerie(v)).c_str(), ios::out | ios::trunc );
		if(!labelFile)
		{
			cerr<<labelSerie<<endl;
			throw invalid_argument("No such file as "+path+labelSerie(v));
		}
    	//fill with label data
    	labelFile << "diameter\tlabel"<<endl;

    	size_t p=0;
    	for(map<size_t,double>::const_iterator l = labels[t].begin();l!=labels[t].end();++l)
    	{
    		while(p<(*l).first)
    		{
				labelFile<<0<<"\t"<<0<<endl;
				p++;
    		}
    		p++;
    		if((*l).second>threshold)
				labelFile<<2.0*radius<<"\t"<<(*l).second<<endl;
			else
				labelFile<<0<<"\t"<<(*l).second<<endl;
    	}
    	while(p<dim)
    	{
    		labelFile<<0<<"\t"<<0<<endl;
			p++;
		}
    	//close label file
    	labelFile.close();
    	//write the reference to the data in the main file
    	/*fld<<varFormat %1 % labelSerie(v) %1;
    	fld<<varFormat %2 % labelSerie(v) %0;
    	for(size_t d=0;d<3;++d)
			fld<<(coordFormat %(d+1) % labelSerie(v) %d);*/
    	fld<<"variable 1  file=./"<<labelSerie(v)<<"  filetype=ascii  skip=1  offset=1  stride=2"<<endl;
    	fld<<"variable 2  file=./"<<labelSerie(v)<<"  filetype=ascii  skip=1  offset=0  stride=2"<<endl;
    	fld<<"coord    1  file=./"<<positionSerie(v)<<"  filetype=ascii  skip=2  offset=0  stride=3"<<endl;
    	fld<<"coord    2  file=./"<<positionSerie(v)<<"  filetype=ascii  skip=2  offset=1  stride=3"<<endl;
    	fld<<"coord    3  file=./"<<positionSerie(v)<<"  filetype=ascii  skip=2  offset=2  stride=3"<<endl;
    	fld<<"EOT"<<endl;
    	//++show_progress;
    	v[0]+=stepSize;
    }
	fld.close();
}

/** @brief export positions, bonds, N time depemdant scalar fields and M time depemdant vector fields to a VTK file serie (Polydata)
	Fields are instantaneous, linking position number to values NOT trajectory numbers to values
  */
void DynamicParticles::exportToVTK(
	const std::vector< scalarDynamicField > &scalars,
	const std::vector< vectorDynamicField > &vectors,
	const size_t &stepSize,
	const boost::ptr_vector< std::vector< std::set<size_t> > > &ngbList
) const
{
	if(scalars.front().second->empty() && vectors.front().second->empty())
        throw invalid_argument("DynamicParticles::exportToVTK : neither scalar field nor vector field");
	/*if(scalars.front().second->size() != vectors.front().second->size())
		throw invalid_argument("DynamicParticles::exportToVTK : incoherent time length between scalar field and vector field");*/
	size_t maxtime;
	if(!scalars.front().second->empty() && !vectors.front().second->empty())
		maxtime = min(scalars.front().second->size(),vectors.front().second->size());
	else
		maxtime = max(scalars.front().second->size(),vectors.front().second->size());

    if(stepSize*maxtime>getNbTimeSteps())
    {
        cerr<<"stepSize="<<stepSize<<"\tnsteps="<<maxtime<<"\tsize="<<getNbTimeSteps()<<endl;
        throw invalid_argument("DynamicParticles::exportToVTK : stepSize*nsteps>size()");
    }
    const size_t halfInterval = stepSize/2;


	//determination of the name patterns
	const string positionFilesPattern = trajectories.tt.getPattern("_dynamic");
	const string vtkFilesPattern = positionFilesPattern.substr(0,positionFilesPattern.find_last_of("."))+".vtk";
	TokenTree vtkSerie(trajectories.tt.tokens,vtkFilesPattern);

	//data export
	size_t avt;
	for(size_t t=0;t<maxtime;t++)
	{
		avt = halfInterval+t*stepSize;
		//convert fields to position dependant
		vector<scalarField> sc(scalars.size());
		for(size_t s=0;s<sc.size();++s)
		{
			sc[s].first = &scalars[s].first;
			sc[s].second = new map<size_t, double>();
			if(!scalars[s].second->empty())
				trajectories.trajToPos(avt,(*scalars[s].second)[t],*sc[s].second);
		}
		vector<vectorField> ve(vectors.size());
		for(size_t v=0;v<ve.size();++v)
		{
			ve[v].first = &vectors[v].first;
			ve[v].second = new map<size_t, Coord>();
			if(!vectors[v].second->empty())
				trajectories.trajToPos(avt,vectors[v].second->at(t),*ve[v].second);
		}
		if(avt<ngbList.size())
		{
			//conversion of the neigbour list into bonds
			deque< pair<size_t,size_t> > bonds;
			for(size_t p=0;p<ngbList[avt].size();++p)
                transform(
                    ngbList[avt][p].lower_bound(p+1),ngbList[avt][p].end(),
                    back_inserter(bonds),
                    boost::bind(make_pair<size_t,size_t>,p,_1)
                );
				/*for(set<size_t>::const_iterator q = ngbList[avt][p].lower_bound(p+1);q!=ngbList[avt][p].end();++q)
					bonds.push_back(make_pair(p,*q));*/

			positions[avt].exportToVTK((vtkSerie % t).str(),bonds,sc,ve);
		}
		else
			positions[avt].exportToVTK((vtkSerie % t).str(),sc,ve);
	}
}


/**
    \brief get the index of the trajectories spanning from t0 to t1 and enclosed inside a given bounding box
*/
set<size_t> DynamicParticles::getSpanning_Enclosed(const TimeBox &b) const
{
    if(!this->hasIndex()) throw logic_error("Set a spatio-temporal index before doing spatio-temporal queries !");
    return (*(this->index))(b);
}

/**
    \brief get the index of the trajectories enclosed inside a given bounding box
    \param b search range
    \return list of the index
*/
set<size_t> DynamicParticles::getEnclosed(const BoundingBox &b) const
{
    if(!this->hasIndex()) throw logic_error("Set a spatio-temporal index before doing spatio-temporal queries !");
    return (*(this->index))(b);
}

/** \brief index of trajectories spanning from t0 to t1 */
set<size_t> DynamicParticles::getSpanning(const Interval &in) const
{
    if(!this->hasIndex()) throw logic_error("Set a spatio-temporal index before doing spatio-temporal queries !");
    return (*(this->index))(in);
}

/**
    \brief get the index of the trajectories enclosed inside a reduction of the minimum bounding box
*/
set<size_t> DynamicParticles::getSpanningInside(const Interval &in,const double &margin) const
{
    if(!this->hasIndex()) throw logic_error("Set a spatio-temporal index before doing spatio-temporal queries !");
    return this->index->getSpanningInside(in, margin);
}

/** \brief get the difference vector between two positions */
Coord DynamicParticles::getDiff(const size_t &tr_from,const size_t &t_from,const size_t &tr_to,const size_t &t_to) const
{
    Coord diff;
    diff = (*this)(tr_to,t_to)-(*this)(tr_from,t_from);
    return diff;
}

/** \brief overall drift between t0 and t1 */
Coord DynamicParticles::getDrift(const set<size_t>&selection,const size_t &t0,const size_t &t1) const
{
    Coord drift(0.0,3);
    for(set<size_t>::iterator tr=selection.begin();tr!=selection.end();++tr)
        drift += getDiff(*tr,t0,*tr,t1);

    drift/=(double)selection.size();
    return drift;
}

/** \brief remove the overall drift between each time step */
void DynamicParticles::removeDrift()
{
    vector<Coord> drifts(positions.size(), Coord(0.0,3));
    Coord maxNegativeDrift(0.0,3);
    for(size_t t0=0;t0+1<getNbTimeSteps();++t0)
    {
        drifts[t0+1] = drifts[t0] - getDrift(getSpanningInside(Interval(t0,t0+1), 2.0*radius),t0,t0+1);
        for(size_t i=0;i<3;++i)
            if(drifts[t0+1][i] < maxNegativeDrift[i])
                maxNegativeDrift[i] = drifts[t0+1][i];
    }
    //the smallest value for origin coordinates is set to 0
    Coord dr;
    for(size_t t0=0;t0<getNbTimeSteps();++t0)
    {
        dr = drifts[t0] - maxNegativeDrift;
        positions[t0] += dr;
    }

	//STindex is now completely wrong and has to be made anew
	this->index = 0;
}

/**    \brief Sum of the square displacement of particles referenced in selection between t0 and t1 */
double DynamicParticles::getSD(const set<size_t>&selection,const size_t &t0,const size_t &t1) const
{
    double result=0.0;
    Coord diff;

    for(set<size_t>::iterator tr=selection.begin();tr!=selection.end();++tr)
    {
        diff = getDiff(*tr,t0,*tr,t1);
        result += dot(diff, diff);
    }
    //cout << result << endl;
    return result;
}
/** \brief Mean square displacement function of time between t0 and t1 for a selection of trajectories */
vector<double> DynamicParticles::getMSD(const set<size_t> &selection,const size_t &t0,const size_t &t1,const size_t &t3) const
{
    const size_t nb_selection = selection.size();
    vector<double> sumSD(t1-t0+1,0.0), nbSD(t1-t0+1,0.0);
    if(selection.empty())
        return sumSD;

    nbSD[0]=1.0;
    if(t3==0)
    {
		for(size_t start=t0;start<t1;++start)
			for(size_t stop=start+1;stop<=t1;++stop)
			{
				sumSD[stop-start] += getSD(selection,start,stop);
				nbSD[stop-start] += 1.0;
			}
		for(size_t t=0;t<sumSD.size();++t)
			sumSD[t]/=nbSD[t]* pow(2.0*radius,2.0)*selection.size();
    }
    else
    {
    	for(size_t Dt=1;Dt<sumSD.size();++Dt)
    	{
			for(size_t start=0;start<t3;++start)
				sumSD[Dt] += getSD(selection,start,start+Dt);
			sumSD[Dt] /= t3 * nb_selection;
    	}
    }
    return sumSD;
}

/** \brief Mean square displacement function of time between t0 and t1 */
vector<double> DynamicParticles::getMSD(const size_t &t0,const size_t &t1,const size_t &t3) const
{
    return getMSD(getSpanning(Interval(t0,t1+t3)),t0,t1,t3);
}

/** \brief Intermediate scatering function of time between t0 and t1 for a selection of trajectories */
vector<double> DynamicParticles::getISF(const set<size_t> &selection,const Coord &q,const size_t &t0,const size_t &t1) const
{
    vector<double> sumISF(t1-t0+1,0.0), nb_per_interval(t1-t0+1,0.0);
    if(selection.empty())
        return sumISF;

    //fill in the basic data used for calculation
    //boost::progress_display show_progress(2*(t1-t0));
    valarray<double> A(0.0,t1-t0+1),B(0.0,t1-t0+1);
    double innerProd;
    for(size_t t=t0;t<=t1;++t)
    {
        for(set<size_t>::iterator tr=selection.begin();tr!=selection.end();++tr)
        {
            innerProd = dot((*this)(*tr,t),q);
            A[t]+=cos(innerProd);
            B[t]+=sin(innerProd);
        }
        //++show_progress;
    }

    //cout << "valarray created" << endl;
    nb_per_interval[0]=1.0;
    for(size_t start=t0;start<t1;++start)
    {
        for(size_t stop=start+1;stop<=t1;++stop)
        {
            for(set<size_t>::iterator tr=selection.begin();tr!=selection.end();++tr)
                sumISF[stop-start] += A[stop]*A[start]+B[stop]*B[start];
            nb_per_interval[stop-start] += 1.0;
        }
        //++show_progress;
    }
    for(size_t t=0;t<sumISF.size();++t)
        sumISF[t]/=nb_per_interval[t]*selection.size();
    sumISF[0]=1.0;
    return sumISF;
}

/** \brief Intermediate scatering function of time between t0 and t1 */
vector<double> DynamicParticles::getISF(const Coord &q,const size_t &t0,const size_t &t1) const
{
    return getISF(getSpanning(Interval(t0,t1)),q,t0,t1);
}

/**
	\brief Self part of intermediate scatering function
	\param selection The trajectories to consider. They should at least span [t0,t1+t3]
	\param q Wave vector
	\param t0 Begining
	\param t1 End
	\param t3 Averaging interval
	\return The self ISF as a t1-t0+1 sized vector of double
	If t3 is 0 (Default), the calculation will act greedily, averaging over all the avilabe intervals of a given length.
		Example : t0=1 t1=4 t3=0
			ISF[0] = 1
			ISF[1] = ( isf([1,2]) + isf([2,3]) + isf([3,4]))/3
			ISF[2] = ( isf([1,3]) + isf([2,4]) )/2
			ISF[3] = isf([1,4])
	If t3>0, the average will be done over t3 time intervals starting from t0
		Example : t0=1 t1=4 t3=2
			ISF[0] = 1
			ISF[1] = ( isf([1,2]) + isf([2,3]) )/2
			ISF[2] = ( isf([1,3]) + isf([2,4]) )/2
			ISF[3] = ( isf([1,4]) + isf([2,5]) )/2
*/
vector<double> DynamicParticles::getSelfISF(const set<size_t> &selection,const Coord &q,const size_t &t0,const size_t &t1,const size_t &t3) const
{
    const size_t nb_selection = selection.size();
    //boost::progress_display show_progress(2*(t1-t0));

    //fill in the basic data used for calculation
    vector< vector<double> > A(t1+t3-t0+1,vector<double>(nb_selection,0.0)),B=A;
    double innerProd;
    size_t p;
    for(size_t t=0;t+t0<=t1+t3;++t)
    {
        p=0;
        for(set<size_t>::iterator tr=selection.begin();tr!=selection.end();++tr)
        {
            innerProd = dot((*this)(*tr,t+t0), q);
            A[t][p]=cos(innerProd);
            B[t][p]=sin(innerProd);
            p++;
        }
        //++show_progress;
    }
    vector<double> sumISF(t1-t0+1,0.0);
    if(t3==0)
    {
		vector<double> nb_per_interval(sumISF.size(),0.0);
		nb_per_interval[0]=1.0;

		for(size_t start=t0;start<t1;++start)
			for(size_t stop=start+1;stop<=t1;++stop)
			{
				for(size_t p=0;p<nb_selection;++p)
					sumISF[stop-start] += A[stop-t0][p]*A[start-t0][p]+B[stop-t0][p]*B[start-t0][p];
				nb_per_interval[stop-start] += 1.0;
			}
		for(size_t t=0;t<sumISF.size();++t)
			sumISF[t]/=nb_per_interval[t]*nb_selection;
    }
    else
    {
    	for(size_t Dt=1;Dt<sumISF.size();++Dt)
    	{
			for(size_t start=0;start<t3;++start)
				for(size_t p=0;p<nb_selection;++p)
					sumISF[Dt] += A[start+Dt][p]*A[start][p]+B[start+Dt][p]*B[start][p];
			sumISF[Dt] /= t3 * nb_selection;
    	}
    }
    sumISF[0]=1.0;
    return sumISF;
}

/** \brief Self part of intermediate scatering function of time between t0 and t1 */
vector<double> DynamicParticles::getSelfISF(const Coord &q,const size_t &t0,const size_t &t1,const size_t &t3) const
{
    return getSelfISF(getSpanning(Interval(t0,t1+t3)),q,t0,t1,t3);
}

/** \brief Get Self ISF averaged over the three axis */
vector<double> DynamicParticles::getSelfISF(const size_t &t0,const size_t &t1,const size_t &t3) const
{
	set<size_t> sp = getSpanning(Interval(t0,t1+t3));
	vector< vector<double> >ISF(4,vector<double>(t1-t0+1));
	vector< Coord > q(3, Coord(0.0,3));
    for(size_t d=0;d<3;++d)
        q[d][d] = M_PI/radius;

	for(size_t d=0;d<3;++d)
		ISF[d] = getSelfISF(sp, q[d], t0, t1, t3);
	for(size_t t=0;t<ISF[3].size();++t)
	{
		for(size_t d=0;d<3;++d)
			ISF[3][t] += ISF[d][t];
		ISF[3][t] /= 3.0;
	}
	return ISF[3];
}

/**
    @brief make MSD and Self ISF for various sets of trajectories
  */
void DynamicParticles::makeDynamics(const std::vector< std::set<size_t> >&sets,std::vector< std::vector<double> > &MSD, std::vector< std::vector<double> > &ISF) const
{
	const size_t stop = getNbTimeSteps()-1;
	MSD.assign(sets.size(),vector<double>(stop+1));
	ISF.assign(sets.size()*4,vector<double>(stop+1));
	cout << "get Mean Square Displacement and Self Intermediate Scattering Function"<<endl;

	vector<Coord> q(3, Coord(0.0,3));
    for(size_t d=0;d<3;++d)
        q[d][d] = M_PI/radius;

	for(size_t s=0;s<sets.size();++s)
	{
        MSD[s] = getMSD(sets[s],0,stop);
        for(size_t d=0;d<3;++d)
        	ISF[4*s+d] = getSelfISF(sets[s],q[d],0,stop);
		for(size_t t=0;t<ISF[4*s+3].size();++t)
		{
			for(size_t d=0;d<3;++d)
                ISF[4*s+3][t] += ISF[4*s+d][t];
            ISF[4*s+3][t] /= 3.0;
		}
	}
}

/** @brief make MSD and Self ISF (along x,y,z + average) for the set of all spanning trajectories */
void DynamicParticles::makeDynamics(vector<double> &MSD,vector< vector<double> > &ISF) const
{
	vector< set<size_t> > sets(1, getSpanning(Interval(0, getNbTimeSteps()-1)));
	vector< vector<double> > MSDs;
	makeDynamics(sets,MSDs,ISF);
	MSD.swap(MSDs.front());
}

/** @brief make MSD and Self ISF (average only) for the set of all spanning trajectories */
void DynamicParticles::makeDynamics(vector<double> &MSD,vector<double> &ISF) const
{
	vector< set<size_t> > sets(1, getSpanning(Interval(0, getNbTimeSteps()-1)));
	vector< vector<double> > MSDs, ISFs;
	makeDynamics(sets,MSDs,ISFs);
	MSD.swap(MSDs.front());
	ISF.swap(ISFs.back());
}

/**
    @brief make and export MSD and Self ISF for various sets of trajectories
  */
void DynamicParticles::exportDynamics(const std::vector< std::set<size_t> >&sets,const std::vector<std::string>&setsNames,const std::string &inputPath) const
{
    vector< vector<double> > MSD;
    vector< vector<double> > ISF;
    makeDynamics(sets,MSD,ISF);

    ostringstream headMSD,headISF;
    string xyz[3] = {"x","y","z"};

    for(size_t s=0;s<sets.size();++s)
    {
        headMSD <<"\t"<<setsNames[s];
        for(size_t d=0;d<3;++d)
			headISF<<"\t"<<setsNames[s]<<"_"<<xyz[d];
		headISF<<"\t"<<setsNames[s]<<"_av";
    }

    saveTable(MSD.begin(),MSD.end(),inputPath + ".msd","t"+headMSD.str(),dt);
    saveTable(ISF.begin(),ISF.end(),inputPath + ".isf","t"+headISF.str(),dt);
}

/** @brief make and export MSD and Self ISF  */
void DynamicParticles::exportDynamics(const string &inputPath) const
{
    vector< set<size_t> > sets(1, getSpanning(Interval(0, getNbTimeSteps()-1)));
    vector<string> setsNames(1,"");
    exportDynamics(sets,setsNames,inputPath);
}

/** @brief averageVelocities
  *	\param
  * @todo: document this function
  */
/*vectorDynamicField DynamicParticles::averageVelocities(const std::set<size_t> &selection,const size_t &displInterval,const size_t &avgInterval) const
{

	if(displInterval==0 || displInterval>=getNbTimeSteps())
	{
		cout<<"\nInvalid displInterval value "<<displInterval<<endl;
		return make_pair("velocities", vectorDynamicField().second);
	}

	const double num = dt*displInterval/radius;

	//case where the full time range is averaged down to one frame
	if(avgInterval==getNbTimeSteps())
	{
		vectorDynamicField avgVel = make_pair("velocities",new vector< map<size_t, Coord> >(1));
		for(set<size_t>::const_iterator tr = selection.begin();tr!=selection.end();++tr)
			avgVel.second->front().insert(
				avgVel.second->front().end(),
				make_pair(
					*tr,
					getDiff(*tr,0,*tr,getNbTimeSteps()-1)/num
				)
			);
		return avgVel;
	}

	//get velocity at each position of each trajectory of the selection
	vector< map<size_t, Coord > > vel(getNbTimeSteps());//-displInterval);
	for(set<size_t>::const_iterator tr = selection.begin();tr!=selection.end();++tr)
	{
		const size_t t_max = min(getNbTimeSteps(),trajectories[*tr].last_time()+1);
		for(size_t t=max((size_t)0,trajectories[*tr].start_time);
					t<t_max-displInterval;
					++t)
			vel[t].insert(
				vel[t].end(),
				make_pair(
					trajectories[*tr][t],
					getDiff(*tr,t,*tr,t+displInterval)/num
				)
			);

		for(size_t t=t_max-displInterval;t<t_max;++t)
			vel[t].insert(
				vel[t].end(),
				make_pair(
					trajectories[*tr][t],
					getDiff(*tr,t-displInterval,*tr,t)/num
				)
			);
	}
	//time average
	vectorDynamicField avgVel = make_pair("velocities",new vector< map<size_t, Coord> >());
	makeSlidingTimeAverage(selection,avgInterval,vel,*avgVel.second);

	if(avgVel.second->empty())
		avgVel.second->push_back(vel.back());

	return avgVel;

	//change the indices form them to correspond to particle positions
	/*vectorDynamicField posVel = make_pair("velocities",vector< map<size_t, valarray<double> > >());
	const size_t halfinterval = avgInterval/2;
	for(size_t t=0;t<avgVel.second.size();++t)
		trajectories.trajToPos(t+halfInterval,avgVel.second[t],posVel.second[t]);

	return posVel;*/
//}*/

/** @brief get the neighbours lost between t_from and t_to by the trajectory tr  */
set<size_t> DynamicParticles::getLostNgbs(const size_t &tr,const size_t &t_from,const size_t &t_to) const
{
	//convert the position index of the neighbours in time t_from to trajectory index
	set<size_t> ngb_from, ngb_to, ngb_diff;
	const size_t p_from = trajectories[tr][t_from];
	transform(
		positions[t_from].getNgbList()[p_from].begin(),
		positions[t_from].getNgbList()[p_from].end(),
		inserter(ngb_from, ngb_from.end()),
		TrajIndex::Inverser(t_from, trajectories)
		);
	//same for t_to
	const size_t p_to = trajectories[tr][t_to];
	transform(
		positions[t_to].getNgbList()[p_to].begin(),
		positions[t_to].getNgbList()[p_to].end(),
		inserter(ngb_to,ngb_to.end()),
		TrajIndex::Inverser(t_to,trajectories)
		);
	//make the difference
	set_difference(
		ngb_from.begin(),ngb_from.end(),
		ngb_to.begin(),ngb_to.end(),
		inserter(ngb_diff,ngb_diff.end())
		);
	return ngb_diff;
}

/** @brief get at each time step the number of Lost Neigbours during interval  */
scalarDynamicField DynamicParticles::getNbLostNgbs(const std::set<size_t> &selection, const size_t &interval) const
{
	if(interval<1 || interval>getNbTimeSteps())
	{
		cout<<"\nInvalid interval value "<<interval<<endl;
		return make_pair("LostNgbs", new vector< map<size_t, double> >());
	}

	scalarDynamicField res;
	res.first = "LostNgbs";
	res.second = new vector< map<size_t,double> >(getNbTimeSteps()-interval+1);

	for(set<size_t>::const_iterator tr=selection.begin();tr!=selection.end();++tr)
		for(size_t t=0;t<res.second->size();++t)
			res.second->at(t).insert(
				res.second->at(t).end(),
				make_pair(
					*tr,
					(double)getLostNgbs(*tr,t,t+interval-1).size()
				)
			);
	return res;
}

/** \brief  Time averaged bond angle distribution */
/*boost::array<double,180> DynamicParticles::getMeanAngularDistribution(const DynNgbList &selection) const
{
    boost::array<double,180> angD;
    fill(angD.begin(), angD.end(), 0.0);
    for(size_t t=0;t<selection.size();++t)
        transform(
            angD.begin(), angD.end(),
            positions[t].getMeanAngularDistribution(selection[t]).begin(), angD.begin(),
            plus<double>()
            );

    transform(
        angD.begin(), angD.end(),
        angD.begin(),
        bind2nd(divides<double>(), (double)selection.size())
        );

    return angD;
}*/



/** \brief import q4 and q6 from file, return spanning trajectories */
/*set<size_t> DynamicParticles::getBooFromFile(const string &prefix,vector< map< size_t, tvmet::Vector<double, 4> > >&qw) const
{
    //cout<<"initializing the returned containers"<<endl;
    qw.assign(getNbTimeSteps(), map<size_t, tvmet::Vector<double, 4> >());

    set<size_t> ret;

    //cout<<"constructing the series of the filenames giving (p q4 q6) for each frame"<<endl;
    vector<size_t> v(1,trajectories.t_offset);
    const string positionFile = trajectories.tt(v);
    const string base = positionFile.substr(0,positionFile.find_last_of("."));
    TokenTree cloud_tt(trajectories.tt.tokens,base+prefix+".cloud");

    //cout<<"reading the (p q4 q6) files"<<endl;
    size_t t,tr;
    //double v4,v6;
    //string trash;
    while(v[0]<trajectories.t_size+trajectories.t_offset)
    {
        t = v[0]-trajectories.t_offset;
        try
        {
			positions[t].getBooFromFile(cloud_tt(v), qw[t]);
        }
        catch(invalid_argument &e)
        {
			cerr<<cloud_tt << endl;
            throw;
        }
        v[0]++;
    }
    cout << "Q4,Q6 data read" << endl;

    bool spans;
    //for each position having a value in the first frame
    for(map<size_t, tvmet::Vector<double, 4> >::const_iterator it=qw.front().begin();it!=qw.front().end();++it)
    {
    	//trajectory index of the particle
    	tr = trajectories.inverse.front()[(*it).first];
    	//does this trajectory span till the end ?
    	spans = trajectories[tr].span(0, getNbTimeSteps()-1);
    	t=1;
    	//if so, does this trajectory have a value in each time step ?
    	while(spans && t<getNbTimeSteps())
    	{
			spans = qw[t].count(trajectories[tr][t]);
			t++;
    	}
    	//if so, the trajectory can be added to the set of interesting trajectories
    	if(spans)
			ret.insert(ret.end(),tr);
    }
	return ret;
}*/

/** @brief Calculate local Bond Orientational Order for each trajectory of the selection at the given time step
  */
void DynamicParticles::makeBoo(const size_t &t, const std::set<size_t> &selection, std::map<size_t,BooData> &allBoo) const
{
    allBoo.clear();
    size_t p;
    for(set<size_t>::const_iterator tr=selection.begin();tr!=selection.end();++tr)
        if(trajectories[*tr].exist(t))
        {
            p = trajectories[*tr][t];
            allBoo.insert(allBoo.end(),make_pair(p,positions[t].getBOO(p)));
        }
}

/** @brief Calculate coarse grained Bond Orientational Order for each position of each trajectory of the selection
  */
void DynamicParticles::makeSBoo(const size_t &t, const std::set<size_t> &selection, const std::map<size_t,BooData> &allBoo, std::map<size_t,BooData> &SallBoo) const
{
    SallBoo.clear();
    size_t p;
    map<size_t,BooData>::iterator it;
    for(set<size_t>::const_iterator tr=selection.begin();tr!=selection.end();++tr)
        if(trajectories[*tr].exist(t))
        {
            p = trajectories[*tr][t];
            it = SallBoo.insert(SallBoo.end(),make_pair(p,BooData()));
            set<size_t> EuNgb = positions[t].getEuclidianNeighbours(positions[t][p],1.3*2.0*radius);
            //sum up the contribution of each neighbour including the particle itself.
            for(set<size_t>::iterator n=EuNgb.begin();n!=EuNgb.end();++n)
                (*it).second += (*allBoo.find(*n)).second;

            (*it).second /= (double)(EuNgb.size());
        }
            //SallBoo[t][*tr] = positions[t].getAvBOO(allBoo[t],trajectories[*tr][t],1.3*2.0*radius);
}


/** @brief Average over time a time dependant and trajectory dependant value.
  *
  * \param selection The trajectories to treat
  * \param timeDependant The input values, function of time and of trajectories
  * \param timeAverage The output, function of the trajectories
  */
void DynamicParticles::makeTimeAverage(const std::set<size_t> &selection, const size_t &avgInterval, const std::vector< std::map<size_t,double> > &timeDependant, std::vector< std::map<size_t,double> > &timeAveraged) const
{
    map<size_t,double>::iterator it;
    timeAveraged.assign((size_t)max(getNbTimeSteps()/(double)avgInterval,1.0),map<size_t,double>());
    for(size_t avt=0;avt<timeAveraged.size();++avt)
        for(set<size_t>::const_iterator tr=selection.begin();tr!=selection.end();++tr)
            if(trajectories[*tr].span(avt*avgInterval,(avt+1)*avgInterval-1))
            {
                it = timeAveraged[avt].insert(timeAveraged[avt].end(),std::make_pair(*tr,0.0));
                for(size_t t=avt*avgInterval;t<(avt+1)*avgInterval;++t)
                    (*it).second += (*timeDependant[t].find(trajectories[*tr][t])).second;

                (*it).second /= (double)avgInterval;
            }
}

/** @brief Average over time a time dependant and trajectory dependant value.
  *
  * \param selection The trajectories to treat
  * \param timeDependant The input values, function of time and of trajectories
  * \param timeAverage The output, function of the trajectories
  */
void DynamicParticles::makeSlidingTimeAverage(
	const std::set<size_t> &selection,
	const size_t &avgInterval,
	const std::vector< std::map<size_t,double> > &timeDependant,
	std::vector< std::map<size_t,double> > &timeAveraged
) const
{
    map<size_t,double>::iterator it;
    map<size_t,double>::const_iterator td;
    timeAveraged.assign(timeDependant.size()-(avgInterval-1),std::map<size_t,double>());
    //cout<<timeAveraged.size()<<" steps left"<<endl;
    for(size_t avt=0;avt<timeAveraged.size();++avt)
    {
    	//cout<<"avt="<<avt<<" keep trajectories spanning between "<<avt<<" and "<<avt+avgInterval-1<<endl;
        for(std::set<size_t>::const_iterator tr=selection.begin();tr!=selection.end();++tr)
            if(trajectories[*tr].span(avt,avt+avgInterval-1))
            {
                it = timeAveraged[avt].insert(timeAveraged[avt].end(),std::make_pair(*tr,0.0));
                for(size_t t=avt;t<avt+avgInterval;++t)
				{
					td = timeDependant[t].find(trajectories[*tr][t]);
					if(td == timeDependant[t].end())
					{
						std::cerr<<"avt="<<avt<<"\tt="<<t<<"\ttr="<<*tr<<"\tstart="<<trajectories[*tr].start_time<<"\tlast="<<trajectories[*tr].last_time()<<std::endl;
						/*for(size_t i=trajectories[*tr].start_time;i<=trajectories[*tr].last_time();++i)
							std::cerr<<(*timeDependant[i].find(trajectories[*tr][i])).second<<"\t";
						std::cerr<<std::endl;*/
						throw std::invalid_argument("the trajectory tr has no assigned value at time step t");
					}
                    (*it).second += (*td).second;
				}

                (*it).second /= (double)avgInterval;
            }
    }
}

