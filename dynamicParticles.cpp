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

#include "dynamicParticles.hpp"
#include "files_series.hpp"

#include <boost/progress.hpp>
#include <boost/bind.hpp>

#include <ctime>
#include <numeric>




using namespace std;
using namespace Colloids;
//using namespace tvmet;

/** @brief Constructor from data. Take the ownership of the positions. No linking necessary  */
DynamicParticles::DynamicParticles(const TrajMap &trajs, boost::ptr_vector<Particles>& positions, const double &rad,const double &time_step) :
    trajectories(trajs), radius(rad), dt(time_step)
{
    this->positions.swap(positions);
}

/** @brief Constructor from data. Take the ownership of the positions and link them into trajectories.  */
DynamicParticles::DynamicParticles(boost::ptr_vector<Particles>& positions, const double &rad,const double &time_step, const string &displFile, const size_t &offset) :
    radius(rad), dt(time_step)
{
    this->positions.swap(positions);
    removeDrift(displFile, offset);
    link();
}

/** @brief Constructor from files. Load the positions. No linking necessary  */
DynamicParticles::DynamicParticles(const TrajMap &trajs, FileSerie &files, const double &rad,const double &time_step) :
    trajectories(trajs), radius(rad), dt(time_step)
{
    fill(files);
}

/** @brief Constructor from files. Load the positions and link them into trajectories.  */
DynamicParticles::DynamicParticles(FileSerie &files, const double &rad,const double &time_step) :
    radius(rad), dt(time_step)
{
    fill(files);
    removeDrift(files.head()+".displ", files.get_offset());
    link();
}

/** \brief constructor from file */
DynamicParticles::DynamicParticles(const string &filename)
{
    //extract the path from the filename
    const size_t endfolder = filename.find_last_of("/\\");
    const string folder = filename.substr(0,endfolder+1);
    size_t t_offset, t_size;

    ifstream input(filename.c_str(), ios::in);
    if(!input.good())
        throw invalid_argument((filename+" doesn't exist").c_str() );

    //header
    input >> radius >> dt;
    input.get(); //escape the endl

    //data of the file serie containing the positions
    string base_name,token;
    getline(input,base_name);
    getline(input,token);
    input >> t_offset >> t_size;

    FileSerie files(folder+base_name, token, t_size, t_offset);
    this->fill(files);

    //construct the TrajIndex from file stream
    input >> trajectories;
    trajectories.makeInverse(this->getFrameSizes());
}

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
    vector<size_t> sizes=getFrameSizes();
    return *max_element(sizes.begin(), sizes.end());
}

/** @brief get the number of particle in each time step  */
vector<size_t> DynamicParticles::getFrameSizes() const
{
    vector<size_t> sizes(positions.size());
    transform(positions.begin(), positions.end(), sizes.begin(), mem_fun_ref(&Particles::size));
    return sizes;
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
/* DynamicParticles::DynamicParticles(const double &rad,const double &time_step,const size_t &t_size)
{
    radius = rad;
    dt = time_step;
    positions.resize(t_size,new Particles(0,rad));
}*/

/** \brief constructor from the first time step */
/*DynamicParticles::DynamicParticles(Particles *parts,const double &time_step)
{
    positions.resize(1,parts);
    for(size_t p=0;p<parts->size();++p)
        trajectories.push_back(Traj(0,p));
    dt=time_step;
    radius = parts->radius;
    return;
}*/



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
/*DynamicParticles::DynamicParticles(
    const double &rad, const double &time_step,
    const string &base_name, const string &token,
    const size_t &t_offset, const size_t &t_size
    )
{
    radius = rad;
    dt = time_step;

    vector<string> tokens(1,token);
    TokenTree tt(tokens,base_name);

     //load all time steps
    try
    {
        vector<size_t> v(1,t_offset);
        positions.reserve(t_size);
        while(v[0]<t_offset+t_size)
            positions.push_back(new Particles(tt(v), radius);
    }
    catch(const exception &e)
    {
        //If any error in this loop, it must be due to the TokenTree
        cerr<<tt<<endl;
        throw;
    }

    //look for a displacement file
    ifstream displFile((base_name.substr(0,base_name.rfind(token))+".displ").c_str(), ios::in);
    if(displFile.good())
    {
        vector<Coord> displ(t_size, Coord(0));
        double trash;
        //get the relative displacements form file
		cout<<"Using "<< base_name.substr(0,base_name.rfind(token)) << ".displ as an hint for global displacements between frames"<<endl;
		for(size_t t=0;t<t_offset;++t)
			displFile >> trash >>trash;
		for(size_t t=0;t+t_offset<t_size;++t)
			displFile >> displ[t][0] >> displ[t][1]; //only 2D displacements
		displFile.close();

		//compute objective displacement vector
		partial_sum(displ.begin(),displ.end(),displ.begin());

		//compensate the displacement
        for(size_t t=0;t<t_size;++t)
            positions[t] += displ[t-t_offset];
    }


    return;
}*/

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
        output << trajectories;

        output.close();
    }
    else cerr << " cannot open the file";
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
/*void DynamicParticles::exportToFLD(const std::string &postfix,const std::vector<std::map<size_t,double> > &labels,const size_t &stepSize,const double &threshold) const
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
}*/

/** @brief export positions, bonds, N time depemdant scalar fields and M time depemdant vector fields to a VTK file serie (Polydata)
	Fields are instantaneous, linking position number to values NOT trajectory numbers to values
  */
void DynamicParticles::exportToVTK(
	FileSerie &files,
	std::vector< ScalarDynamicField > &scalars,
	std::vector< VectorDynamicField > &vectors
) const
{
	//determination of the name patterns
	/*const string positionFilesPattern = trajectories.tt.getPattern("_dynamic");
	const string vtkFilesPattern = positionFilesPattern.substr(0,positionFilesPattern.find_last_of("."))+".vtk";
	TokenTree vtkSerie(trajectories.tt.tokens,vtkFilesPattern);*/

	//data export
	for(size_t t=0;t<getNbTimeSteps();t++)
	{
		vector<ScalarField> sc;
		for(size_t s=0; s<scalars.size();++s)
			sc.push_back(scalars[s][t]);

		vector<VectorField> vec;
		for(size_t v=0; v<vectors.size();++v)
			vec.push_back(vectors[v][t]);
		positions[t].exportToVTK(files % t,sc, vec);
	}
}


/**
    \brief get the index of the trajectories spanning from t0 to t1 and enclosed inside a given bounding box
*/
vector<size_t> DynamicParticles::selectSpanning_Enclosed(const TimeBox &b) const
{
    if(!this->hasIndex()) throw logic_error("Set a spatio-temporal index before doing spatio-temporal queries !");
    return (*(this->index))(b);
}

/**
    \brief get the index of the trajectories enclosed inside a given bounding box
    \param b search range
    \return list of the index
*/
vector<size_t> DynamicParticles::selectEnclosed(const BoundingBox &b) const
{
    if(!this->hasIndex()) throw logic_error("Set a spatio-temporal index before doing spatio-temporal queries !");
    return (*(this->index))(b);
}

/** \brief index of trajectories spanning from t0 to t1 */
vector<size_t> DynamicParticles::selectSpanning(const Interval &in) const
{
    if(this->hasIndex())
		return (*(this->index))(in);
	else
	{
		list<size_t> a, b;
    	copy(
			trajectories.inverse[in.first].begin(),
			trajectories.inverse[in.first].end(),
			back_inserter(a)
			);
		a.sort();
		a.unique();
		copy(
			trajectories.inverse[in.second].begin(),
			trajectories.inverse[in.second].end(),
			back_inserter(b)
			);
		b.sort();
		b.unique();
		vector<size_t> c;
		c.reserve(min(a.size(), b.size()));
		set_intersection(
			a.begin(), a.end(),
			b.begin(), b.end(),
			back_inserter(c)
			);
		return c;
	}
}

/**
    \brief get the index of the trajectories enclosed inside a reduction of the minimum bounding box
*/
vector<size_t> DynamicParticles::selectSpanningInside(const Interval &in,const double &margin) const
{
    if(!this->hasIndex()) throw logic_error("Set a spatio-temporal index before doing spatio-temporal queries !");
    return this->index->getSpanningInside(in, margin);
}

/** \brief get the difference vector between two positions */
Coord DynamicParticles::getDiff(const size_t &tr_from,const size_t &t_from,const size_t &tr_to,const size_t &t_to) const
{
    Coord diff(3);
    diff = (*this)(tr_to,t_to)-(*this)(tr_from,t_from);
    return diff;
}

/** \brief overall drift between t0 and t1 */
Coord DynamicParticles::getDrift(const vector<size_t>&selection,const size_t &t0,const size_t &t1) const
{
    Coord drift(0.0,3);
    for(ssize_t tr=0;tr<selection.size();++tr)
        drift += getDiff(selection[tr],t0,selection[tr],t1);

    drift/=(double)selection.size();
    return drift;
}
/** @brief getDrift. With or without indexing  */
Coord DynamicParticles::getDrift(const size_t &t0,const size_t &t1) const
{
	if(hasIndex())
		return getDrift(selectSpanningInside(Interval(t0,t0+1), 2.0*radius),t0,t0+1);
	else
		return getDrift(selectSpanning(Interval(t0,t0+1)), t0, t0+1);
}


/** \brief remove the overall drift between each time step */
void DynamicParticles::removeDrift()
{
    vector<Coord> relative_drifts(positions.size(), Coord(0.0,3)),
            drifts(positions.size(), Coord(0.0,3));

    //#pragma omp parallel for shared(relative_drifts)
    for(ssize_t t=1;t<(ssize_t)getNbTimeSteps();++t)
        relative_drifts[t] = - getDrift(t-1, t);

    partial_sum(
        relative_drifts.begin(), relative_drifts.end(),
        drifts.begin()
        );

    Coord maxNegativeDrift(0.0,3);
    //#pragma omp parallel for shared(drifts, maxNegativeDrift)
    for(ssize_t t=1;t<(ssize_t)getNbTimeSteps();++t)
        for(size_t i=0;i<3;++i)
            if(drifts[t][i] < maxNegativeDrift[i])
                maxNegativeDrift[i] = drifts[t][i];

    //the smallest value for origin coordinates is set to 0
    //#pragma omp parallel for shared(drifts, maxNegativeDrift)
    for(ssize_t t0=0;t0<(ssize_t)positions.size();++t0)
    {
        Coord dr(3);
        dr = drifts[t0] - maxNegativeDrift;
        positions[t0] += dr;
    }

	//STindex is now completely wrong and has to be made anew
	this->index.reset();
}

/** @brief remove drift indicated in displFile  */
void DynamicParticles::removeDrift(const std::string &displFile, const size_t &t_offset)
{
	vector<Coord> displ(getNbTimeSteps(),Coord(0.0,3));
	double trash;
    ifstream f(displFile.c_str(), ios::in);
    if(f.good())
    {
		//get the relative displacements form file
		cout<<"Using "<< displFile << " as an hint for global displacements between frames"<<endl;
		for(size_t t=0;t<t_offset;++t)
			f >> trash >>trash;
		for(size_t t=0;t+t_offset<getNbTimeSteps();++t)
			f >> displ[t][0] >> displ[t][1]; //only 2D displacements
		f.close();

		//compute objective displacement vector
		partial_sum(displ.begin(),displ.end(),displ.begin());
		//largest absolute negative displacement
		valarray<double> maxd(0.0,3);
		for(int t=0;t<(int)displ.size();++t)
		{
			maxd[0] = max(maxd[0],displ[t][0]);
			maxd[1] = max(maxd[1],displ[t][1]);
		}
		transform(
			displ.begin(),displ.end(),displ.begin(),
			bind2nd(minus< valarray<double> >(),maxd)
			);
		for(int t=0;t<(int)displ.size();++t)
			positions[t] += displ[t];
    }
    else
		cout<<"impossible to find "<<displFile<<endl;


}



/**    \brief Sum of the square displacement of particles referenced in selection between t0 and t1 */
double DynamicParticles::getSD(const vector<size_t>&selection,const size_t &t0,const size_t &t1) const
{
    double result=0.0;
    Coord diff(3);

    for(ssize_t tr=0;tr<selection.size();++tr)
    {
        diff = getDiff(selection[tr],t0,selection[tr],t1);
        result += dot(diff, diff);
    }
    //cout << result << endl;
    return result;
}

/** @brief get square displacement of all particles at time t
  *
  * Using centered scheme except at the begining or the end of the trajectory
  */
vector<double> DynamicParticles::getSD(const size_t &t, const size_t &halfInterval) const
{
	vector<Coord> vel = velocities(t, halfInterval);
	vector<double> sd(vel.size());
	transform(
		vel.begin(), vel.end(),
		vel.begin(), sd.begin(),
		dot
		);
	transform(
		sd.begin(), sd.end(), sd.begin(),
		bind2nd(multiplies<double>(), halfInterval*2+1)
		);
	return sd;
}



/** \brief Mean square displacement function of time between t0 and t1 for a selection of trajectories */
vector<double> DynamicParticles::getMSD(const vector<size_t> &selection,const size_t &t0,const size_t &t1,const size_t &t3) const
{
    const size_t nb_selection = selection.size();
    vector<double> sumSD(t1-t0+1,0.0), nbSD(t1-t0+1,0.0);
    if(selection.empty())
        return sumSD;

    nbSD[0]=1.0;
    if(t3==0)
    {
    	#pragma omp parallel for schedule(runtime) shared(sumSD, nbSD, t0, t1, t3, selection)
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
    	    #pragma omp parallel for schedule(runtime) shared(sumSD, t3, Dt)
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
    return getMSD(selectSpanning(Interval(t0,t1+t3)),t0,t1,t3);
}

/** @brief get Non Gaussian parameter function of time between t0 and t1  */
vector<double> DynamicParticles::getNonGaussian(const std::vector<size_t> &selection,const size_t &t0,const size_t &t1,const size_t &t3) const
{
	const size_t nb_selection = selection.size();
    vector<double> alpha(t1-t0+1,0.0), sumSD(t1-t0+1,0.0), sumQD(t1-t0+1,0.0), nb(t1-t0+1,0.0);
    if(selection.empty())
        return sumSD;

    nb[0]=1.0;
    if(t3==0)
    {
    	#pragma omp parallel for schedule(runtime) shared(sumSD, nbSD, t0, t1, t3, selection)
		for(size_t start=t0;start<t1;++start)
			for(size_t stop=start+1;stop<=t1;++stop)
			{
				for(ssize_t tr=0;tr<selection.size();++tr)
				{
					Coord diff(3);
					diff = getDiff(selection[tr],t0,selection[tr],t1);
					const double sd = dot(diff, diff);
					sumSD[stop-start] += sd;
					sumQD[stop-start] += sd*sd;
				}
				nb[stop-start] += 1.0;
			}
		for(size_t t=0;t<sumSD.size();++t)
			alpha[t] = nb[t]*selection.size()/3.0 * sumQD[t]/sumSD[t]/sumSD[t];
    }
    else
    {
    	for(size_t Dt=1;Dt<sumSD.size();++Dt)
    	{
    	    #pragma omp parallel for schedule(runtime) shared(sumSD, t3, Dt)
			for(size_t start=0;start<t3;++start)
				for(ssize_t tr=0;tr<selection.size();++tr)
				{
					Coord diff(3);
					diff = getDiff(selection[tr],t0,selection[tr],t1);
					const double sd = dot(diff, diff);
					sumSD[Dt] += sd;
					sumQD[Dt] += sd*sd;
				}
			alpha[Dt] = (t3*nb_selection)/3.0 * sumQD[Dt]/sumSD[Dt]/sumSD[Dt];
    	}
    }
    return alpha;
}



/** \brief Intermediate scatering function of time between t0 and t1 for a selection of trajectories */
vector<double> DynamicParticles::getISF(const vector<size_t> &selection,const Coord &q,const size_t &t0,const size_t &t1) const
{
    vector<double> sumISF(t1-t0+1,0.0), nb_per_interval(t1-t0+1,0.0);
    if(selection.empty())
        return sumISF;

    //fill in the basic data used for calculation
    //boost::progress_display show_progress(2*(t1-t0));
    valarray<double> A(0.0,t1-t0+1),B(0.0,t1-t0+1);
    double innerProd;
    #pragma omp parallel for schedule(runtime) shared(t0, t1, selection) private(innerProd)
    for(size_t t=t0;t<=t1;++t)
    {
        for(ssize_t tr=0;tr<selection.size();++tr)
        {
            innerProd = dot((*this)(selection[tr],t),q);
            A[t]+=cos(innerProd);
            B[t]+=sin(innerProd);
        }
        //++show_progress;
    }

    //cout << "valarray created" << endl;
    nb_per_interval[0]=1.0;
    #pragma omp parallel for schedule(runtime) shared(A, B, sumISF, t0, t1, selection, nb_per_interval)
    for(size_t start=t0;start<t1;++start)
    {
        for(size_t stop=start+1;stop<=t1;++stop)
        {
            for(ssize_t tr=0;tr<selection.size();++tr)
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
    return getISF(selectSpanning(Interval(t0,t1)),q,t0,t1);
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
vector<double> DynamicParticles::getSelfISF(const vector<size_t> &selection,const Coord &q,const size_t &t0,const size_t &t1,const size_t &t3) const
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
        for(ssize_t tr=0;tr<selection.size();++tr)
        {
            innerProd = dot((*this)(selection[tr],t+t0), q);
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
		for(int t=0;t<(int)sumISF.size();++t)
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
    return getSelfISF(selectSpanning(Interval(t0,t1+t3)),q,t0,t1,t3);
}

/** \brief Get Self ISF averaged over the three axis */
vector<double> DynamicParticles::getSelfISF(const size_t &t0,const size_t &t1,const size_t &t3) const
{
	vector<size_t> sp = selectSpanning(Interval(t0,t1+t3));
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
void DynamicParticles::makeDynamics(const std::vector< std::vector<size_t> >&sets,std::vector< std::vector<double> > &MSD, std::vector< std::vector<double> > &ISF) const
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
	vector< vector<size_t> > sets(1, selectSpanning(Interval(0, getNbTimeSteps()-1)));
	vector< vector<double> > MSDs;
	makeDynamics(sets,MSDs,ISF);
	MSD.swap(MSDs.front());
}

/** @brief make MSD and Self ISF (average only) for the set of all spanning trajectories */
void DynamicParticles::makeDynamics(vector<double> &MSD,vector<double> &ISF) const
{
	vector< vector<size_t> > sets(1, selectSpanning(Interval(0, getNbTimeSteps()-1)));
	vector< vector<double> > MSDs, ISFs;
	makeDynamics(sets,MSDs,ISFs);
	MSD.swap(MSDs.front());
	ISF.swap(ISFs.back());
}

/**
    @brief make and export MSD and Self ISF for various sets of trajectories
  */
void DynamicParticles::exportDynamics(const std::vector< std::vector<size_t> >&sets,const std::vector<std::string>&setsNames,const std::string &inputPath) const
{
    vector< vector<double> > MSD;
    vector< vector<double> > ISF;
    makeDynamics(sets,MSD,ISF);

    string xyz[3] = {"x","y","z"};
    ofstream msd_f((inputPath + ".msd").c_str());
    ofstream isf_f((inputPath + ".isf").c_str());
    msd_f << "#t";
    isf_f << "#t";

    for(size_t s=0;s<sets.size();++s)
    {
        msd_f <<"\t"<<setsNames[s];
        for(size_t d=0;d<3;++d)
			isf_f<<"\t"<<setsNames[s]<<"_"<<xyz[d];
		isf_f<<"\t"<<setsNames[s]<<"_av";
    }

	for(size_t t=0; t<MSD.front().size(); ++t)
	{
		msd_f << t*dt;
		isf_f << t*dt;
		for(size_t s=0;s<sets.size();++s)
		{
			msd_f << "\t"<< MSD[s][t];
			for(size_t d=0; d<4; ++d)
				isf_f<<"\t"<< ISF[4*s+d][t];
		}
		msd_f<<"\n";
		isf_f<<"\n";
	}
}

/** @brief make and export MSD and Self ISF  */
void DynamicParticles::exportDynamics(const string &inputPath) const
{
    vector< vector<size_t> > sets(1, selectSpanning(Interval(0, getNbTimeSteps()-1)));
    vector<string> setsNames(1,"");
    exportDynamics(sets,setsNames,inputPath);
}

/** @brief velocities of every particle at time t

	Using centered scheme except at begining and end of trajectory.
 */
vector<Coord> DynamicParticles::velocities(const size_t &t, const size_t &halfInterval) const
{
	vector<Coord> vel(trajectories.inverse[t].size(), Coord(0.0, 3));
	size_t start, stop;
	for(size_t p=0; p<trajectories.inverse[t].size(); ++p)
	{
		const Traj &tr = trajectories[trajectories.inverse[t][p]];
		if(tr.steps.size()>1)
		{
			start = max((int)(t-halfInterval), (int)tr.start_time);
			stop = min(t+halfInterval, tr.last_time());
			vel[p] = (positions[stop][tr[stop]] - positions[start][tr[start]]) / (double)(stop-start);
		}
	}
	return vel;
}



/** @brief get the neighbours lost between t_from and t_to by the trajectory tr  */
vector<size_t> DynamicParticles::getLostNgbs(const size_t &tr,const size_t &t_from,const size_t &t_to) const
{
	if(t_from==t_to)
		return vector<size_t>();
	//convert the position index of the neighbours in time t_from to trajectory index
	list<size_t> ngb_from, ngb_to;
	const vector<size_t>& n_from = positions[t_from].getNgbList()[trajectories[tr][t_from]];
	for(ssize_t p=0; p<n_from.size(); ++p)
		ngb_from.push_back(trajectories.inverse[t_from][n_from[p]]);
	ngb_from.sort();
	ngb_from.unique();

	//same for t_to
	const vector<size_t>& n_to = positions[t_to].getNgbList()[trajectories[tr][t_to]];
	for(ssize_t p=0; p<n_to.size(); ++p)
		ngb_to.push_back(trajectories.inverse[t_to][n_to[p]]);
	ngb_to.sort();
	ngb_to.unique();

	//make the difference
	vector<size_t> ngb_diff;
	ngb_diff.reserve(ngb_from.size());
	set_difference(
		ngb_from.begin(),ngb_from.end(),
		ngb_to.begin(),ngb_to.end(),
		back_inserter(ngb_diff)
		);
	return ngb_diff;
}

/** @brief Number of lost neighbours of every particle at time t

	Using centered scheme except at begining and end of trajectory.
 */
vector<double> DynamicParticles::getNbLostNgbs(const size_t &t, const size_t &halfInterval) const
{
	vector<double> nb(trajectories.inverse[t].size());
	size_t tr, start, stop;
	for(size_t p=0; p<trajectories.inverse[t].size(); ++p)
	{
		tr = trajectories.inverse[t][p];
		if(trajectories[tr].steps.size()>1)
		{
			start = max((int)(t-halfInterval), (int)trajectories[tr].start_time);
			stop = min(t+halfInterval, trajectories[tr].last_time());
			nb[p] = getLostNgbs(tr, start, stop).size();
		}
	}
	return nb;
}

/** \brief split recursively a succession of positions into well separated regions */
void splitByCageJump(const vector<Coord> &positions, const double &threshold, const size_t &resolution, const Interval &in, list<size_t> &jumps)
{
    vector<double> separation(in.second-in.first+1);
    for(size_t t=in.first; t<in.second; ++t)
    {
        //center of mass of both sub intervals
        Coord c1(3), c2(3);
        c1 = accumulate(positions.begin()+in.first, positions.begin()+t+1, c1)/(double)(t+1-in.first);
        c2 = accumulate(positions.begin()+t+1, positions.begin()+in.second+1, c1)/(double)(in.second-t);
        //distances to the center of mass of the sub interval
        double d1=0.0, d2=0.0;
        for(size_t i=in.first; i<t+1; ++t)
            d1 += norm2(positions[i]-c1);
        d1 /= (double)(t+1-in.first);
        for(size_t i=t+1; i<in.second+1; ++t)
            d2 += norm2(positions[i]-c2);
        d2 /= (double)(in.second-t);

        separation[t-in.first] = sqrt(
            (t-in.first)/(double)in.second*(1.0-(t-in.first)/(double)in.second)
            * d1 * d2
            );
    }
    const size_t tc = (max_element(separation.begin(), separation.end()) - separation.begin()) + in.first;
    //continue only if the maximum separation is larger than the threshold
    if(separation[tc-in.first] >= threshold)
    {
        jumps.push_back(tc);
        Interval in1(in.first, tc), in2(tc+1, in.second);
        if(in1.second-in1.first > resolution)
            splitByCageJump(positions, threshold, resolution, in1, jumps);
        if(in2.second-in2.first > resolution)
            splitByCageJump(positions, threshold, resolution, in2, jumps);
    }
}

/** @brief split trajectories by cages jump
  *
  * \param threshold minimum distance of a jump, typically diameter/10
  * \param resolution minimum number of time steps between two jumps, >0
  */
TrajIndex DynamicParticles::getCages(const size_t &resolution) const
{
    TrajIndex cages;
    const double thrsq = pow(radius*0.2, 2.0);
    for(TrajIndex::const_iterator tr=trajectories.begin(); tr!=trajectories.end(); ++tr)
    {
        //trajectory shorter than the resolution
        if(tr->last_time()+1-tr->start_time < resolution)
        {
            cages.push_back(*tr);
            continue;
        }
        //prepare to split the trajectory by cage jump
        vector<Coord> pos((tr->last_time())+1-(tr->start_time), Coord(3));
        for(size_t t=tr->start_time; t<tr->last_time()+1; ++t)
            pos[t-tr->start_time] = positions[t][(*tr)[t]];
        Interval in(0, pos.size()-1);
        list<size_t> jumps;
        //find the jumps
        splitByCageJump(pos, thrsq, resolution, in, jumps);
        assert(!jumps.empty());
        //split
        jumps.sort();
        assert(jumps.front()>0);
        jumps.push_back(pos.size()-1);

        cages.push_back(tr->subtraj(tr->start_time, jumps.front()+(tr->start_time)-1));
        list<size_t>::const_iterator i=jumps.begin(), j=jumps.begin();
        j++;
        while(j!=jumps.end())
        {
            cages.push_back(tr->subtraj((*i)+tr->start_time, (*j)+tr->start_time-1));
            j++;
            i++;
        }
    }
    //make TrajIndex inverse
    cages.makeInverse(this->getFrameSizes());
    return cages;
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
void DynamicParticles::makeBoo(const size_t &t, const std::vector<size_t> &selection, std::map<size_t,BooData> &allBoo) const
{
    allBoo.clear();
    size_t p;
    for(ssize_t tr=0;tr<selection.size();++tr)
        if(trajectories[selection[tr]].exist(t))
        {
            p = trajectories[selection[tr]][t];
            allBoo.insert(allBoo.end(),make_pair(p,positions[t].getBOO(p)));
        }
}

/** @brief Calculate coarse grained Bond Orientational Order for each position of each trajectory of the selection
  */
void DynamicParticles::makeSBoo(const size_t &t, const std::vector<size_t> &selection, const std::map<size_t,BooData> &allBoo, std::map<size_t,BooData> &SallBoo) const
{
    SallBoo.clear();
    size_t p;
    map<size_t,BooData>::iterator it;
    for(ssize_t tr=0;tr<selection.size();++tr)
        if(trajectories[selection[tr]].exist(t))
        {
            p = trajectories[selection[tr]][t];
            it = SallBoo.insert(SallBoo.end(),make_pair(p,BooData()));
            vector<size_t> EuNgb = positions[t].getEuclidianNeighbours(positions[t][p],1.3*2.0*radius);
            //sum up the contribution of each neighbour including the particle itself.
            for(ssize_t n=0;n<EuNgb.size();++n)
                (*it).second += (*allBoo.find(EuNgb[n])).second;

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
void DynamicParticles::makeTimeAverage(const std::vector<size_t> &selection, const size_t &avgInterval, const std::vector< std::map<size_t,double> > &timeDependant, std::vector< std::map<size_t,double> > &timeAveraged) const
{
    map<size_t,double>::iterator it;
    timeAveraged.assign((size_t)max(getNbTimeSteps()/(double)avgInterval,1.0),map<size_t,double>());
    for(size_t avt=0;avt<timeAveraged.size();++avt)
        for(ssize_t tr=0;tr<selection.size();++tr)
            if(trajectories[selection[tr]].span(avt*avgInterval,(avt+1)*avgInterval-1))
            {
                it = timeAveraged[avt].insert(timeAveraged[avt].end(),std::make_pair(selection[tr],0.0));
                for(size_t t=avt*avgInterval;t<(avt+1)*avgInterval;++t)
                    (*it).second += (*timeDependant[t].find(trajectories[selection[tr]][t])).second;

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
	const std::vector<size_t> &selection,
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
        for(ssize_t tr=0;tr<selection.size();++tr)
            if(trajectories[selection[tr]].span(avt,avt+avgInterval-1))
            {
                it = timeAveraged[avt].insert(timeAveraged[avt].end(),std::make_pair(selection[tr],0.0));
                for(size_t t=avt;t<avt+avgInterval;++t)
				{
					td = timeDependant[t].find(trajectories[selection[tr]][t]);
					if(td == timeDependant[t].end())
					{
						std::cerr<<"avt="<<avt<<"\tt="<<t<<"\ttr="<<selection[tr]<<"\tstart="<<trajectories[selection[tr]].start_time<<"\tlast="<<trajectories[selection[tr]].last_time()<<std::endl;
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

/** @brief load each file into the positions  */
void DynamicParticles::fill(FileSerie &files)
{
    positions.reserve(files.size());
    for(size_t t=0; t<files.size();++t)
        positions.push_back(new Particles(files%t, radius));
}

/** @brief link positions into trajectories  */
void DynamicParticles::link()
{
	const double range = this->radius * 2.0;
    //spatially index each unindexed frame by a RTreeIndex. Needed for the linking
    cout<<"index"<<endl;
    #pragma omp parallel for schedule(runtime)
    for(size_t t=0; t<positions.size(); ++t)
        if(!positions[t].hasIndex())
            positions[t].makeRTreeIndex();

    //link the positions into trajectories
    boost::progress_display show_progress(positions.size()-1);
    double Error=0, maxError=0, sumError=0;
    TrajMap tm(positions[0].size());
    for(size_t t=0; t<positions.size()-1; ++t)
    {
        size_t nbTraj = tm.getNbTraj();
        vector< multimap<double,size_t> > followersByDist(positions[t].size());

        #pragma omp parallel for schedule(runtime) shared(t, followersByDist)
        for(ssize_t p=0;p<(ssize_t)positions[t].size();++p)
            followersByDist[p] = positions[t+1].getEuclidianNeighboursBySqDist(positions[t][p], range);

        tm.push_back(followersByDist, positions[t+1].size());

        Error = (tm.getNbTraj() - nbTraj)/(double)tm.getNbTraj();
        sumError+=Error;
        if(maxError<Error) maxError=Error;

        ++show_progress;
    }
    cout<<"Trajectory creation rate : mean="<<100.0*sumError/(positions.size()-1)<<"%\tmax="<<100.0*maxError<<"%"<<endl;

    //create the trajIndex from the trajMap
    trajectories = TrajIndex(tm);
}



