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

#include "particles.hpp"

using namespace std;

/**
    \brief Makes the bounding box centered on a particle
    \param r radius of the box
*/
BoundingBox Particles::bounds(const valarray<double> &center,const double &r)
{
	BoundingBox bb;

	for(size_t i=0;i<3;++i)
	{
        bb.edges[i].first  = center[i]-r;
        bb.edges[i].second = center[i]+r;
	}

	return bb;
}

/**    \brief empty list constructor */
Particles::Particles(const size_t &n, const double &d=0.0) : deque< valarray<double> >(n,valarray<double>(d,3)){return;}

/**    \brief constructor from DAT file */
Particles::Particles(const string &filename) : deque< valarray<double> >(0,valarray<double>(0.0,0))
{
    size_t listSize=0, trash;
    string line;

	ifstream file(filename.c_str(), ios::in);
    if(!file)
        throw invalid_argument("No such file as "+filename);

    //Header
	file >> trash >> listSize >> trash;
	assign(listSize,valarray<double>(0.0,3));

    for(size_t i=0;i<3;++i)
	{
        bb.edges[i].first  = 0.0;
        file >> bb.edges[i].second;
	}

    //Data
    for(size_t i=0;i<listSize;++i)
        file >> at(i)[0] >> at(i)[1] >> at(i)[2];

    file.close();
    //cout <<"done"<< endl;
    return;
}

/**    \brief constructor from GRV file */
Particles::Particles(const size_t &Nb, const BoundingBox &b, const string &filename) : deque< valarray<double> >(0,valarray<double>(0.0,0))
{
	ifstream file(filename.c_str(), ios::in);
    if(!file)
        throw invalid_argument("No such file as "+filename);

	assign(Nb,valarray<double>(0.0,3));

    bb=b;

    //Data
    for(iterator p=begin();p!=end();++p)
        file >> (*p)[0] >> (*p)[1] >> (*p)[2];

    file.close();
    return;
}

/** \brief resizing the box and rescaling the coordinates */
Particles& Particles::operator*=(const valarray<double> &v)
{
    for(size_t i=0; i<3;++i)
        bb.edges[i].second*=v[i];

    for(iterator p=begin();p!=end();++p)
        (*p)*=v;

    return *this;
}
/** \brief resizing the box, rescaling the coordinates and the radius */
Particles& Particles::operator*=(const double &mul)
{
    const valarray<double> v(mul,3);
    operator*=(v);
    radius*=mul;

    return *this;
}

/** \brief translation of the box and of the coordinates */
Particles& Particles::operator+=(const valarray<double> &v)
{
    bb+=v;

    for(iterator p=begin();p!=end();++p)
        (*p)+=v;

    return *this;
}

/** \brief get the difference vector between a position and one of the particles */
valarray<double> Particles::getDiff(const valarray<double> &from,const size_t &to) const
{
    valarray<double> diff = at(to)-from;
    return diff;
}
/** \brief get the difference vector between two particles */
valarray<double> Particles::getDiff(const size_t &from,const size_t &to) const
{
    valarray<double> diff = at(to)-at(from);
    return diff;
}
/**
    \brief get the angle between two vectors joining origin on one hand and respectively a and b on the other hand
    \return An angle in radian between 0 and Pi
*/
double Particles::getAngle(const size_t &origin,const size_t &a,const size_t &b) const
{
    valarray<double> va(0.0,3),vb(0.0,3);
    va = getDiff(origin,a);
    vb = getDiff(origin,b);
    return acos((va*vb).sum()/sqrt((va*va).sum()*(vb*vb).sum()));
}

/** \brief export the data to a dat file */
void Particles::exportToFile(const string &filename) const
{
    //cout << "export to " << filename << endl;

    ofstream output(filename.c_str(), ios::out | ios::trunc);
    if(output)
    {
      //DAT header
      output << "1\t" << size() << "\t1" << endl;
      output << bb.edges[0].second << "\t" << bb.edges[1].second << "\t" << bb.edges[2].second << endl;

      for(const_iterator p=begin();p!=end();++p)
      {
        for(size_t i=0;i<3;++i)
          output << (*p)[i] << "\t";
        output << endl;
      }
      output.close();
    }
    else
		throw invalid_argument("Cannot write on "+filename);
}

/** @brief dummy function. Use IndexedParticles for efficient implementation  */
deque< pair<size_t,size_t> > Particles::getBonds(const double &bondLength) const
{
	return deque< pair<size_t,size_t> >();
}

/** @brief Most general export to VTK Polydata format
	\param filename Name of the file to export to
	\param bonds The explicit unoriented bonds between particles
	\param scalars N Scalar fields with a name and mapping particle numbers to scalar (double) values
	\param vectors N Vector fields with a name and mapping particle numbers to vector (valarray<double>) values
	\param dataName The name of the full dataset
*/
void Particles::exportToVTK(
	const std::string &filename,
	const std::deque< std::pair<size_t,size_t> > &bonds,
	const std::vector<scalarField> &scalars,
	const std::vector<vectorField> &vectors,
	const std::string &dataName
) const
{
	ofstream output(filename.c_str(), ios::out | ios::trunc);
    if(!output)
		throw invalid_argument("Cannot write on "+filename);

	output<<"# vtk DataFile Version 3.0\n"
			<<dataName<<"\n"
			"ASCII\n"
			"DATASET POLYDATA\n"
			"POINTS "<<size()<<" double\n";
	for(const_iterator p=begin();p!=end();++p)
	{
		for(size_t d=0;d<3;++d)
			output<<(*p)[d]<<" ";
		output<<endl;
	}

	output << "LINES "<<bonds.size()<<" "<<bonds.size()*3<<endl;
	for(deque< pair<size_t,size_t> >::const_iterator b= bonds.begin();b!=bonds.end();++b)
		output<<"2 "<< (*b).first<<" "<<(*b).second<<endl;

	output<<"POINT_DATA "<<size()<<endl;
	for(size_t s=0;s<scalars.size();++s)
	{
		output<<"SCALARS "<< *scalars[s].first<<" double\n"
				"LOOKUP_TABLE default\n";

		size_t p=0;
    	for(map<size_t,double>::const_iterator l = scalars[s].second->begin();l!=scalars[s].second->end();++l)
    	{
    		while(p<(*l).first)
    		{
				output<<0<<endl;
				p++;
    		}
    		p++;
			output<<(*l).second<<endl;
    	}
    	while(p<size())
    	{
    		output<<0<<endl;
			p++;
		}
	}

	for(size_t v=0;v<vectors.size();++v)
	{
		output<<"VECTORS "<<*vectors[v].first<<" double\n";

		size_t p=0;
    	for(map<size_t, valarray<double> >::const_iterator l = vectors[v].second->begin();l!=vectors[v].second->end();++l)
    	{
    		while(p<(*l).first)
    		{
				for(size_t d=0;d<3;++d)
					output<<0<<" ";
				output<<endl;
				p++;
    		}
    		p++;
			for(size_t d=0;d<3;++d)
				output<<(*l).second[d]<<" ";
			output<<endl;
    	}
    	while(p<size())
    	{
    		for(size_t d=0;d<3;++d)
				output<<0<<" ";
			output<<endl;
			p++;
		}
	}
	output.close();
}

/** @brief exportToVTK without bonds  */
void Particles::exportToVTK(const std::string &filename, const std::vector<scalarField> &scalars, const std::vector<vectorField> &vectors, const std::string &dataName) const
{
	exportToVTK(filename,getBonds(1.3*2.0*radius),scalars,vectors,dataName);
}

/** @brief export only positions and scalar fields to VTK	*/
void Particles::exportToVTK(const std::string &filename, const std::vector<scalarField> &scalars, const std::string &dataName) const
{
	exportToVTK(filename,scalars,std::vector<vectorField>(),dataName);
}


/** \brief return the minimum dimension of the bounding box */
double Particles::getMinDim() const
{
    return min(bb.edges[0].second,min(bb.edges[1].second,bb.edges[2].second));
}

/** \brief return the number density */
double Particles::getNumberDensity() const
{
    //get the volume accessible to the particles (edge-2*radius)
    double vol =1;
    for(size_t i=0;i<3;++i)
        vol *= bb.edges[i].second-bb.edges[i].first-2*radius;
    //calculate the number density (number of particles per unit size^3)
    return size()/vol;
}

/** \brief return the volume fraction, considering a margin equal to the radius */
double Particles::getVF() const
{
    return 4*M_PI*pow(radius,3.0)/3.0 * getNumberDensity();
}

/** \brief return true if the two particles are closer together than Sep */
bool Particles::areTooClose(const valarray<double> &c, const valarray<double> &d,const double &Sep)
{
     const valarray<double> diff = c-d;
     return (diff*diff).sum()<Sep*Sep ;
}

/** @brief get q4 and q6 from a cloud file  */
void Particles::getBooFromFile(const string &filename,map< size_t,valarray<double> >&qw) const
{
	ifstream cloud(filename.c_str(), ios::in);
	if(!cloud)
		throw invalid_argument("no such file as "+filename);

	size_t p=0;
	valarray<double> v(0.0,4);
    string trash;
	//trashing the header
	cloud >> trash >> trash >> trash >> trash >> trash;
	while(!cloud.eof())
	{
		cloud >> p;
		for(size_t i=0;i<4;++i)
			cloud >> v[i];
		qw.insert(qw.end(),make_pair(p,v));
	}
	cloud.close();
}

/*deque<unsigned char>& Particles::makeLabels(unsigned char &labelingFunction(const size_t&,Particles&))
{
    labels = deque<unsigned char>(size(),0);
    for(size_t p=0;p<size();++p)
        labels[p]=labelingFunction(p,*this);

    return labels;
}
deque<unsigned char>& Particles::makeLabels(unsigned char &labelingFunction(iterator))
{
    labels = deque<unsigned char>(size(),0);
    size_t i=0;
    for(iterator p=begin();p!=end();++p)
        labels[i++]=labelingFunction(p);

    return labels;
}*/

/** \brief return the cross product between two 3D vectors (=valarray of dim 3) */
valarray<double> cross_prod(const valarray<double> &u,const valarray<double> &v)
{
    valarray<double>w(0.0,3);
    w[0]=u[1]*v[2]-u[2]*v[1];
    w[1]=u[2]*v[0]-u[0]*v[2];
    w[2]=u[0]*v[1]-u[1]*v[0];
    return w;
}

