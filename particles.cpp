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
using namespace Colloids;
using namespace tvmet;



/**    \brief empty list constructor */
Particles::Particles(const size_t &n, const double &d, const double &r) : vector<Coord>(n,Coord(d)){radius=r;}

/**    \brief constructor from DAT file */
Particles::Particles(const string &filename, const double &r) : vector<Coord>(0,Coord(0.0))
{
    radius = r;
    size_t listSize=0, trash;
    string line;

	ifstream file(filename.c_str(), ios::in);
    if(!file)
        throw invalid_argument("No such file as "+filename);

    //Header
	file >> trash >> listSize >> trash;
	this->assign(listSize,Coord(0.0));

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
Particles::Particles(const size_t &Nb, const BoundingBox &b, const string &filename, const double &r) : vector<Coord>(0,Coord(0.0))
{
	radius=r;
	ifstream file(filename.c_str(), ios::in);
    if(!file)
        throw invalid_argument("No such file as "+filename);

	this->assign(Nb, Coord(0.0));

    bb=b;

    //Data
    for(iterator p=begin();p!=end();++p)
        file >> (*p)[0] >> (*p)[1] >> (*p)[2];

    file.close();
    return;
}

void Particles::push_back(const Coord &p)
{
    if(index)
        index->insert(size(),bounds(p));
    vector<Coord>::push_back(p);
}

/** @brief return a copy with no particle closer than sep.
    First in first served
    The copy is indexed by a R*Tree
  */
Particles Particles::cut(const double &sep)
{
    Particles out;
    out.bb = this->bb;
    out.reserve(this->size());
    out.setIndex(new RStarIndex_S(vector<BoundingBox>()));
    for(iterator p = this->begin(); p!=this->end();++p)
        if(out.getEuclidianNeighbours(*p,sep).empty())
            out.push_back(*p);
    return out;
}



/** \brief resizing the box and rescaling the coordinates */
Particles& Particles::operator*=(const Coord &v)
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
    const Coord v(mul);
    operator*=(v);
    radius*=mul;

    return *this;
}

/** \brief translation of the box and of the coordinates */
Particles& Particles::operator+=(const Coord &v)
{
    bb+=v;
    for(iterator p=begin(); p!=end(); ++p)
        *p += v;

    (*index)+=v;
    return *this;
}

/** \brief get the difference vector between a position and one of the particles */
Coord Particles::getDiff(const Coord &from,const size_t &to) const
{
    Coord diff;
    diff = at(to)-from;
    return diff;
}
/** \brief get the difference vector between two particles */
Coord Particles::getDiff(const size_t &from,const size_t &to) const
{
    Coord diff;
    diff = at(to)-at(from);
    return diff;
}
/**
    \brief get the angle between two vectors joining origin on one hand and respectively a and b on the other hand
    \return An angle in radian between 0 and Pi
*/
double Particles::getAngle(const size_t &origin,const size_t &a,const size_t &b) const
{
    Coord va,vb;
    va = getDiff(origin,a);
    vb = getDiff(origin,b);
    return acos(dot(va,vb)/sqrt(dot(va,va)*dot(vb,vb)));
}

/**
    \brief Makes the bounding box centered on a particle
    \param r radius of the box
*/
BoundingBox Particles::bounds(const Coord &center,const double &r)
{
	BoundingBox bb;

	for(size_t i=0;i<3;++i)
	{
        bb.edges[i].first  = center[i]-r;
        bb.edges[i].second = center[i]+r;
	}

	return bb;
}

/** @brief make a RTree spatial index for the present particles set  */
void Particles::makeRTreeIndex()
{
    vector<BoundingBox> boxes;
    boxes.reserve(this->size());
    for(iterator p = this->begin(); p!=this->end();++p)
        boxes.push_back(bounds(*p));

    this->index = new RStarIndex_S(boxes);
}

/** @brief get the indices of the particles enclosed by a query box  */
set<size_t> Particles::getEnclosed(const BoundingBox &b) const
{
    if(!this->hasIndex()) throw logic_error("Set a spatial index before doing spatial queries !");
    return (*index)(b);
}

/** @brief getOverallBox  */
BoundingBox Particles::getOverallBox() const
{
    if(this->hasIndex())
        return index->getOverallBox();
    else
        return bb;
}



/**
    \brief get the index of the particles closer than range to center (Euclidian norm)
*/
set<size_t> Particles::getEuclidianNeighbours(const Coord &center, const double &range) const
{
    set<size_t> NormOneNeighbours = getEnclosed(bounds(center,range));
    set<size_t> NormTwoNeighbours;
    Coord diff;
    double rSq = range*range;
    for(set<size_t>::const_iterator p=NormOneNeighbours.begin();p!=NormOneNeighbours.end();++p)
    {
        diff = getDiff(center,*p);
        if(dot(diff,diff)<rSq) NormTwoNeighbours.insert(NormTwoNeighbours.end(),*p);
    }
    return NormTwoNeighbours;
}

/**
    \brief get the index of the particles closer than range to center sorted by Sqare distance to the center (Euclidian norm)
*/
multimap<double,size_t> Particles::getEuclidianNeighboursBySqDist(const Coord &center, const double &range) const
{
    set<size_t> NormOneNeighbours = getEnclosed(bounds(center,range));
    multimap<double,size_t> NormTwoNeighbours;
    Coord diff;
    double rSq = range*range, distSq;
    for(set<size_t>::const_iterator p=NormOneNeighbours.begin();p!=NormOneNeighbours.end();++p)
    {
        diff = getDiff(center,*p);
        distSq = dot(diff, diff);
        if(distSq<rSq) NormTwoNeighbours.insert(make_pair(distSq,*p));
    }
    return NormTwoNeighbours;
}

/**
    \brief get the index of the closest particle to center (Euclidian norm)
    \param range Guess of the distance to the nearest neighbour
*/
size_t Particles::getNearestNeighbour(const Coord &center, const double &range) const
{
    double rg = range;
    set<size_t> ngb = getEuclidianNeighbours(center,rg);
    //seeking for an acceptable range
    while(ngb.empty())
    {
        rg*=1.1;
        ngb = getEuclidianNeighbours(center,rg);
    }
    //if(rg!=range) cout << "you should increase the range by " << rg/range << endl;

    if (ngb.size()==1) return *(ngb.begin());

    size_t nN=size();
    double dist=0.0,mindist=rg*rg;
    Coord diff;
    for(set<size_t>::const_iterator it=ngb.begin();it!=ngb.end();++it)
    {
        diff = getDiff(center,*it);
        dist = dot(diff, diff);
        if(dist<mindist)
        {
            mindist = dist;
            nN=*it;
        }
    }
    return nN;
}

/**
    \brief get the avearge number of neigbours of the particles.
*/
double Particles::getMeanNbNeighbours(const double &range) const
{
    double res=0.0;
    set<size_t> ngb;
    for(size_t p=0;p<size();++p)
    {
        ngb=getEuclidianNeighbours(at(p),range);
        res+=(double)ngb.size();
    }
    res/=size();
    return res;
}


/** @brief get all couples of particles closer than bondLength  */
deque< pair<size_t,size_t> > Particles::getBonds(const double &bondLength) const
{
	deque< pair<size_t,size_t> > bonds;
	for(size_t p=0;p<size();++p)
	{
		set<size_t> ngb = getEuclidianNeighbours(at(p),bondLength);
		transform(ngb.lower_bound(p+1),ngb.end(),back_inserter(bonds),boost::bind(make_pair<size_t,size_t>,p,_1));
	}
	return bonds;
}

/** @brief get the neighbours of each particle
  * \param bondLength The maximum separation between two neighbour particles. In diameter units
  * \param ngbList The output
  */
void Particles::getNgbList(const double &bondLength,  NgbList &ngbs) const
{
    const double sep = 2.0*bondLength*radius;
    for(size_t p=0; p<size(); ++p)
        ngbs.insert(
            ngbs.end(),
            make_pair(p, this->getEuclidianNeighbours((*this)[p], sep))
            );
}


/** \brief return the value of the spherical harmonics for the bound between two particles */
BooData Particles::sphHarm_OneBond(const size_t &center, const size_t &neighbour) const
{
    return BooData(getDiff(center,neighbour));
}

/**
    \brief get the orientational order around a given particle
    \param numPt Index of the reference particle
    \param range maximum distance to consider a particle as a neighbour
*/
BooData Particles::getBOO(const size_t &center,const double &range) const
{
    BooData boo;
    set<size_t> EuNgb = getEuclidianNeighbours(at(center),range);
    const size_t nb = EuNgb.size()-1;
    if(nb > 0)
    {
        //sum up the contribution of each neighbour to every spherical harmonic.
        for(set<size_t>::const_iterator p=EuNgb.begin();p!=EuNgb.end();++p)
            if( *p != center)
                boo+=sphHarm_OneBond(center,*p);

        boo/=(double)nb;
    }
    return boo;
}

/**
    \brief get the averaged orientational order around a given particle
    \param BOO Array of the non-averaged orientational orders around each particle
    \param numPt Index of the reference particle
    \param range maximum distance to consider a particle as a neighbour
*/
BooData Particles::getAvBOO(const map<size_t,BooData> &BOO, const size_t &center,const double &range) const
{
    BooData avBoo;
    set<size_t> EuNgb = getEuclidianNeighbours(at(center),range);
    //sum up the contribution of each neighbour including the particle itself.
    map<size_t,BooData>::const_iterator it;
    for(set<size_t>::const_iterator p=EuNgb.begin();p!=EuNgb.end();++p)
    {
        it=BOO.find(*p);
        if(it!=BOO.end());
            avBoo += (*it).second;
    }

    avBoo/=(double)(EuNgb.size());
    return avBoo;
}

/** \brief get the orientational order around a given particle
    \param numPt Index of the reference particle
    \param ngbList List of the center's neighbours
  */
BooData Particles::getBOO(const size_t &center,const std::set<size_t> &ngbList) const
{
	BooData boo;
    const size_t nb = ngbList.size()-1;
    if(nb > 0)
    {
        //sum up the contribution of each neighbour to every spherical harmonic.
        for(set<size_t>::const_iterator p=ngbList.begin();p!=ngbList.end();++p)
            if( *p != center)
                boo+=sphHarm_OneBond(center,*p);

        boo/=(double)nb;
    }
    return boo;
}


/** \brief get the averaged orientational order around a given particle
    \param BOO Array of the non-averaged orientational orders around each particle
    \param numPt Index of the reference particle
    \param ngbList List of the center's neighbours
  */
BooData Particles::getAvBOO(const std::map<size_t,BooData> &BOO, const size_t &center,const std::set<size_t> &ngbList) const
{
	BooData avBoo;
    //sum up the contribution of each neighbour including the particle itself.
    map<size_t,BooData>::const_iterator it;
    for(set<size_t>::const_iterator p=ngbList.begin();p!=ngbList.end();++p)
    {
        it=BOO.find(*p);
        if(it!=BOO.end());
            avBoo += (*it).second;
    }

    avBoo/=(double)(ngbList.size());
    return avBoo;
}




/**
    \brief get the bond orientational order for all usefull particles
*/
map<size_t,BooData> Particles::getBOOs() const
{
    if(!this->hasIndex()) throw logic_error("Set a spatial index before looking for neighbours !");
    map<size_t,BooData> res;
    set<size_t> inside = index->getInside(1.3*2.0*radius);
    for(set<size_t>::const_iterator p=inside.begin();p!=inside.end();++p)
        res.insert(res.end(),make_pair(*p,getBOO(*p,1.3*2.0*radius)));
    return res;
}

/**
    \brief get the averaged bond orientational order for all usefull particles
*/
map<size_t,BooData> Particles::getavBOOs(const map<size_t,BooData> &BOO) const
{
    if(!this->hasIndex()) throw logic_error("Set a spatial index before looking for neighbours !");
    map<size_t,BooData> res;
    set<size_t> inside = index->getInside(2.0*1.3*2.0*radius);
    for(set<size_t>::const_iterator p=inside.begin();p!=inside.end();++p)
        res.insert(res.end(),make_pair(*p,getAvBOO(BOO,*p,1.3*2.0*radius)));
    return res;
}

/** @brief export qlm in binary  */
void Particles::exportQlm(const std::map<size_t,BooData> &BOO, const std::string &outputPath) const
{
    ofstream qlm;
    qlm.open(outputPath.c_str(), ios::binary | ios::trunc);
    if(!qlm.is_open())
        throw invalid_argument("No such file as "+outputPath);

    double buffer[32];
    for(map<size_t,BooData>::const_iterator p=BOO.begin();p!=BOO.end();++p)
    {
        qlm.write((char*)&(*p).first,sizeof(size_t));
        qlm.write((*p).second.toBinary(&buffer[0]),32*sizeof(double));
    }
    qlm.close();
}
/** @brief export qlm for l==6 in ascii  */
void Particles::exportQ6m(const std::map<size_t,BooData> &BOO, const std::string &outputPath) const
{
    ofstream q6m;
    q6m.open(outputPath.c_str(), std::ios::out | ios::trunc);
    if(!q6m.is_open())
        throw invalid_argument("No such file as "+outputPath);

    for(map<size_t,BooData>::const_iterator p=BOO.begin();p!=BOO.end();++p)
    {
    	q6m <<p->first;
    	for(size_t m=0;m<=6;++m)
			q6m <<"\t"<<p->second(6,m);
		q6m<<endl;
    }
    q6m.close();
}

/** @brief load q6m from file as BooData  */
void Particles::load_q6m(const string &filename,map<size_t,BooData> &allBoo) const
{
	ifstream f(filename.c_str(), ios::in);
	if(!f)
		throw invalid_argument("no such file as "+filename);

	size_t p=0;
	BooData boo;
	while(!f.eof())
	{
		f >> p;
		for(size_t m=0;m<=6;++m)
			f>> boo(6,m);
		allBoo.insert(allBoo.end(),make_pair(p,boo));
	}
	f.close();
}

/** \brief Get the bond angle distribution around one particle given the list of the particles it is bounded with    */
boost::array<double,180> Particles::getAngularDistribution(const size_t &numPt, const std::set<size_t> &ngbs) const
{
    boost::array<double,180> angD;
    fill(angD.begin(), angD.end(), 0.0);
    const size_t nb = ngbs.size()-1;
    if(nb > 1)
    {
        //histogram is scaled by the number of bond angles
        const double scale = nb>2 ? 1.0 / ((nb-1)*(nb-2)/2) : 1.0;
        //sum up the contribution of each bond angle.
        for(set<size_t>::const_iterator a=ngbs.begin();a!=ngbs.end();++a)
            if( numPt != *a)
            {
                set<size_t>::const_iterator b=a;
                b++;
                while(b!=ngbs.end())
                {
                    if( numPt != *b)
                        angD[(size_t)(getAngle(numPt,*a,*b)* 180.0 / M_PI)] = scale;
                    b++;
                }
            }
    }
    return angD;
}

/** \brief Get the bond angle distribution around one particle given the maximum bound length    */
boost::array<double,180> Particles::getAngularDistribution(const size_t &numPt,const double &range) const
{
    return getAngularDistribution(numPt, getEuclidianNeighbours(at(numPt), range));
}

/** \brief get the mean angular distribution of a given set of particles */
boost::array<double,180> Particles::getMeanAngularDistribution(const NgbList &selection) const
{
    boost::array<double,180> angD;
    fill(angD.begin(), angD.end(), 0.0);
    for(NgbList::const_iterator p=selection.begin();p!=selection.end();++p)
        transform(
            angD.begin(), angD.end(),
            getAngularDistribution(p->first,p->second).begin(), angD.begin(),
            plus<double>()
            );

    transform(
        angD.begin(), angD.end(),
        angD.begin(),
        bind2nd(divides<double>(), (double)selection.size())
        );
    return angD;
}

/** @brief Grows recursively a cluster of neighbouring particles
  * \param population The indicies of the particles than can be added to the cluster
  * \param cluster Collecting the indicies of the particles belonging to the cluster
  * \param center The particle we are presently looking at
  * \param ngbList The List of each particle's neighbours
  */
void Particles::growCluster(std::set<size_t> &population, std::set<size_t> &cluster, size_t center, const NgbList &ngbs)
{
    NgbList::const_iterator c = ngbs.find(center);
    if(c != ngbs.end())
    {
        for(set<size_t>::const_iterator n=c->second.begin();n!=c->second.end();++n)
        {
            //are we able to use this particle ?
            set<size_t>::iterator toBeRemoved = population.find(*n);
            if(toBeRemoved != population.end())
            {
                cluster.insert(cluster.end(),*n);
                //this particle will be used now so it must be removed from the population
                // to prevent infinite recursion
                population.erase(toBeRemoved);
                //recursion
                growCluster(population,cluster,*n,ngbs);
            }
        }
    }
}

/** @brief Segregate a population of particles into clusters of recursively neighbouring particles
  * \param population The indicies of the particles than can be added to the clusters
  * \param clusters The list of clusters with the indicies of the particles belonging to each (output)
  * \param ngbList The List of each particle's neighbours
  */
void Particles::segregate(std::set<size_t> &population, std::vector< std::set<size_t> > &clusters, const NgbList &ngbs)
{
    size_t center = 0;
    while(!population.empty())
    {
        center = *population.begin();
        clusters.push_back(set<size_t>());
        clusters.back().insert(center);
        growCluster(population,clusters.back(),center,ngbs);
    }
}

/** @brief Segregate all particles into clusters of recursively neighbouring particles
    Prefer this version if you have already computed a neighbour list
*/
void Particles::segregateAll(std::vector< std::set<size_t> > &clusters, const NgbList &ngbs)
{
    set<size_t> all;
    for(size_t p=0;p<size();++p)
        all.insert(all.end(),p);
    segregate(all,clusters,ngbs);
}

/** @brief Segregate all particles into clusters of recursively neighbouring particles
    This version will compute the neighbour list from range (in diameter unit).
*/
void Particles::segregateAll(std::vector< std::set<size_t> > &clusters, const double &range)
{
    NgbList ngbs;
    getNgbList(range,ngbs);
    segregateAll(clusters,ngbs);
}

Particles::Binner::~Binner(void){};

/**	\brief Bin the particles given by selection (coupled to their neighbours). */
void Particles::Binner::operator<<(const std::set<size_t> &selection)
{
	//for each inside particle, select all particles having their center within the cutoff (norm infinity) and bin them
    for(std::set<size_t>::const_iterator p=selection.begin();p!=selection.end();++p)
    {
        std::set<size_t> around = parts.getEuclidianNeighbours(parts[*p],cutoff);
        for(std::set<size_t>::const_iterator q=around.begin();q!=around.end();++q)
			(*this)(*p,*q);
    }
}

/**	\brief Bin a couple of particles into the histogram. */
void Particles::RdfBinner::operator()(const size_t &p, const size_t &q)
{
	g[(size_t)(norm2(parts.getDiff(p,q)) * scale)]++;
	count++;
};

/**	\brief Normalize the histogram. Do not bin afterward */
void Particles::RdfBinner::normalize(const size_t &n)
{
    g[0]=0.0;
    const double norm = 4.0 * M_PI * parts.getNumberDensity() / pow(scale,3.0) *n;
    for(size_t r=0;r<g.size();++r)
        g[r]/=norm;
    for(size_t r=1;r<g.size();++r)
        g[r]/=r*r;
}

/**	\brief Make and export the rdf of the selection */
std::vector<double> Particles::getRdf(const std::set<size_t> &selection, const size_t &n, const double &nbDiameterCutOff) const
{
	RdfBinner b(*this,n,nbDiameterCutOff);
	b<<selection;
	b.normalize(selection.size());
	return b.g;
}

/**	\brief Make and export the rdf */
std::vector<double> Particles::getRdf(const size_t &n, const double &nbDiameterCutOff) const
{
	return getRdf(index->getInside(2.0*radius*nbDiameterCutOff),n,nbDiameterCutOff);
}

/**	\brief Bin a couple of particles into the g and g6 histogram. */
void Particles::G6Binner::operator()(const size_t &p, const size_t &q)
{
	count++;
	const size_t r = (size_t)(norm2(parts.getDiff(p, q)) * scale);
	g[r]++;

	std::map<size_t,BooData>::const_iterator plm = boo.find(p), qlm = boo.find(q);
	if(plm==boo.end() ||qlm==boo.end())
		throw domain_error("Bond orientational order is not calculated for all concerned particles");
	double prod[7];
	for(size_t m=0;m<=6;++m)
		prod[m] = real(plm->second(6,m) * conj(qlm->second(6,m)));
	g6[r] += prod[0];
	for(size_t m=1; m<=6;++m)
		g6[r] += 2*prod[m];
};

/**	\brief Normalize the histogram. Do not bin afterward */
void Particles::G6Binner::normalize(const size_t &n)
{
    g6[0]=0.0;
    const double norm = 13.0/(4.0*M_PI);
    for(size_t r=1;r<g.size();++r)
        g6[r] /= norm * g[r];
	RdfBinner::normalize(n);
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

/** @brief Most general export to VTK Polydata format
	\param filename Name of the file to export to
	\param bonds The explicit unoriented bonds between particles
	\param scalars N Scalar fields with a name and mapping particle numbers to scalar (double) values
	\param vectors N Vector fields with a name and mapping particle numbers to vector (Coord) values
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
    	for(map<size_t, Coord>::const_iterator l = vectors[v].second->begin();l!=vectors[v].second->end();++l)
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
/*bool Particles::areTooClose(const Coord &c, const Coord &d,const double &Sep)
{
     const valarray<double> diff = c-d;
     return dot(diff*diff).sum()<Sep*Sep ;
}*/

/** @brief get q4 and q6 from a cloud file  */
void Particles::getBooFromFile(const string &filename, map<size_t, tvmet::Vector<double, 4> >&qw) const
{
	ifstream cloud(filename.c_str(), ios::in);
	if(!cloud)
		throw invalid_argument("no such file as "+filename);

	size_t p=0;
	tvmet::Vector<double, 4> v;
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

/** \brief return the cross product between two 3D vectors (=valarray of dim 3) */
/*valarray<double> cross_prod(const valarray<double> &u,const valarray<double> &v)
{
    valarray<double>w(0.0,3);
    w[0]=u[1]*v[2]-u[2]*v[1];
    w[1]=u[2]*v[0]-u[0]*v[2];
    w[2]=u[0]*v[1]-u[1]*v[0];
    return w;
}*/
