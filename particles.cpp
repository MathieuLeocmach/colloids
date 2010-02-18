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
//#include <boost/progress.hpp>

using namespace std;
using namespace Colloids;
//using namespace tvmet;



/**    \brief empty list constructor */
Particles::Particles(const size_t &n, const double &d, const double &r) : vector<Coord>(n,Coord(d,3)){radius=r;}

/**    \brief constructor from DAT file */
Particles::Particles(const string &filename, const double &r) : vector<Coord>(0,Coord(0.0,3))
{
    radius = r;
    size_t listSize=0, trash;
    string line;

	ifstream file(filename.c_str(), ios::in);
    if(!file)
        throw invalid_argument("No such file as "+filename);

    //Header
	file >> trash >> listSize >> trash;
	this->assign(listSize, Coord(0.0,3));

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
Particles::Particles(const size_t &Nb, const BoundingBox &b, const string &filename, const double &r) : vector<Coord>(0,Coord(0.0,3))
{
	radius=r;
	ifstream file(filename.c_str(), ios::in);
    if(!file)
        throw invalid_argument("No such file as "+filename);

	this->assign(Nb, Coord(0.0,3));

    bb=b;

    //Data
    for(iterator p=begin();p!=end();++p)
        file >> (*p)[0] >> (*p)[1] >> (*p)[2];

    file.close();
    return;
}

void Particles::push_back(const Coord &p)
{
    if(hasIndex())
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
    const Coord v(mul,3);
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

    if(hasIndex())
		(*index)+=v;
    return *this;
}


/**
    \brief get the angle between two vectors joining origin on one hand and respectively a and b on the other hand
    \return An angle in radian between 0 and Pi
*/
double Particles::getAngle(const size_t &origin,const size_t &a,const size_t &b) const
{
    Coord va(3),vb(3);
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

    setIndex(new RStarIndex_S(boxes));
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
    \brief get the indices of the particles closer than range to center (Euclidian norm)
*/
set<size_t> Particles::getEuclidianNeighbours(const Coord &center, const double &range) const
{
    set<size_t> NormOneNeighbours = getEnclosed(bounds(center,range));
    set<size_t> NormTwoNeighbours;
    Coord diff(3);
    double rSq = range*range;
    for(set<size_t>::const_iterator p=NormOneNeighbours.begin();p!=NormOneNeighbours.end();++p)
    {
        diff = getDiff(center,*p);
        if(dot(diff,diff)<rSq) NormTwoNeighbours.insert(NormTwoNeighbours.end(),*p);
    }
    return NormTwoNeighbours;
}

/**
    \brief get the indices of the particles closer than range to center (Euclidian norm), discarding center itself
*/
set<size_t> Particles::getEuclidianNeighbours(const size_t &center, const double &range) const
{
    set<size_t> NormOneNeighbours = getEnclosed(bounds((*this)[center],range));
    NormOneNeighbours.erase(center);
    set<size_t> NormTwoNeighbours;
    Coord diff(3);
    double rSq = range*range;
    for(set<size_t>::const_iterator p=NormOneNeighbours.begin();p!=NormOneNeighbours.end();++p)
    {
        diff = getDiff((*this)[center],*p);
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
    Coord diff(3);
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
    Coord diff(3);
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

/** @brief get the neighbours of each particle
  * \param bondLength The maximum separation between two neighbour particles. In diameter units
  A particle is not it's own neighbour.
  */
NgbList & Particles::makeNgbList(const double &bondLength)
{
    this->neighboursList.reset(new NgbList(size()));
    const double sep = 2.0*bondLength*radius;
    for(size_t p=0;p<size();++p)
        (*neighboursList)[p] = getEuclidianNeighbours(p, sep);

    return *this->neighboursList;
}

/** @brief make the neighbour list using a list of bonds  */
NgbList & Particles::makeNgbList(const BondList &bonds)
{
    this->neighboursList.reset(new NgbList(size()));
    for(BondList::const_iterator b=bonds.begin(); b!=bonds.end();++b)
    {
        (*neighboursList)[b->first].insert((*neighboursList)[b->first].end(), b->second);
        (*neighboursList)[b->second].insert((*neighboursList)[b->second].end(), b->first);
    }
    return *this->neighboursList;
}


/** \brief return the value of the spherical harmonics for the bound between two particles */
BooData Particles::sphHarm_OneBond(const size_t &center, const size_t &neighbour) const
{
    return BooData(getDiff(center,neighbour));
}

/** \brief get the orientational order around a given particle
    \param numPt Index of the reference particle
    \param ngbList List of the center's neighbours
  */
BooData Particles::getBOO(const size_t &center) const
{
	BooData boo;
	const set<size_t> & ngbList = getNgbList()[center];
    const size_t nb = ngbList.size();
    if(nb > 0)
    {
        //sum up the contribution of each neighbour to every spherical harmonic.
        for(set<size_t>::const_iterator p=ngbList.begin();p!=ngbList.end();++p)
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
BooData Particles::getCgBOO(const std::vector<BooData> &BOO, const size_t &center) const
{
    //sum up the contribution of each neighbour including the particle itself.
	BooData avBoo = BOO[center];
    const std::set<size_t> &ngbList = getNgbList()[center];
    for(set<size_t>::const_iterator p=ngbList.begin();p!=ngbList.end();++p)
            avBoo += BOO[*p];

    avBoo/=(double)(1+ngbList.size());
    return avBoo;
}

/**
    \brief get the bond orientational order for all particles
*/
void Particles::getBOOs(const set<size_t> &selection, std::vector<BooData> &BOO) const
{
    BOO.resize(size());
    for(set<size_t>::const_iterator p=selection.begin();p!=selection.end();++p)
        BOO[*p] = getBOO(*p);
}

/**
    \brief get the coarse grained bond orientational order for all particles
*/
void Particles::getCgBOOs(const set<size_t> &selection, const std::vector<BooData> &BOO, std::vector<BooData> &cgBOO) const
{
    cgBOO.resize(size());
    for(set<size_t>::const_iterator p=selection.begin();p!=selection.end();++p)
        cgBOO[*p] = getCgBOO(BOO, *p);
}

/** @brief export qlm in binary  */
void Particles::exportQlm(const std::vector<BooData> &BOO, const std::string &outputPath) const
{
    ofstream qlm;
    qlm.open(outputPath.c_str(), ios::binary | ios::trunc);
    if(!qlm.is_open())
        throw invalid_argument("No such file as "+outputPath);

    double buffer[32];
    for(vector<BooData>::const_iterator p=BOO.begin();p!=BOO.end();++p)
    {
        qlm.write(p->toBinary(&buffer[0]),32*sizeof(double));
    }
    qlm.close();
}
/** @brief export qlm for l==6 in ascii  */
void Particles::exportQ6m(const std::vector<BooData> &BOO, const std::string &outputPath) const
{
    ofstream q6m;
    q6m.open(outputPath.c_str(), std::ios::out | ios::trunc);
    if(!q6m.is_open())
        throw invalid_argument("No such file as "+outputPath);

    for(vector<BooData>::const_iterator p=BOO.begin();p!=BOO.end();++p)
    {
    	for(size_t m=0;m<=6;++m)
			q6m <<"\t"<<(*p)(6,m);
		q6m<<"\n";
    }
    q6m.close();
}

/** @brief load q6m from file as BooData  */
void Particles::load_q6m(const string &filename, vector<BooData> &BOO) const
{
    BOO.resize(size());
	ifstream f(filename.c_str(), ios::in);
	if(!f)
		throw invalid_argument("no such file as "+filename);

	size_t p=0;
	while(!f.eof())
	{
		for(size_t m=0;m<=6;++m)
			f>> BOO[p](6,m);
        p++;
	}
	f.close();
}

/** \brief Get the bond angle distribution around one particle given the list of the particles it is bounded with    */
boost::array<double,180> Particles::getAngularDistribution(const size_t &numPt) const
{
    boost::array<double,180> angD;
    const std::set<size_t> &ngbs = getNgbList()[numPt];
    fill(angD.begin(), angD.end(), 0.0);
    const size_t nb = ngbs.size();
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

/** \brief get the mean angular distribution of a given set of particles */
/*boost::array<double,180> Particles::getMeanAngularDistribution(const NgbList &selection) const
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
}*/

Particles::Binner::~Binner(void){};

/**	\brief Bin the particles given by selection (coupled to their neighbours). */
void Particles::Binner::operator<<(const std::set<size_t> &selection)
{
	//for each inside particle, select all particles having their center within the cutoff (norm infinity) and bin them
    for(std::set<size_t>::const_iterator p=selection.begin();p!=selection.end();++p)
    {
        std::set<size_t> around = parts.getEuclidianNeighbours(*p,cutoff);
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
	return getRdf(index->getInside(2.0*radius*nbDiameterCutOff), n, nbDiameterCutOff);
}

/**	\brief Bin a couple of particles into the g and g6 histogram. */
void Particles::G6Binner::operator()(const size_t &p, const size_t &q)
{
	count++;
	const size_t r = (size_t)(norm2(parts.getDiff(p, q)) * scale);
	g[r]++;

	double prod[7];
	for(size_t m=0;m<=6;++m)
		prod[m] = real(boo[p](6,m) * conj(boo[q](6,m)));
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
        output << "\n";
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
	const BondList &bonds,
	const std::vector<ScalarField> &scalars,
	const std::vector<VectorField> &vectors,
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
		output<<"\n";
	}

	output << "LINES "<<bonds.size()<<" "<<bonds.size()*3<<endl;
	for(deque< pair<size_t,size_t> >::const_iterator b= bonds.begin();b!=bonds.end();++b)
		output<<"2 "<< b->first<<" "<<b->second<<"\n";

	output<<"POINT_DATA "<<size()<<endl;
	copy(
		scalars.begin(), scalars.end(),
		ostream_iterator<ScalarField>(output)
		);
	copy(
		vectors.begin(), vectors.end(),
		ostream_iterator<VectorField>(output)
		);

	output.close();
}

/** @brief exportToVTK without bonds  */
void Particles::exportToVTK(const std::string &filename, const std::vector<ScalarField> &scalars, const std::vector<VectorField> &vectors, const std::string &dataName) const
{
	exportToVTK(filename,getBonds(),scalars,vectors,dataName);
}

/** @brief export only positions and scalar fields to VTK	*/
void Particles::exportToVTK(const std::string &filename, const std::vector<ScalarField> &scalars, const std::string &dataName) const
{
	exportToVTK(filename,scalars,std::vector<VectorField>(),dataName);
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

/** @brief get q4, q6, w4, w6 from a cloud file  */
void Particles::loadBoo(const string &filename, boost::multi_array<double,2>&qw) const
{
	ifstream cloud(filename.c_str(), ios::in);
	if(!cloud)
		throw invalid_argument("no such file as "+filename);

	string trash;
	//trashing the header
	getline(cloud, trash);

	boost::array<size_t, 2> shape = {{size(), 4}};
	qw.resize(shape);
	copy(
		istream_iterator<double>(cloud), istream_iterator<double>(),
		qw.origin()
		);

	cloud.close();
}

/** @brief from a neighbour list to a bond list  */
BondList Colloids::ngb2bonds(const NgbList& ngbList)
{
    BondList bonds;
	for(size_t p=0;p<ngbList.size();++p)
		transform(
            ngbList[p].lower_bound(p+1), ngbList[p].end(),
            back_inserter(bonds),
            boost::bind(make_pair<size_t,size_t>,p,_1)
            );
	return bonds;
}

/** @brief load bonds from file  */
BondList Colloids::loadBonds(const std::string &filename)
{
	BondList bonds;
	pair<size_t, size_t> b;
	ifstream f(filename.c_str());
	while(f.good())
	{
		f >> b.first >> b.second;
		bonds.push_back(b);
	}
	return bonds;
}



/** @brief Grows recursively a cluster of neighbouring particles
  * \param population The indicies of the particles than can be added to the cluster
  * \param cluster Collecting the indicies of the particles belonging to the cluster
  * \param center The particle we are presently looking at
  * \param ngbList The List of each particle's neighbours
  */
void Colloids::growCluster(std::set<size_t> &population, std::set<size_t> &cluster, size_t center, const NgbList &ngbs)
{
    for(set<size_t>::const_iterator n=ngbs[center].begin();n!=ngbs[center].end();++n)
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

/** @brief Segregate a population of particles into clusters of recursively neighbouring particles
  * \param population The indicies of the particles than can be added to the clusters
  * \param clusters The list of clusters with the indicies of the particles belonging to each (output)
  * \param ngbList The List of each particle's neighbours
  */
void Colloids::segregate(std::set<size_t> &population, std::vector< std::set<size_t> > &clusters, const NgbList &ngbs)
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

/** @brief Segregate all particles into clusters of recursively neighbouring particles */
void Colloids::segregateAll(std::vector< std::set<size_t> > &clusters, const Particles& parts)
{
    set<size_t> all;
    for(size_t p=0;p<parts.size();++p)
        all.insert(all.end(),p);
    segregate(all, clusters, parts.getNgbList());
}
