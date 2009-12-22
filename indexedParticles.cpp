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

#include "indexedParticles.hpp"
#include <boost/bind.hpp>

using namespace std;

/**
    \brief Tests if a particle can be inserted in the gap left by it's neigbours
    \param particle The particle to insert
    \param neighbours Indexes of the nearest neighbours
    \return True if the particle can be inserted
*/
bool IndexedParticles::noOverlap(const valarray<double> &p,const set<size_t> &neighbours,const double &minSep)
{
    if(minSep==0 || empty() || neighbours.empty()) return true;

    for(set<size_t>::const_iterator it=neighbours.begin();it!=neighbours.end();++it)
      if(areTooClose(p,at(*it),minSep)) return false;

    return true;
}

/** \brief Force the indexation in a R* Tree */
void IndexedParticles::makeIndex()
{
    tree = RTree();
    for(size_t i=0;i<size();++i)
        tree.Insert(i,bounds(at(i)));
}

/**
    \brief constructor from data
    \param input data
    \param radius
*/
IndexedParticles::IndexedParticles(const deque< valarray<double> > &input,const double &rad) : Particles(input)
{
    radius = rad;
    makeIndex();
    return;
}

/**
    \brief constructor from DAT file
    \param filename data
    \param rad radius
*/
IndexedParticles::IndexedParticles(const string &filename,const double &rad) : Particles(filename)
{
    radius = rad;
    makeIndex();
    return;
}

/**
    \brief constructor from DAT file
    \param filename data
    \param rad radius
*/
IndexedParticles::IndexedParticles(const size_t &Nb, const BoundingBox &b,const double &rad, const std::string &filename) : Particles(Nb,b,filename)
{
    radius = rad;
    makeIndex();
    return;
}

/**
    \brief constructor from data with minimum separation condition
    \param input data
    \param rad radius
    \param minSep minimum separation between two particles in diameter unit
    The particles closer than minSep from another particle are not taken into account.
    First in, first served.
    Automatic R* Tree indexation.
*/
IndexedParticles::IndexedParticles(const deque< valarray<double> > &input,const double &rad,const double &minSep) : Particles(0,0.0)
{
    radius = rad;
    for(size_t i=0;i<input.size();++i)
    {
        if(getEuclidianNeighbours(input[i],2.0*radius*minSep).empty())
        {
            push_back(input[i]);
            tree.Insert(size()-1,bounds(input[i]));
        }
    }
    bb=tree.getOverallBox();
    return;
}

/**
    \brief constructor from Particle with minimum separation condition
    \param input data
    \param minSep minimum separation between two particles in diameter unit
    The particles closer than minSep from another particle are not taken into account.
    First in, first served.
    Automatic R* Tree indexation.
*/
IndexedParticles::IndexedParticles(const Particles &input,const double &minSep)
{
    radius = input.radius;
    bb = input.bb;
    for(const_iterator p=begin();p!=end();++p)
    {
        if(getEuclidianNeighbours(*p,2.0*radius*minSep).empty())
        {
            push_back(*p);
            tree.Insert(size()-1,bounds(*p));
        }
    }
    return;
}

/**
    \brief constructor from DAT file with R* Tree indexing
    \param filename
    \param rad radius
    \param minSep minimum separation between two particles in diameter unit
    The particles closer than minSep from another particle are not taken into account.
    First in, first served.
    Automatic R* Tree indexation.
*/
IndexedParticles::IndexedParticles(const string &filename,const double &rad,const double &minSep) : Particles(0,0.0)
{
    radius = rad;
    size_t listSize=0, trash;
    string line;

	ifstream file(filename.c_str(), ios::in);
	if(!file)
        throw invalid_argument("No such file as "+filename);

    //Header
	file >> trash >> listSize >> trash;
    for(size_t i=0;i<3;++i)
	{
        bb.edges[i].first  = 0.0;
        file >> bb.edges[i].second;
	}

    //Data
    for(size_t i=0;i<listSize;++i)
    {
        valarray<double> p(0.0,3);
        file >> p[0] >> p[1] >> p[2];
        if(getEuclidianNeighbours(p,2.0*radius*minSep).empty())
        {
            push_back(p);
            tree.Insert(size()-1,bounds(p));
        }
    }
    file.close();
    return;
}

/** \brief translation of the box, the coordinates and the R*Tree */
IndexedParticles& IndexedParticles::operator+=(const valarray<double> &v)
{
    bb+=v;

    for(iterator p=begin();p!=end();++p)
        (*p)+=v;
	//((Particles)(*this))+=v;

	tree+=v;
	//tree.Query(RStarAcceptAny<RTree::Node,RTree::Leaf>(),TranslateBoundingBox<RTree::BoundedItem>(&v));

    return *this;
}


/**
    \brief get the index of the particle contained inside a reduction of the bounding box
    \param cutoff range to exclude from each side of the box
    \return list of the index
*/
set<size_t> IndexedParticles::getInside(const double &cutoff) const
{
    //cout<<endl;
    //get the box containing totally all the particles having their center further than the cutoff distence from each edge
    BoundingBox insideBox = bb;
    for(size_t i=0;i<3;++i)
        if(bb.edges[i].second-bb.edges[i].first>2*cutoff)
        {
            insideBox.edges[i].first  += cutoff;
            insideBox.edges[i].second -= cutoff;
            //cout << "resized dim " << i << " from " << insideBox.edges[i].first <<" to " << insideBox.edges[i].second << endl;
        }
    //cout <<" ... ";
    //gets the indexes of the particles totally contained inside this volume
    return getEnclosed(insideBox);
}
/**
    \brief get the index of the particle contained inside a reduction of the R*Tree root bounding box
    \param cutoff range to exclude from each side of the box
    \return list of the index
    The R*Tree root bounding box can be smaller than the bonding box of the "particles" structure
*/
set<size_t> IndexedParticles::getRealInside(const double &cutoff) const
{
    //cout<<endl;
    //get the box containing totally all the particles having their center further than the cutoff distence from each edge
    BoundingBox insideBox = tree.getOverallBox();
    for(size_t i=0;i<3;++i)
        if(bb.edges[i].second-bb.edges[i].first>2*cutoff)
        {
            insideBox.edges[i].first  += cutoff;
            insideBox.edges[i].second -= cutoff;
            //cout << "resized dim " << i << " from " << insideBox.edges[i].first <<" to " << insideBox.edges[i].second << endl;
        }
    //cout <<" ... ";
    //gets the indexes of the particles totally contained inside this volume
    return getEnclosed(insideBox);
}

/**
    \brief get the index of the particles enclosed inside a given bounding box
    \param b search range
    \return list of the index
*/
set<size_t> IndexedParticles::getEnclosed(const BoundingBox &b) const
{
    return tree.Query(RTree::AcceptEnclosing(b), Gatherer()).gathered;
}

/**
    \brief get the index of the particles closer than range to center (Euclidian norm)
*/
set<size_t> IndexedParticles::getEuclidianNeighbours(const valarray<double> &center, const double &range) const
{
    set<size_t> NormOneNeighbours = getEnclosed(bounds(center,range));
    set<size_t> NormTwoNeighbours;
    valarray<double> diff(0.0,3);
    double rSq = range*range;
    for(set<size_t>::const_iterator p=NormOneNeighbours.begin();p!=NormOneNeighbours.end();++p)
    {
        diff = getDiff(center,*p);
        if((diff*diff).sum()<rSq) NormTwoNeighbours.insert(NormTwoNeighbours.end(),*p);
    }
    return NormTwoNeighbours;
}

/**
    \brief get the index of the particles closer than range to center sorted by Sqare distance to the center (Euclidian norm)
*/
multimap<double,size_t> IndexedParticles::getEuclidianNeighboursBySqDist(const valarray<double> &center, const double &range) const
{
    set<size_t> NormOneNeighbours = getEnclosed(bounds(center,range));
    multimap<double,size_t> NormTwoNeighbours;
    valarray<double> diff(0.0,3);
    double rSq = range*range, distSq;
    for(set<size_t>::const_iterator p=NormOneNeighbours.begin();p!=NormOneNeighbours.end();++p)
    {
        diff = getDiff(center,*p);
        distSq = (diff*diff).sum();
        if(distSq<rSq) NormTwoNeighbours.insert(make_pair(distSq,*p));
    }
    return NormTwoNeighbours;
}

/**
    \brief get the index of the closest particle to center (Euclidian norm)
    \param range Guess of the distance to the nearest neighbour
*/
size_t IndexedParticles::getNearestNeighbour(const valarray<double> &center, const double &range) const
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
    valarray<double> diff(0.0,3);
    for(set<size_t>::const_iterator it=ngb.begin();it!=ngb.end();++it)
    {
        diff = getDiff(center,*it);
        dist = (diff*diff).sum();
        if(dist<mindist)
        {
            mindist = dist;
            nN=*it;
        }
    }
    return nN;
}

/**
    \brief get the avearge number of neigbour of the particles.
*/
double IndexedParticles::getMeanNbNeighbours(const double &range) const
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

/** \brief return the value of the spherical harmonics for the bound between two particles */
BooData IndexedParticles::sphHarm_OneBond(const size_t &center, const size_t &neighbour) const
{
    valarray<double> diff(0.0,3);
    diff = getDiff(center,neighbour);
    return BooData(diff);
}

/**
    \brief get the orientational order around a given particle
    \param numPt Index of the reference particle
    \param range maximum distance to consider a particle as a neighbour
*/
BooData IndexedParticles::getBOO(const size_t &center,const double &range) const
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
BooData IndexedParticles::getAvBOO(const map<size_t,BooData> &BOO, const size_t &center,const double &range) const
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
BooData IndexedParticles::getBOO(const size_t &center,const std::set<size_t> &ngbList) const
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
BooData IndexedParticles::getAvBOO(const std::map<size_t,BooData> &BOO, const size_t &center,const std::set<size_t> &ngbList) const
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
map<size_t,BooData> IndexedParticles::getBOOs() const
{
    map<size_t,BooData> res;
    set<size_t> inside = getRealInside(1.3*2.0*radius);
    for(set<size_t>::const_iterator p=inside.begin();p!=inside.end();++p)
        res.insert(res.end(),make_pair(*p,getBOO(*p,1.3*2.0*radius)));
    return res;
}

/**
    \brief get the averaged bond orientational order for all usefull particles
*/
map<size_t,BooData> IndexedParticles::getavBOOs(const map<size_t,BooData> &BOO) const
{
    map<size_t,BooData> res;
    set<size_t> inside = getRealInside(2.0*1.3*2.0*radius);
    for(set<size_t>::const_iterator p=inside.begin();p!=inside.end();++p)
        res.insert(res.end(),make_pair(*p,getAvBOO(BOO,*p,1.3*2.0*radius)));
    return res;
}

/** @brief export qlm in binary  */
void IndexedParticles::exportQlm(const std::map<size_t,BooData> &BOO, const std::string &outputPath) const
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
void IndexedParticles::exportQ6m(const std::map<size_t,BooData> &BOO, const std::string &outputPath) const
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
void IndexedParticles::load_q6m(const string &filename,map<size_t,BooData> &allBoo) const
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

/** @brief get all couples of particles closer than bondLength  */
deque< pair<size_t,size_t> > IndexedParticles::getBonds(const double &bondLength) const
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
void IndexedParticles::getNgbList(const double &bondLength, std::vector< std::set<size_t> > &ngbList) const
{
	ngbList.assign(size(),set<size_t>());
	transform(begin(),end(),ngbList.begin(),boost::bind(&IndexedParticles::getEuclidianNeighbours,this,_1,2.0*bondLength*radius));
}




/**
    @brief make radial and angular distribution functions for various sets of particles and export the results to file
*/
void IndexedParticles::rdf_angD(const std::vector< std::set<size_t> >&sets,const std::vector<std::string>&setsNames,const std::string &inputPath) const
{
    vector< vector<double> > g(sets.size(),vector<double>(200,0.0));
    vector< AngularDistrib > angD(sets.size());
    ostringstream head;

    const double shell =1.30;

    cout<<"Bond angle distribution for various sets of particles: ";
    for(size_t s=0;s<sets.size();++s)
    {
        cout<<setsNames[s]<<" ... ";
        head <<"\t"<<setsNames[s];
        angD[s] = getMeanAngularDistribution(sets[s],shell*2.0*radius);
    }
    cout<<endl;
    //export bond angle distribution
    saveTable(angD.begin(),angD.end(),inputPath + ".angular","probability"+head.str());
    cout<<endl;

    cout<<"Radial distribution function for various sets of particles: ";
    for(size_t s=0;s<sets.size();++s)
    {
        cout<<setsNames[s]<<" ... ";
        g[s] = getRdf(sets[s],g[s].size(),shell*2.0);
    }
    cout<<endl;
    //export g(r)
    saveTable(g.begin(),g.end(),inputPath + "_q6.rdf","r"+head.str(),2.0*shell/200.0);
    cout<<endl;
}



/** @brief Grows recursively a cluster of neighbouring particles
  * \param population The indicies of the particles than can be added to the cluster
  * \param cluster Collecting the indicies of the particles belonging to the cluster
  * \param center The particle we are presently looking at
  * \param ngbList The List of each particle's neighbours
  */
void IndexedParticles::growCluster(std::set<size_t> &population, std::set<size_t> &cluster, size_t center, const std::vector< std::set<size_t> > &ngbList)
{
    //set<size_t> ngb = getEuclidianNeighbours(at(center),range*2.0*radius);
    for(set<size_t>::const_iterator n=ngbList[center].begin();n!=ngbList[center].end();++n)
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
            growCluster(population,cluster,*n,ngbList);
        }
    }
}
/** @brief Segregate a population of particles into clusters of recursively neighbouring particles
  * \param population The indicies of the particles than can be added to the clusters
  * \param clusters The indicies of the particles belonging to each cluster
  * \param ngbList The List of each particle's neighbours
  */
void IndexedParticles::segregate(std::set<size_t> &population, std::deque< std::set<size_t> > &clusters, const std::vector< std::set<size_t> > &ngbList)
{
    size_t center = 0;
    while(!population.empty())
    {
        center = *population.begin();
        clusters.push_back(set<size_t>());
        clusters.back().insert(center);
        growCluster(population,clusters.back(),center,ngbList);
    }
}

/** @brief Segregate all particles into clusters of recursively neighbouring particles
    Prefer this version if you have already computed a neighbour list
*/
void IndexedParticles::segregateAll(std::deque< std::set<size_t> > &clusters, const std::vector< std::set<size_t> > &ngbList)
{
    set<size_t> all;
    for(size_t p=0;p<size();++p)
        all.insert(all.end(),p);
    segregate(all,clusters,ngbList);
}

/** @brief Segregate all particles into clusters of recursively neighbouring particles
    This version will compute the neighbour list from range (in diameter unit).
*/
void IndexedParticles::segregateAll(std::deque< std::set<size_t> > &clusters, const double &range)
{
    set<size_t> all;
    for(size_t p=0;p<size();++p)
        all.insert(all.end(),p);
    std::vector< std::set<size_t> > ngbList;
    getNgbList(range,ngbList);
    segregate(all,clusters,ngbList);
}

IndexedParticles::Binner::~Binner(void){};

/**	\brief Bin the particles given by selection (coupled to their neighbours). */
void IndexedParticles::Binner::operator<<(const std::set<size_t> &selection)
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
void IndexedParticles::RdfBinner::operator()(const size_t &p, const size_t &q)
{
	std::valarray<double> diff = parts.getDiff(p,q);
	g[(size_t)(sqrt((diff * diff).sum()) * scale)]++;
	count++;
};

/**	\brief Normalize the histogram. Do not bin afterward */
void IndexedParticles::RdfBinner::normalize(const size_t &n)
{
    g[0]=0.0;
    const double norm = 4.0 * M_PI * parts.getNumberDensity() / pow(scale,3.0) *n;
    for(size_t r=0;r<g.size();++r)
        g[r]/=norm;
    for(size_t r=1;r<g.size();++r)
        g[r]/=r*r;
}

/**	\brief Make and export the rdf of the selection */
std::vector<double> IndexedParticles::getRdf(const std::set<size_t> &selection, const size_t &n, const double &nbDiameterCutOff) const
{
	RdfBinner b(*this,n,nbDiameterCutOff);
	b<<selection;
	b.normalize(selection.size());
	return b.g;
}

/**	\brief Make and export the rdf */
std::vector<double> IndexedParticles::getRdf(const size_t &n, const double &nbDiameterCutOff) const
{
	return getRdf(getInside(2.0*radius*nbDiameterCutOff),n,nbDiameterCutOff);
}

/**	\brief Bin a couple of particles into the g and g6 histogram. */
void IndexedParticles::G6Binner::operator()(const size_t &p, const size_t &q)
{
	count++;
	std::valarray<double> diff = parts.getDiff(p,q);
	const size_t r = (size_t)(sqrt((diff * diff).sum()) * scale);
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
void IndexedParticles::G6Binner::normalize(const size_t &n)
{
    g6[0]=0.0;
    const double norm = 13.0/(4.0*M_PI);
    for(size_t r=1;r<g.size();++r)
        g6[r] /= norm * g[r];
	RdfBinner::normalize(n);
}

/** \brief saving radial distribution function to file */
void saveRDF(const vector<double>&g,const string &filename,const double &rscale)
{
     vector< vector<double> > temp(1,g);
     saveTable(temp.begin(),temp.end(),filename,"r\tg(r)",1/rscale);
}

