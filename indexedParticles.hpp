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


 * \file indexedParticles.hpp
 * \brief Defines class for spatially indexed particles
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 26 Novemeber 2008
 *
 * R* Tree indexation
 *
 */

#ifndef indexed_particles_H
#define indexed_particles_H

#include "particles.hpp"
#include <set>
#include <map>
#include <list>
#include "boo_data.hpp"
#include "angular_distrib.hpp"
#include "saveTable.hpp"

/**
    \brief defines an indexed set of particles
*/
class IndexedParticles : public Particles
{
    public:
        /** \brief R* Tree spatial index */
        RTree tree;

        void makeIndex();

        IndexedParticles(const double &rad) : Particles(){radius=rad;return;}
        IndexedParticles(const std::deque< std::valarray<double> > &input,const double &rad);
        IndexedParticles(const Particles &input):Particles(input){makeIndex();return;};
        IndexedParticles(const std::string &filename,const double &rad);
        IndexedParticles(const size_t &Nb, const BoundingBox &b,const double &rad, const std::string &filename);
        IndexedParticles(const std::deque< std::valarray<double> > &input,const double &rad,const double &minSep);
        IndexedParticles(const Particles &input,const double &minSep);
        IndexedParticles(const std::string &filename,const double &rad,const double &minSep);

        virtual IndexedParticles& operator+=(const std::valarray<double> &v);


        bool noOverlap(const std::valarray<double> &p,const std::set<size_t> &neighbours,const double &minSep);
        virtual std::set<size_t> getInside(const double &cutoff) const;
        virtual std::set<size_t> getRealInside(const double &cutoff) const;
        virtual std::set<size_t> getEnclosed(const BoundingBox &b) const;
        std::set<size_t> getEuclidianNeighbours(const std::valarray<double> &center, const double &range) const;
        size_t getNearestNeighbour(const std::valarray<double> &center, const double &range=1.0) const;
        std::multimap<double,size_t> getEuclidianNeighboursBySqDist(const std::valarray<double> &center, const double &range) const;
        double getMeanNbNeighbours(const double &range) const;
        virtual std::deque< std::pair<size_t,size_t> > getBonds(const double &bondLength) const;
        void getNgbList(const double &bondLength, std::vector< std::set<size_t> > &ngbList) const;

        BooData sphHarm_OneBond(const size_t &center, const size_t &neighbour) const;
        BooData getBOO(const size_t &center,const double &range) const;
        BooData getAvBOO(const std::map<size_t,BooData> &BOO, const size_t &center,const double &range) const;
        BooData getBOO(const size_t &center,const std::set<size_t> &ngbList) const;
        BooData getAvBOO(const std::map<size_t,BooData> &BOO, const size_t &center,const std::set<size_t> &ngbList) const;
        std::map<size_t,BooData> getBOOs() const;
        std::map<size_t,BooData> getavBOOs(const std::map<size_t,BooData> &BOO) const;
        void exportQlm(const std::map<size_t,BooData> &BOO, const std::string &outputPath) const;
        void exportQ6m(const std::map<size_t,BooData> &BOO, const std::string &outputPath) const;
        void load_q6m(const std::string &filename, std::map<size_t,BooData> &allBoo) const;

        AngularDistrib TwoBondsAngle(const size_t &origin,const size_t &a,const size_t &b) const;
        AngularDistrib getAngularDistribution(const size_t &numPt,const double &range) const;
        AngularDistrib getMeanAngularDistribution(const std::set<size_t> &considered,const double &range) const;

        void rdf_angD(const std::vector< std::set<size_t> >&sets,const std::vector<std::string>&setsNames,const std::string &inputPath) const;

        /**cluster */
        void segregateAll(std::deque< std::set<size_t> > &clusters, const double &range);
        void segregateAll(std::deque< std::set<size_t> > &clusters, const std::vector< std::set<size_t> > &ngbList);
        void segregate(std::set<size_t> &population, std::deque< std::set<size_t> > &clusters, const std::vector< std::set<size_t> > &ngbList);
        void growCluster(std::set<size_t> &population, std::set<size_t> &cluster, size_t center, const std::vector< std::set<size_t> > &ngbList);

        /** histograms*/

        struct Binner : public std::binary_function<const size_t &,const size_t &,void>
        {
        	const IndexedParticles & parts;
        	size_t count;
        	double cutoff;

        	Binner(const IndexedParticles &p, const double &nbDiameterCutOff) : parts(p)
        	{
        		count = 0;
        		cutoff = 2.0 * parts.radius * nbDiameterCutOff;
			};
			virtual ~Binner(void);
			virtual void operator()(const size_t &p, const size_t &q){};
			void operator<<(const std::set<size_t> &selection);
        };

        struct RdfBinner : public Binner
        {
        	std::vector<double> g;
        	double scale;

        	RdfBinner(const IndexedParticles &p, size_t n, const double &nbDiameterCutOff) : Binner(p,nbDiameterCutOff)
        	{
        		g = std::vector<double>(n,0.0);
        		scale = n / cutoff;
			};
        	void operator()(const size_t &p, const size_t &q);
			void normalize(const size_t &n);
        };

        std::vector<double> getRdf(const std::set<size_t> &selection, const size_t &n, const double &nbDiameterCutOff) const;
        std::vector<double> getRdf(const size_t &n, const double &nbDiameterCutOff) const;

        struct G6Binner : public RdfBinner
        {
        	std::vector<double> g6;
        	const std::map<size_t,BooData> &boo;

        	G6Binner(const IndexedParticles &p, size_t n, const double &nbDiameterCutOff, const std::map<size_t,BooData> &BOO)
        	: RdfBinner(p,n,nbDiameterCutOff),boo(BOO)
        	{
        		g6 = std::vector<double>(n,0.0);
			};
        	void operator()(const size_t &p, const size_t &q);
			void normalize(const size_t &n);
        };



};

/** \brief Visitor gathering particles indexes */
struct Gatherer {
	std::set<size_t> gathered;
	bool ContinueVisiting;

	Gatherer() : gathered(), ContinueVisiting(true) {};

	void operator()(const RTree::Leaf * const leaf)
	{
		gathered.insert(leaf->leaf);
	}
};

//void addToRDF(std::vector<double>&g,const std::valarray<double>&diff,const double&scale);
void saveRDF(const std::vector<double>&g,const std::string &filename,const double &rscale);

#endif

