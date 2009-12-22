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

 * \file particles.hpp
 * \brief Defines classes for particles
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 26 Novemeber 2008
 *
 * Define the mother of all Particles classes
 *
 */


#ifndef particles_H
#define particles_H

#include <algorithm>
#include <deque>
#include <valarray>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
//#include <boost/operators.hpp>
#include <boost/bind.hpp>
#include <tvmet/Vector.h>

#include "index.hpp"
#include "boo_data.hpp"

namespace Colloids
{
    typedef RStarIndex_S::RTree                                                         RTree;
    typedef std::map<size_t, std::set<size_t> >                                         NgbList;
    typedef std::pair< const std::string*,std::map<size_t, double>* >					scalarField;
    typedef std::pair< const std::string*,std::map<size_t, Coord>* >	                vectorField;

    /**
        \brief defines a set of particles having the same radius
    */
    class Particles : public std::vector<Coord>
    {
        /** \brief A spatial index of the particles */
        SpatialIndex *index;

        public:
            /** \brief overall bounding box */
            BoundingBox bb;
            /** \brief (mean) radius of all the particles */
            double radius;


            /** \brief constructors and destructor */
            Particles(void) : std::vector<Coord>(0,Coord(0.0)){radius=0;return;};
            Particles(const std::vector<Coord> &data, const double &r=0.0) : std::vector<Coord>(data){radius=r;};
            Particles(const size_t &n, const double &d=0.0, const double &r=0.0);
            Particles(const std::string &filename, const double &r=0.0);
            Particles(const size_t &Nb, const BoundingBox &b, const std::string &filename, const double &r=0.0);
            virtual ~Particles(){return;}

            void push_back(const Coord &p);

            Particles cut(const double &sep);

            /** Geometric transforms    */
            Particles& operator*=(const Coord &v);
            Particles& operator*=(const double &mul);
            virtual Particles& operator+=(const Coord &v);

            /** Geometry related */
            virtual Coord getDiff(const Coord &from,const size_t &to) const;
            virtual Coord getDiff(const size_t &from,const size_t &to) const;
            virtual double getAngle(const size_t &origin,const size_t &a,const size_t &b) const;

            /** Index related   */
            static BoundingBox bounds(const Coord &center,const double &r=0.0);
            bool hasIndex() const {return !!index;};
            void setIndex(SpatialIndex *I) {index = I;}
            void makeRTreeIndex();
            virtual std::set<size_t> getEnclosed(const BoundingBox &b) const;
            BoundingBox getOverallBox() const;

            /** Spatial query and neighbours. Depends on both geometry and spatial index */
            std::set<size_t> getEuclidianNeighbours(const Coord &center, const double &range) const;
            size_t getNearestNeighbour(const Coord &center, const double &range=1.0) const;
            std::multimap<double,size_t> getEuclidianNeighboursBySqDist(const Coord &center, const double &range) const;
            double getMeanNbNeighbours(const double &range) const;
            std::deque< std::pair<size_t,size_t> > getBonds(const double &bondLength) const;
            void getNgbList(const double &bondLength, NgbList &ngbs) const;

            /**Bond Orientational Order related */
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

            /**Bond angle distribution related  */
            boost::array<double,180> getAngularDistribution(const size_t &numPt, const std::set<size_t> &ngbs) const;
            boost::array<double,180> getAngularDistribution(const size_t &numPt,const double &range) const;
            boost::array<double,180> getMeanAngularDistribution(const NgbList &selection) const;

            /**cluster */
            void growCluster(std::set<size_t> &population, std::set<size_t> &cluster, size_t center, const NgbList &ngbs);
            void segregate(std::set<size_t> &population, std::vector< std::set<size_t> > &clusters, const NgbList &ngbs);
            void segregateAll(std::vector< std::set<size_t> > &clusters, const NgbList &ngbs);
            void segregateAll(std::vector< std::set<size_t> > &clusters, const double &range);




            /** histograms*/
            struct Binner : public std::binary_function<const size_t &,const size_t &,void>
            {
                const Particles & parts;
                size_t count;
                double cutoff;

                Binner(const Particles &p, const double &nbDiameterCutOff) : parts(p)
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

                RdfBinner(const Particles &p, size_t n, const double &nbDiameterCutOff) : Binner(p,nbDiameterCutOff)
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

                G6Binner(const Particles &p, size_t n, const double &nbDiameterCutOff, const std::map<size_t,BooData> &BOO)
                : RdfBinner(p,n,nbDiameterCutOff),boo(BOO)
                {
                    g6 = std::vector<double>(n,0.0);
                };
                void operator()(const size_t &p, const size_t &q);
                void normalize(const size_t &n);
            };

            /** file outputs */
            void exportToFile(const std::string &filename) const;
            void exportToVTK(
                const std::string &filename,const std::deque< std::pair<size_t,size_t> > &bonds,
                const std::vector<scalarField> &scalars,	const std::vector<vectorField> &vectors,
                const std::string &dataName = "particles"
            ) const;
            void exportToVTK(const std::string &filename,
                const std::vector<scalarField> &scalars,	const std::vector<vectorField> &vectors,
                const std::string &dataName = "particles"
            ) const;
            void exportToVTK(const std::string &filename,
                const std::vector<scalarField> &scalars,
                const std::string &dataName = "particles"
            ) const;

            double getMinDim() const;
            virtual double getNumberDensity() const;
            double getVF() const;

            void getBooFromFile(const std::string &filename,std::map<size_t, tvmet::Vector<double, 4> >&qw) const;
            //static bool areTooClose(const std::valarray<double> &c, const Coord &d,const double &Sep);

    };

    //std::valarray<double> cross_prod(const std::valarray<double> &u,const std::valarray<double> &v);

};
#endif

