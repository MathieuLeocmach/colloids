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
#include "saveTable.hpp"

namespace Colloids
{
    /**
        \brief defines an indexed set of particles
    */
    class IndexedParticles : public Particles
    {
        public:
            /** \brief R* Tree spatial index */
            RStarIndex_S::RTree tree;

            void makeIndex();

            IndexedParticles(const double &rad) : Particles(){radius=rad;return;}
            IndexedParticles(const std::vector< std::valarray<double> > &input,const double &rad);
            IndexedParticles(const Particles &input):Particles(input){makeIndex();return;};
            IndexedParticles(const std::string &filename,const double &rad);
            IndexedParticles(const size_t &Nb, const BoundingBox &b,const double &rad, const std::string &filename);
            IndexedParticles(const std::vector< std::valarray<double> > &input,const double &rad,const double &minSep);
            IndexedParticles(const Particles &input,const double &minSep);
            IndexedParticles(const std::string &filename,const double &rad,const double &minSep);

            virtual IndexedParticles& operator+=(const std::valarray<double> &v);


            bool noOverlap(const std::valarray<double> &p,const std::set<size_t> &neighbours,const double &minSep);

            void rdf_angD(const std::vector< std::set<size_t> >&sets,const std::vector<std::string>&setsNames,const std::string &inputPath) const;

    };

    /** \brief Visitor gathering particles indexes */
    struct Gatherer {
        std::set<size_t> gathered;
        bool ContinueVisiting;

        Gatherer() : gathered(), ContinueVisiting(true) {};

        void operator()(const RStarIndex_S::RTree::Leaf * const leaf)
        {
            gathered.insert(leaf->leaf);
        }
    };

    //void addToRDF(std::vector<double>&g,const std::valarray<double>&diff,const double&scale);
    void saveRDF(const std::vector<double>&g,const std::string &filename,const double &rscale);
};
#endif

