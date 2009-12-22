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


 * \file periodic.hpp
 * \brief Defines classes for particles spatially indexed with periodic boundary conditions
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 29 January 2009
 *
 */

#ifndef periodic_particles_H
#define periodic_particles_H

#include "indexedParticles.hpp"


/** \brief indexed particles with periodic boundary conditions */
class PeriodicParticles : public IndexedParticles
{
    public:

        PeriodicParticles(const double &rad) : IndexedParticles(rad){return;};
        PeriodicParticles(const std::deque< std::valarray<double> > &input,const double &rad) : IndexedParticles(input,rad){return;};
        PeriodicParticles(const Particles &input):IndexedParticles(input){return;};
        PeriodicParticles(const std::string &filename,const double &rad) : IndexedParticles(filename,rad){return;};
        PeriodicParticles(const std::deque< std::valarray<double> > &input,const double &rad,const double &minSep) : IndexedParticles(input,rad,minSep){return;};
        PeriodicParticles(const std::string &filename,const double &rad,const double &minSep) : IndexedParticles(filename,rad,minSep){return;};
        PeriodicParticles(const size_t &Nb, const BoundingBox &b,const double &rad, const std::string &filename) : IndexedParticles(Nb,b,rad,filename){return;};

        inline double getPeriod(const size_t &i) const;

        void periodify(std::valarray<double>&v) const;
        std::valarray<double> getDiff(const std::valarray<double> &from,const size_t &to) const;
        std::valarray<double> getDiff(const size_t &from,const size_t &to) const;
        double getNumberDensity() const;
        std::set<size_t> getInside(const double &cutoff) const;
        std::set<size_t> getRealInside(const double &cutoff) const;
        std::set<size_t> getEnclosed(const BoundingBox &b) const;
        //set<size_t> getEuclidianNeighbours(const valarray<double> &center, const double &range) const;
};

/** \brief get the periodicity according to the bounding box */
inline double PeriodicParticles::getPeriod(const size_t &i) const
{
    return bb.edges[i].second-bb.edges[i].first;
}
#endif
