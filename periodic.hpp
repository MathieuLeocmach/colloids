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

#include "particles.hpp"

namespace Colloids
{
    /** \brief indexed particles with periodic boundary conditions */
    class PeriodicParticles : public Particles
    {
        public:

            PeriodicParticles(const double &rad) : Particles(0, rad){return;};
            PeriodicParticles(const std::vector<Coord> &input,const double &rad) : Particles(input,rad){return;};
            PeriodicParticles(const Particles &input) : Particles(input){return;};
            PeriodicParticles(const std::string &filename,const double &rad) : Particles(filename,rad){return;};
            PeriodicParticles(const size_t &Nb, const BoundingBox &b, const std::string &filename, const double &rad) : Particles(Nb,b,filename,rad){return;};

            inline double getPeriod(const size_t &d) const;

            void periodify(Coord &v) const;
            Coord getDiff(const Coord &from,const size_t &to) const;
            Coord getDiff(const size_t &from,const size_t &to) const;
            double getNumberDensity() const;
            std::vector<size_t> selectInside(const double &margin) const;
            std::vector<size_t> selectEnclosed(const BoundingBox &b) const;
            std::vector<size_t> selectInside_noindex(const double &margin) const{return this->selectInside(margin);};
            //vector<size_t> getEuclidianNeighbours(const valarray<double> &center, const double &range) const;
    };

    /** \brief get the periodicity according to the bounding box */
    inline double PeriodicParticles::getPeriod(const size_t &d) const
    {
        return bb.edges[d].second-bb.edges[d].first;
    }
};
#endif
