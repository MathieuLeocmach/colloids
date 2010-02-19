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

 * \file dynamicClusters.hpp
 * \brief Defines classes for time tracked clusters
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 23 April 2009
 *
 */

#ifndef dynamic_clusters_H
#define dynamic_clusters_H

#include "dynamicParticles.hpp"
namespace Colloids
{

    void growCluster(std::set<size_t> &population, std::set<size_t> &cluster, size_t center, const NgbList &ngbs);
    void segregate(std::set<size_t> &population, std::vector< std::set<size_t> > &clusters, const NgbList &ngbs);
    void segregateAll(std::vector< std::set<size_t> > &clusters, const Particles &parts);

    /**
        \brief Object representing clusters evolving in time
    */
    class DynamicClusters
    {

        public:
            /** \brief The DynamicParticles of which the clusters are made of */
            DynamicParticles *parts;

            /**
                \brief Clusters' members, time step by time step.
                The cluster given by members[t][i] is NOT (in general) the same cluster as members[t+1][i]
            */
            std::deque< std::deque<std::set<size_t> > > members;
            /**
                \brief The list of the clusters.
                The cluster given by members[t][trajectories[i][t]] IS the same cluster as members[t+1][trajectories[i][t+1]]
            */
            TrajIndex trajectories;

            DynamicClusters(DynamicParticles &dynParts, std::set<size_t> &population);

            DynamicClusters& assign(DynamicParticles &dynParts, std::set<size_t> &population);

            void save(FileSerie &serie) const;

            ScalarDynamicField getLabels() const;

            BoundingBox bounds(const std::set<size_t> &cluster,const size_t &time);

            std::valarray<double> getLargestDelta(const size_t &time);
    };
};
#endif
