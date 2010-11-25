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

#include "dynamicClusters.hpp"
#include <ctime>
#include <boost/bind.hpp>

using namespace std;
using namespace Colloids;

/** @brief Grows recursively a cluster of neighbouring particles
  * \param population The indicies of the particles than can be added to the cluster
  * \param cluster Collecting the indicies of the particles belonging to the cluster
  * \param center The particle we are presently looking at
  * \param ngbList The List of each particle's neighbours
  */
void Colloids::growCluster(std::set<size_t> &population, std::set<size_t> &cluster, size_t center, const NgbList &ngbs)
{
    for(size_t n=0;n<ngbs[center].size();++n)
    {
        //are we able to use this particle ?
        set<size_t>::iterator toBeRemoved = population.find(ngbs[center][n]);
        if(toBeRemoved != population.end())
        {
            cluster.insert(cluster.end(),ngbs[center][n]);
            //this particle will be used now so it must be removed from the population
            // to prevent infinite recursion
            population.erase(toBeRemoved);
            //recursion
            growCluster(population,cluster,ngbs[center][n],ngbs);
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
        population.erase(population.begin());
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

/** @brief segregate at each time step a population of trajectories into clusters of recursively neighbouring particles
  * \param dynParts The DynamicParticles the clusters are made of
  * \param population The indicies of the trajectories than can be added to the clusters
  * \param ngbList The list of the neighbours each particle of each time step.
  * The cluster k of time step t0+1 is the cluster having the maximum common particles with cluster k of time t0.
  */
DynamicClusters::DynamicClusters(DynamicParticles &dynParts, std::set<size_t> &population)
{
    assign(dynParts,population);
    return;
}

/** @brief segregate at each time step a population of trajectories into clusters of recursively neighbouring particles
  * \param dynParts The DynamicParticles the clusters are made of
  * \param population The indicies of the trajectories than can be added to the clusters
  * \param ngbList The list of the neighbours each particle of each time step.
  * The cluster k of time step t0+1 is the cluster having the maximum common particles with cluster k of time t0.
  */
DynamicClusters& DynamicClusters::assign(DynamicParticles &dynParts, std::set<size_t> &population)
{
    parts = &dynParts;

    // retreive the unsorted cluster list at each time step
    vector< vector< set<size_t> > > unsorted_clusters(parts->getNbTimeSteps());
    for(size_t t=0;t<parts->getNbTimeSteps();++t)
    {
        set<size_t> popul_t;
        for(set<size_t>::const_iterator tr=population.begin();tr!=population.end();++tr)
			if(parts->trajectories[*tr].exist(t))
				popul_t.insert(popul_t.end(),parts->trajectories[*tr][t]);

        segregate(popul_t, unsorted_clusters[t], parts->positions[t].getNgbList());
    }

    //translate in terms of trajectories, removing single particle clusters
    members.resize(parts->getNbTimeSteps());
    for(size_t t=0;t<parts->getNbTimeSteps();++t)
        for(vector< set<size_t> >::const_iterator k=unsorted_clusters[t].begin();k!=unsorted_clusters[t].end();++k)
            if((*k).size()>1)
            {
                members[t].push_back(set<size_t>());
                transform(
                    (*k).begin(),(*k).end(),
                    inserter(members[t].back(),members[t].back().end()),
                    TrajIndex::Inverser(t,parts->trajectories)
                    );
            }

    //initial clusters
    TrajMap tm(members[0].size());
    double Error=0, maxError=0, sumError=0;

    //advancing time
    for(size_t t=0;t<parts->getNbTimeSteps()-1;++t)
    {
    	size_t nbTraj = tm.getNbTraj();
        //each possible follower S of a given cluster K are ranked by a distence
        // defined from the extend of their intersection with K
        vector< multimap<double,size_t> > followerByDist(members[t].size());
        for(size_t K=0; K<members[t].size(); ++K)
			for(size_t S=0;S<members[t+1].size();++S)
			{
				vector<size_t> Intersection;
				set_intersection(
					members[t][K].begin(),members[t][K].end(),
					members[t+1][S].begin(),members[t+1][S].end(),
					back_inserter(Intersection)
					);

				if(!Intersection.empty())
				{
					vector<size_t> Union;
					set_union(
						members[t][K].begin(),members[t][K].end(),
						members[t+1][S].begin(),members[t+1][S].end(),
						back_inserter(Union)
						);

					const double dist = Union.size()/(double)Intersection.size() - 1.0;
					if(dist<2)
						followerByDist[K].insert(make_pair(dist,S));
				}
            }
        tm.push_back(followerByDist, members[t+1].size());

        Error = (tm.getNbTraj() - nbTraj)/(double)tm.getNbTraj();
        sumError+=Error;
        if(maxError<Error) maxError=Error;
    }
    //cout<<"Trajectory creation rate : mean="<<100.0*sumError/(parts->getNbTimeSteps()-1)<<"%\tmax="<<100.0*maxError<<"%"<<endl;
    this->trajectories = TrajIndex(tm);
    return *this;
}

/** @brief save cluster information into a file serie  */
void DynamicClusters::save(FileSerie &serie) const
{
	for(size_t t=0; t<trajectories.nbFrames(); ++t)
	{
		ofstream output((serie%t).c_str(), ios::out | ios::trunc);
		output<<trajectories.getInverse(t).size()<<"\n";
		for(size_t K=0; K<trajectories.getInverse(t).size(); ++K)
		{
			output<<trajectories.getInverse(t)[K]<<"\t"<<members[t][K].size()<<"\t";
			transform(
				members[t][K].begin(), members[t][K].end(),
				ostream_iterator<size_t>(output, "\t"),
				TrajIndex::Converter(t, parts->trajectories)
				);
			output<<"\n";
		}
	}
}



/** @brief label the clusters on the reference DynamicParticles  */
/*ScalarDynamicField DynamicClusters::getLabels() const
{
    ScalarDynamicField labels;
    labels.first = "clusters";
    labels.second = new vector< map<size_t,double> > (parts->getNbTimeSteps());
    double label =1;
    for(TrajIndex::const_iterator K=trajectories.begin();K!=trajectories.end();++K)
    {
        for(size_t t=K->start_time;t<=K->last_time();++t)
        {
            try
            {
                transform(
                    members[t][(*K)[t]].begin(),members[t][(*K)[t]].end(),
                    inserter((*labels.second)[t],(*labels.second)[t].end()),
                    boost::bind(
                        make_pair<size_t,double>,
                        boost::bind(TrajIndex::Converter(t,parts->trajectories),_1),
                        label
                    )
                );
            }
            catch(const TrajError &e)
            {
                cerr<<"cluster error"<<endl;
                throw IdTrajError(e,K-trajectories.begin());
            }
        }
        label++;
    }
    return labels;
}*/

/** @brief bounds a cluster  */
BoundingBox DynamicClusters::bounds(const std::set<size_t> &cluster,const size_t &time)
{
    if(!parts->trajectories[*cluster.begin()].exist(time))
        throw invalid_argument("no agreement between first cluster element and time step");

    BoundingBox b = Particles::bounds((*parts)(*cluster.begin(),time),parts->radius);

    for(set<size_t>::const_iterator tr=cluster.begin();tr!=cluster.end();++tr)
    {
        if(parts->trajectories[*tr].exist(time))
            b.stretch(Particles::bounds((*parts)(*tr,time),parts->radius));
        else
            throw invalid_argument("no agreement between cluster element and time step");
    }
    return b;
}

/** @brief get Largest dx, Largest dy, Largest dz at a given time step  */
valarray<double> DynamicClusters::getLargestDelta(const size_t &time)
{
    valarray<double> max(0.0,3);
    BoundingBox b;
    double delta;
    for(deque< set<size_t> >::iterator k=members[time].begin();k!=members[time].end();++k)
    {
        b=bounds(*k,time);
        for(size_t d=0;d<3;++d)
        {
            delta = b.edges[d].second - b.edges[d].first;
            if(delta>max[d])
                max[d]=delta;
        }
    }
    return max;
}
