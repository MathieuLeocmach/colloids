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

 * \file dynamicParticles.hpp
 * \brief Defines classes for particle dynamics
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 3 March 2009
 *
 */

#ifndef dynamic_particles_H
#define dynamic_particles_H

#include "indexedParticles.hpp"
#include "traj.hpp"
#include "pv.hpp"

typedef std::pair< std::string,std::vector< std::map<size_t,double> >* >						scalarDynamicField;
typedef std::pair< std::string,std::vector< std::map<size_t, std::valarray<double> > >* >		vectorDynamicField;
typedef boost::ptr_vector< std::vector< std::set<size_t> > >									DynamicNeigbourList;


/**
    \brief Object representing a N particle system in 3D + time with 0 <= t < maxtime
*/
class DynamicParticles
{
    public:
        /**
            \brief Particles' positions, time step by time step.
            The particle given by positions[t]->at(i) is NOT (in general) the same particle as positions[t+1]->at(i)
        */
        boost::ptr_deque<IndexedParticles> positions;
        //std::deque<IndexedParticles*> positions;
        /**
            \brief The list of the particles.
            The particle given by positions[t]->at(trajectories[i][t]) IS the same particle as positions[t+1]->at(trajectories[i][t+1])
        */
        TrajIndex trajectories;

        /** \brief Time step */
        double dt;

        /** \brief radius of the particles */
        double radius;

        /** \brief Spatio-temporal index */
        boost::ptr_vector< boost::ptr_vector<RTree> > STindex;
        //std::vector< std::vector<RTree*> > STindex;

        /** constructors */
        DynamicParticles(IndexedParticles *parts,const double &time_step=1);
        DynamicParticles(const std::string &filename);
        DynamicParticles(const double &rad,const double &time_step,const std::string &base_name,const std::string &token,const size_t &t_offset, const size_t &t_size, const double &minSep=0.0);
        DynamicParticles(const double &rad,const double &time_step,const size_t &t_size);
        virtual ~DynamicParticles(){return;}

        /** access to elements */
        std::valarray<double> &operator()(const size_t &tr, const size_t &t);
        const std::valarray<double> &operator()(const size_t &tr, const size_t &t) const;
        inline size_t getNbTimeSteps() const {return positions.size();};
        BoundingBox getMaxBox() const;
        size_t getMaxSimultaneousParticles() const;

        void push_back(IndexedParticles* parts);

        /** export to various file formats */
        void save(const std::string &filename,const std::string &base_name,const std::string &token,const size_t &t_offset, const size_t &t_size) const;
        void saveAll(const std::string &filename,const std::string &base_name,const std::string &token) const;
        void exportToPV(const std::string &filename,const std::vector<std::map<size_t,unsigned char> > &labels,const size_t &stepSize=1) const;
        void exportToFLD(const std::string &postfix,const std::vector<std::map<size_t,double> > &labels,const size_t &stepSize=1,const double &threshold=0.0) const;
        void exportToVTK(
			const std::vector< scalarDynamicField > &scalars,
			const std::vector< vectorDynamicField > &vectors,
			const size_t &stepSize=1,
			const boost::ptr_vector< std::vector< std::set<size_t> > > &ngbList=boost::ptr_vector< std::vector< std::set<size_t> > >()
		) const;

        //void exportToDb(const std::string &dbname,const size_t &measurement_id) const;


        BoundingBox boundsTrajectory(const size_t &tr) const;

        //unsigned char &label(const size_t &tr, const size_t &t);
        //void setTrajLabel(const size_t &tr, const unsigned char& lab);

        //pv getPV(const size_t &t0,const size_t &nsteps,const size_t &stepSize=1);

        /** spatio-temporal indexation related **/
		void makeSTindex(const bool reindexAllFrames=true);
        template <typename Acceptor>
        std::set<size_t> getSpanning_Accepted(const size_t &t0,const size_t &t1,const Acceptor &accept) const;
        std::set<size_t> getSpanning_Enclosed(const size_t &t0,const size_t &t1,const BoundingBox &b) const;
        std::set<size_t> getEnclosed(const BoundingBox &b) const;
        std::set<size_t> getSpanning(const size_t &t0,const size_t &t1) const;
        std::set<size_t> getSpanningInside(const size_t &t0,const size_t &t1,const double &cutoff) const;
        DynamicNeigbourList getNgbList(const double &range);

        /** geometry and dynamics related **/
        virtual std::valarray<double> getDiff(const size_t &tr_from,const size_t &t_from,const size_t &tr_to,const size_t &t_to) const;
        std::valarray<double> getDrift(const std::set<size_t>&selection,const size_t &t0,const size_t &t1) const;
        void removeDrift();
        double getSD(const std::set<size_t>&selection,const size_t &t0,const size_t &t1) const;
        std::vector<double> getMSD(const std::set<size_t> &selection,const size_t &t0,const size_t &t1,const size_t &t3=0) const;
        std::vector<double> getMSD(const size_t &t0,const size_t &t1,const size_t &t3=0) const;
        std::vector<double> getISF(const std::set<size_t> &selection,const std::valarray<double> &q,const size_t &t0,const size_t &t1) const;
        std::vector<double> getISF(const std::valarray<double> &q,const size_t &t0,const size_t &t1) const;
        std::vector<double> getSelfISF(const std::set<size_t> &selection,const std::valarray<double> &q,const size_t &t0,const size_t &t1,const size_t &t3=0) const;
        std::vector<double> getSelfISF(const std::valarray<double> &q,const size_t &t0,const size_t &t1,const size_t &t3=0) const;
        std::vector<double> getSelfISF(const size_t &t0,const size_t &t1,const size_t &t3=0) const;
        void makeDynamics(std::vector<double> &MSD,std::vector<double> &ISF) const;
        void makeDynamics(std::vector<double> &MSD,std::vector<std::vector<double> >&ISF) const;
        void makeDynamics(const std::vector< std::set<size_t> >&sets,std::vector< std::vector<double> > &MSD,std::vector< std::vector<double> > &ISF) const;
        void exportDynamics(const std::string &inputPath) const;
        void exportDynamics(const std::vector< std::set<size_t> >&sets,const std::vector<std::string>&setsNames,const std::string &inputPath) const;
		vectorDynamicField averageVelocities(const std::set<size_t> &selection,const size_t &displInterval,const size_t &avgInterval) const;

		std::set<size_t> getLostNgbs(const DynamicNeigbourList &ngbList,const size_t &tr,const size_t &t_from,const size_t &t_to) const;
		scalarDynamicField getNbLostNgbs(const std::set<size_t> &selection,const DynamicNeigbourList &ngbList,const size_t &interval) const;

        AngularDistrib getMeanAngularDistribution(const std::set<size_t> &considered,const double &range) const;

        /** bond orientationa order related **/
        std::set<size_t> getBooFromFile(const std::string &prefix,std::vector<std::map< size_t,std::valarray<double> > >&qw) const;
        void makeBoo(const size_t &t, const std::set<size_t> &selection, std::map<size_t,BooData> &allBoo) const;
        void makeSBoo(const size_t &t, const std::set<size_t> &selection, const std::map<size_t,BooData> &allBoo, std::map<size_t,BooData> &SallBoo) const;
        void makeTimeAverage(
			const std::set<size_t> &selection,
			const size_t &avgInterval,
			const std::vector< std::map<size_t,double> > &timeDependant,
			std::vector< std::map<size_t,double> > &timeAveraged
		) const;
        void makeSlidingTimeAverage(
			const std::set<size_t> &selection,
			const size_t &avgInterval,
			const std::vector< std::map<size_t,double> > &timeDependant,
			std::vector< std::map<size_t,double> > &timeAveraged
		) const;
		void makeSlidingTimeAverage(
			const std::set<size_t> &selection,
			const size_t &avgInterval,
			const std::vector< std::map<size_t,std::valarray<double> > > &timeDependant,
			std::vector< std::map<size_t,std::valarray<double> > > &timeAveraged
		) const;

};

/** \brief give a reference to the position of the particle tr at time t */
inline std::valarray<double> &DynamicParticles::operator()(const size_t &tr, const size_t &t)
{
    try
    {
        return positions[t][trajectories[tr][t]];
    }
    catch(const TrajError &e)
    {
        throw IdTrajError(e,tr);
    }
}

/** \brief give a reference to the position of the particle tr at time t */
inline const std::valarray<double> &DynamicParticles::operator()(const size_t &tr, const size_t &t) const
{
    try
    {
        return positions[t][trajectories[tr][t]];
    }
    catch(const TrajError &e)
    {
        throw IdTrajError(e,tr);
    }
}

/**
    \brief get the index of the trajectories spanning from t0 to t1 and enclosed inside a given bounding box
*/
template <typename Acceptor>
std::set<size_t> DynamicParticles::getSpanning_Accepted(const size_t &t0,const size_t &t1,const Acceptor &accept) const
{
    std::set<size_t> ret;
    for(size_t start=0;start<=t0;++start)
        for(size_t stop=t1;stop<positions.size();++stop)
        {
            std::set<size_t>temp = STindex[start][stop].Query(accept, Gatherer()).gathered;
            ret.insert(temp.begin(),temp.end());
        }
    return ret;
}


#endif
