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

#include "particles.hpp"
#include "traj.hpp"

namespace Colloids
{

    typedef std::pair< std::string,std::vector< std::map<size_t,double> >* >	scalarDynamicField;
    typedef std::pair< std::string,std::vector< std::map<size_t, Coord> >* >	vectorDynamicField;


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
            boost::ptr_vector<Particles> positions;
            //std::deque<Particles*> positions;
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
            SpatioTemporalIndex *index;

            /** constructors */
            DynamicParticles(Particles *parts,const double &time_step=1);
            DynamicParticles(const std::string &filename);
            DynamicParticles(
                const double &rad, const double &time_step,
                const std::string &base_name, const std::string &token,
                const size_t &t_offset, const size_t &t_size
                );
            DynamicParticles(const double &rad,const double &time_step,const size_t &t_size);
            virtual ~DynamicParticles(){return;}

            /** access to elements */
            Coord &operator()(const size_t &tr, const size_t &t);
            const Coord &operator()(const size_t &tr, const size_t &t) const;
            inline size_t getNbTimeSteps() const {return positions.size();};
            BoundingBox getMaxBox() const;
            size_t getMaxSimultaneousParticles() const;

            void push_back(Particles &parts);

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
            bool hasIndex() const {return !!index;};
            void setIndex(SpatioTemporalIndex *I) {index = I;}
            template <class ParentSTIndex>
            void sliceIndex(bool force = false);

            std::set<size_t> getSpanning_Enclosed(const TimeBox &b) const;
            std::set<size_t> getEnclosed(const BoundingBox &b) const;
            std::set<size_t> getSpanning(const Interval &in) const;
            std::set<size_t> getSpanningInside(const Interval &in, const double &margin) const;

            /** geometry and dynamics related **/
            virtual Coord getDiff(const size_t &tr_from,const size_t &t_from,const size_t &tr_to,const size_t &t_to) const;
            Coord getDrift(const std::set<size_t>&selection,const size_t &t0,const size_t &t1) const;
            void removeDrift();
            double getSD(const std::set<size_t>&selection,const size_t &t0,const size_t &t1) const;
            std::vector<double> getMSD(const std::set<size_t> &selection,const size_t &t0,const size_t &t1,const size_t &t3=0) const;
            std::vector<double> getMSD(const size_t &t0,const size_t &t1,const size_t &t3=0) const;
            std::vector<double> getISF(const std::set<size_t> &selection,const Coord &q,const size_t &t0,const size_t &t1) const;
            std::vector<double> getISF(const Coord &q,const size_t &t0,const size_t &t1) const;
            std::vector<double> getSelfISF(const std::set<size_t> &selection,const Coord &q,const size_t &t0,const size_t &t1,const size_t &t3=0) const;
            std::vector<double> getSelfISF(const Coord &q,const size_t &t0,const size_t &t1,const size_t &t3=0) const;
            std::vector<double> getSelfISF(const size_t &t0,const size_t &t1,const size_t &t3=0) const;
            void makeDynamics(std::vector<double> &MSD,std::vector<double> &ISF) const;
            void makeDynamics(std::vector<double> &MSD,std::vector<std::vector<double> >&ISF) const;
            void makeDynamics(const std::vector< std::set<size_t> >&sets,std::vector< std::vector<double> > &MSD,std::vector< std::vector<double> > &ISF) const;
            void exportDynamics(const std::string &inputPath) const;
            void exportDynamics(const std::vector< std::set<size_t> >&sets,const std::vector<std::string>&setsNames,const std::string &inputPath) const;
            vectorDynamicField averageVelocities(const std::set<size_t> &selection,const size_t &displInterval,const size_t &avgInterval) const;

            std::set<size_t> getLostNgbs(const size_t &tr,const size_t &t_from,const size_t &t_to) const;
            scalarDynamicField getNbLostNgbs(const std::set<size_t> &selection, const size_t &interval) const;

            //boost::array<double,180> getMeanAngularDistribution(const DynNgbList &selection) const;

            /** bond orientational order related **/
            //std::set<size_t> getBooFromFile(const std::string &prefix, std::vector<std::map<size_t, tvmet::Vector<double, 4> > >&qw) const;
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
            /*template<size_t N>
            void makeSlidingTimeAverage(
                const std::set<size_t> &selection,
                const size_t &avgInterval,
                const std::vector< std::map<size_t, tvmet::Vector<double, N> > > &timeDependant,
                std::vector< std::map<size_t, tvmet::Vector<double, N> > > &timeAveraged
            ) const;*/

    };

    /** \brief give a reference to the position of the particle tr at time t */
    inline Coord &DynamicParticles::operator()(const size_t &tr, const size_t &t)
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
    inline const Coord &DynamicParticles::operator()(const size_t &tr, const size_t &t) const
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

    /** @brief set the spatial index of each time step to a slice of the spatio-temporal index.  */
    template <class ParentSTIndex>
    void DynamicParticles::sliceIndex(bool force)
    {
        if(!this->hasIndex()) throw std::logic_error("Set a spatio-temporal index before !");
        for(size_t t=0; t<positions.size();++t)
            if(force || !positions[t].hasIndex())
                positions[t].setIndex(
                    new SpatioTemporalIndex_slice<ParentSTIndex, TrajIndex::Converter>(
                        this->index, t,
                        TrajIndex::Converter(t, this->trajectories)
                        )
                    );
    }

    /** @brief Average over time a time dependant and trajectory dependant value.
      *
      * \param selection The trajectories to treat
      * \param timeDependant The input values, function of time and of trajectories
      * \param timeAverage The output, function of the trajectories
      */
    /*template<size_t N>
    void DynamicParticles::makeSlidingTimeAverage(
        const std::set<size_t> &selection,
        const size_t &avgInterval,
        const std::vector< std::map<size_t, tvmet::Vector<double, N> > > &timeDependant,
        std::vector< std::map<size_t, tvmet::Vector<double, N> > > &timeAveraged
    ) const
    {
        typedef std::map<size_t, tvmet::Vector<double, N> >    map;
        typename map::iterator it;
        typename map::const_iterator td;
        timeAveraged.assign(timeDependant.size()-(avgInterval-1), map());
        //cout<<timeAveraged.size()<<" steps left"<<endl;
        for(size_t avt=0;avt<timeAveraged.size();++avt)
        {
            //cout<<"avt="<<avt<<" keep trajectories spanning between "<<avt<<" and "<<avt+avgInterval-1<<endl;
            for(std::set<size_t>::const_iterator tr=selection.begin();tr!=selection.end();++tr)
                if(trajectories[*tr].span(avt,avt+avgInterval-1))
                {
                    it = timeAveraged[avt].insert(timeAveraged[avt].end(),std::make_pair(*tr, tvmet::Vector<double, N>(0.0)));
                    for(size_t t=avt;t<avt+avgInterval;++t)
                    {
                        td = timeDependant[t].find(trajectories[*tr][t]);
                        if(td == timeDependant[t].end())
                        {
                            std::cerr<<"avt="<<avt<<"\tt="<<t<<"\ttr="<<*tr<<"\tstart="<<trajectories[*tr].start_time<<"\tlast="<<trajectories[*tr].last_time()<<std::endl;
                            throw std::invalid_argument("the trajectory tr has no assigned vector at time step t");
                        }
                        it->second += td->second;
                    }
                    it->second /= (double)avgInterval;
                }
        }
    }*/
};



#endif
