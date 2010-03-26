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
            std::auto_ptr<SpatioTemporalIndex> index;

            /** constructors */
            explicit DynamicParticles(const TrajMap &trajs, boost::ptr_vector<Particles>& positions, const double &rad=1.0,const double &time_step=1.0);
            explicit DynamicParticles(boost::ptr_vector<Particles>& positions, const double &rad=1.0,const double &time_step=1.0, const std::string &displFile="", const size_t &offset=0);
            explicit DynamicParticles(const TrajMap &trajs, FileSerie &files, const double &rad=1.0,const double &time_step=1.0);
            explicit DynamicParticles(FileSerie &files, const double &rad=1.0, const double &time_step=1.0);
            explicit DynamicParticles(const std::string &filename);

            //DynamicParticles(Particles *parts,const double &time_step=1);
            /*
            explicit DynamicParticles(
                const double &rad, const double &time_step,
                const std::string &base_name, const std::string &token,
                const size_t &t_offset, const size_t &t_size
                );
            DynamicParticles(const double &rad,const double &time_step,const size_t &t_size);*/
            virtual ~DynamicParticles(){return;}

            /** access to elements */
            Coord &operator()(const size_t &tr, const size_t &t);
            const Coord &operator()(const size_t &tr, const size_t &t) const;
            inline size_t getNbTimeSteps() const {return positions.size();};
            BoundingBox getMaxBox() const;
            size_t getMaxSimultaneousParticles() const;
            std::vector<size_t> getFrameSizes() const;

            /** export to various file formats */
            void save(const std::string &filename,const std::string &base_name,const std::string &token,const size_t &t_offset, const size_t &t_size) const;
            void exportToPV(const std::string &filename,const std::vector<std::map<size_t,unsigned char> > &labels,const size_t &stepSize=1) const;
            void exportToFLD(const std::string &postfix,const std::vector<std::map<size_t,double> > &labels,const size_t &stepSize=1,const double &threshold=0.0) const;
            void exportToVTK(
				FileSerie &files,
                std::vector< ScalarDynamicField > &scalars,
                std::vector< VectorDynamicField > &vectors
            ) const;

            //void exportToDb(const std::string &dbname,const size_t &measurement_id) const;


            BoundingBox boundsTrajectory(const size_t &tr) const;

            //unsigned char &label(const size_t &tr, const size_t &t);
            //void setTrajLabel(const size_t &tr, const unsigned char& lab);

            //pv getPV(const size_t &t0,const size_t &nsteps,const size_t &stepSize=1);

            /** spatio-temporal indexation related **/
            bool hasIndex() const {return index.get();};
            void setIndex(SpatioTemporalIndex *I) {index.reset(I);}
            template <class ParentSTIndex>
            void sliceIndex(bool force = false);

            std::vector<size_t> selectSpanning_Enclosed(const TimeBox &b) const;
            std::vector<size_t> selectEnclosed(const BoundingBox &b) const;
            std::vector<size_t> selectSpanning(const Interval &in) const;
            std::vector<size_t> selectSpanningInside(const Interval &in, const double &margin) const;

            /** geometry and dynamics related **/
            virtual Coord getDiff(const size_t &tr_from,const size_t &t_from,const size_t &tr_to,const size_t &t_to) const;
            Coord getDrift(const std::vector<size_t>&selection,const size_t &t0,const size_t &t1) const;
            Coord getDrift(const size_t &t0,const size_t &t1) const;
            void removeDrift();
            void removeDrift(const std::string &displFile, const size_t &t_offset=0);
            double getSD(const std::vector<size_t>&selection,const size_t &t0,const size_t &t1) const;
            std::vector<double> getSD(const size_t &t, const size_t &halfInterval=1) const;
            std::vector<double> getMSD(const std::vector<size_t> &selection,const size_t &t0,const size_t &t1,const size_t &t3=0) const;
            std::vector<double> getMSD(const size_t &t0,const size_t &t1,const size_t &t3=0) const;
            std::vector<double> getISF(const std::vector<size_t> &selection,const Coord &q,const size_t &t0,const size_t &t1) const;
            std::vector<double> getISF(const Coord &q,const size_t &t0,const size_t &t1) const;
            std::vector<double> getSelfISF(const std::vector<size_t> &selection,const Coord &q,const size_t &t0,const size_t &t1,const size_t &t3=0) const;
            std::vector<double> getSelfISF(const Coord &q,const size_t &t0,const size_t &t1,const size_t &t3=0) const;
            std::vector<double> getSelfISF(const size_t &t0,const size_t &t1,const size_t &t3=0) const;
            void makeDynamics(std::vector<double> &MSD,std::vector<double> &ISF) const;
            void makeDynamics(std::vector<double> &MSD,std::vector<std::vector<double> >&ISF) const;
            void makeDynamics(const std::vector< std::vector<size_t> >&sets,std::vector< std::vector<double> > &MSD,std::vector< std::vector<double> > &ISF) const;
            void exportDynamics(const std::string &inputPath) const;
            void exportDynamics(const std::vector< std::vector<size_t> >&sets,const std::vector<std::string>&setsNames,const std::string &inputPath) const;
            std::vector<Coord> velocities(const size_t &t, const size_t &halfInterval=1) const;

            std::vector<size_t> getLostNgbs(const size_t &tr,const size_t &t_from,const size_t &t_to) const;
            std::vector<double> getNbLostNgbs(const size_t &t, const size_t &halfInterval=1) const;
            TrajIndex getCages(const size_t &resolution=1) const;

            //boost::array<double,180> getMeanAngularDistribution(const DynNgbList &selection) const;

            /** bond orientational order related **/
            //std::set<size_t> getBooFromFile(const std::string &prefix, std::vector<std::map<size_t, tvmet::Vector<double, 4> > >&qw) const;
            void makeBoo(const size_t &t, const std::vector<size_t> &selection, std::map<size_t,BooData> &allBoo) const;
            void makeSBoo(const size_t &t, const std::vector<size_t> &selection, const std::map<size_t,BooData> &allBoo, std::map<size_t,BooData> &SallBoo) const;
            void makeTimeAverage(
                const std::vector<size_t> &selection,
                const size_t &avgInterval,
                const std::vector< std::map<size_t,double> > &timeDependant,
                std::vector< std::map<size_t,double> > &timeAveraged
            ) const;
            void makeSlidingTimeAverage(
                const std::vector<size_t> &selection,
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
        private:
            void fill(FileSerie &files);
            void link();

    };

    /** \brief give a reference to the position of the particle tr at time t. Complexity log(P) with P the number of particles at time t */
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

    /** \brief give a constant reference to the position of the particle tr at time t. Complexity log(P) with P the number of particles at time t */
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
