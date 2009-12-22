/**
 * \file pv.hpp
 * \brief Defines classes for particle motion time series
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 17 December 2008
 *
 */

#ifndef pv_H
#define pv_H

#include "particles.hpp"
#include <map>

/** \brief header of 16 bites PV file */
struct header16{
  unsigned short moltype,nmol,nbond,ndata;
  double vlx1,vlx2,vly1,vly2,vlz1,vlz2,time0,dt;
};
/** \brief header of 32 bites PV file */
struct header32{
    long moltype,nmol,nbond,ndata;
    double vlx1,vlx2,vly1,vly2,vlz1,vlz2,time0,dt;
};
/** \brief atom structure recognised by PV file */
struct atom{
  unsigned char attr;
  long x,y,z;
};

std::ostream& operator<< (std::ostream& out, const atom &ato );

std::istream& operator>> (std::istream& is, atom& ato );

/** \brief animation of tracked particles for PVwin*/
class pv : public std::deque<Particles>
{
    public:
        std::string name;
        long moltype,nmol,nbond,ndata;
        BoundingBox bb;
        double time0,dt,diam;
        std::vector< std::map<size_t,unsigned char> > labels;

        /** \brief default constructor */
        pv() : std::deque<Particles>(0,Particles())
        {
            name = "PV-32 /Shoji-Maruyama Laboratory";
            time0=0;
            dt=1;
            moltype=1;
            nmol=0;
            nbond=0;
            diam = 1.0;
            return;
        }

        /** \brief constructor from one frame */
        explicit pv(Particles ps, const std::map<size_t,unsigned char> *lab = false) : std::deque<Particles>(1,ps)
        {
            name = "PV-32 /Shoji-Maruyama Laboratory";
            bb = ps.bb;
            time0=0;
            dt=1;
            moltype=1;
            nmol=ps.size();
            nbond=0;
            diam = 2.0 * ps.radius;
            if(lab)
				labels.assign(1,*lab);
			else
				labels.assign(1,std::map<size_t,unsigned char>());
            return;
        }
        pv(const std::string &filename);

        void push_back(const Particles &ps, const std::map<size_t,unsigned char> *lab = false);
        void operator<<(const pv &a);

        void exportToPV(const std::string &filename);
};

#endif
