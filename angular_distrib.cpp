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

 * \file angular_distrib.cpp
 * \brief Implementation of the bond angular distribution related functions
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 17 December 2008
 * Needs boost library
 */

#include "dynamicParticles.hpp"

using namespace std;

/** \brief constructor from the angle between two vectors */
AngularDistrib::AngularDistrib(const valarray<double> &va,const valarray<double> &vb) : valarray<double>(0.0,180)
{
    (*this)[(size_t)(acos((va*vb).sum()/sqrt((va*va).sum()*(vb*vb).sum()))* 180.0 / M_PI)]=1.0;
    return;
}

/** \brief returns a basic angular distribution, 1.0 for the angle between oa and ob, 0.0 elsewhere */
AngularDistrib IndexedParticles::TwoBondsAngle(const size_t &origin,const size_t &a,const size_t &b) const
{
    valarray<double> va(0.0,3),vb(0.0,3);
    va = getDiff(origin,a);
    vb = getDiff(origin,b);
    return AngularDistrib(va,vb);
}

/**
    \brief get the angular distribution of the bonds around a given particle
    \param numPt Index of the reference particle
    \param range maximum distance to consider a particle as a neighbour
*/
AngularDistrib IndexedParticles::getAngularDistribution(const size_t &numPt,const double &range) const
{
    AngularDistrib angD;
    set<size_t> EuNgb = getEuclidianNeighbours(at(numPt),range);
    const size_t nb = EuNgb.size()-1;
    if(nb > 1)
    {
        //sum up the contribution of each bond angle.
        for(set<size_t>::iterator a=EuNgb.begin();a!=EuNgb.end();++a)
            if( numPt != *a)
            {
                set<size_t>::iterator b=a;
                b++;
                while(b!=EuNgb.end())
                {
                    if( numPt != *b)
                        angD[(size_t)(getAngle(numPt,*a,*b)* 180.0 / M_PI)]=1.0;
                    b++;
                }
            }

        //scale by the number of bond angles
        angD/=(double)((nb-1)*(nb-2)/2);
    }
    return angD;
}
/**    \brief get the mean angular distribution of a given set of particles */
AngularDistrib IndexedParticles::getMeanAngularDistribution(const set<size_t> &considered,const double &range) const
{
    AngularDistrib angD;
    for(set<size_t>::iterator p=considered.begin();p!=considered.end();++p)
        angD+=getAngularDistribution(*p,range);
    angD/=(double)considered.size();
    return angD;
}

/**    \brief get the mean angular distribution of a given set of particles */
AngularDistrib DynamicParticles::getMeanAngularDistribution(const set<size_t> &considered,const double &range) const
{
    AngularDistrib angD;
    for(size_t t=0;t<positions.size();++t)
    {
        set<size_t> selection;
        for(set<size_t>::iterator tr=considered.begin();tr!=considered.end();++tr)
            selection.insert(trajectories[*tr][t]);
        AngularDistrib actualAngD = positions[t].getMeanAngularDistribution(selection,range);
        for(size_t i=0;i<angD.size();++i)
            angD[i]+=actualAngD[i];
    }
    for(size_t i=0;i<angD.size();++i)
            angD[i]/=positions.size();

    return angD;
}
