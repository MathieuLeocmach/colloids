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

#include "indexedParticles.hpp"
#include <boost/bind.hpp>

using namespace std;
using namespace Colloids;
/**
    \brief Tests if a particle can be inserted in the gap left by it's neigbours
    \param particle The particle to insert
    \param neighbours Indexes of the nearest neighbours
    \return True if the particle can be inserted
*/
/*bool IndexedParticles::noOverlap(const valarray<double> &p,const set<size_t> &neighbours,const double &minSep)
{
    if(minSep==0 || empty() || neighbours.empty()) return true;

    for(set<size_t>::const_iterator it=neighbours.begin();it!=neighbours.end();++it)
      if(areTooClose(p,at(*it),minSep)) return false;

    return true;
}*/





/**
    @brief make radial and angular distribution functions for various sets of particles and export the results to file
*/
/*void IndexedParticles::rdf_angD(const std::vector< std::set<size_t> >&sets,const std::vector<std::string>&setsNames,const std::string &inputPath) const
{
    vector< vector<double> > g(sets.size(),vector<double>(200,0.0));
    vector< boost::array<double, 180> > angD(sets.size());
    ostringstream head;

    const double shell =1.30;

    cout<<"Bond angle distribution for various sets of particles: ";
    for(size_t s=0;s<sets.size();++s)
    {
        cout<<setsNames[s]<<" ... ";
        head <<"\t"<<setsNames[s];
        angD[s] = getMeanAngularDistribution(sets[s],shell*2.0*radius);
    }
    cout<<endl;
    //export bond angle distribution
    saveTable(angD.begin(),angD.end(),inputPath + ".angular","probability"+head.str());
    cout<<endl;

    cout<<"Radial distribution function for various sets of particles: ";
    for(size_t s=0;s<sets.size();++s)
    {
        cout<<setsNames[s]<<" ... ";
        g[s] = getRdf(sets[s],g[s].size(),shell*2.0);
    }
    cout<<endl;
    //export g(r)
    saveTable(g.begin(),g.end(),inputPath + "_q6.rdf","r"+head.str(),2.0*shell/200.0);
    cout<<endl;
}*/



/** \brief saving radial distribution function to file */
void saveRDF(const vector<double>&g,const string &filename,const double &rscale)
{
     vector< vector<double> > temp(1,g);
     saveTable(temp.begin(),temp.end(),filename,"r\tg(r)",1/rscale);
}

