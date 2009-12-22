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

 * \file particles.hpp
 * \brief Defines classes for particles
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 26 Novemeber 2008
 *
 * Define the mother of all Particles classes
 *
 */


#ifndef particles_H
#define particles_H

#include "RStarTree/RStarTree.h"
#include <algorithm>
#include <deque>
#include <valarray>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <boost/ptr_container/ptr_container.hpp>
#include <boost/operators.hpp>

typedef RStarTree<size_t, 3, 2, 16,double> 											RTree;

typedef RTree::BoundingBox			                								BoundingBox;

typedef std::pair< const std::string*,std::map<size_t,double>* >					scalarField;
typedef std::pair< const std::string*,std::map<size_t, std::valarray<double> >* >	vectorField;

/**
    \brief defines a set of particles
    Reading and writing from files and data
*/
class Particles : public std::deque< std::valarray<double> >
{
    public:
        /** \brief overall bounding box */
        BoundingBox bb;
        /** \brief (mean) radius of all the particles */
        double radius;
        //std::deque<unsigned char> labels;

        /** \brief constructor from data */
        Particles(void) : std::deque< std::valarray<double> >(0,std::valarray<double>(0.0,3)){return;};
        Particles(const std::deque< std::valarray<double> > &data) : std::deque< std::valarray<double> >(data){return;};
        Particles(const size_t &n, const double &d);
        Particles(const std::string &filename);
        Particles(const size_t &Nb, const BoundingBox &b, const std::string &filename);
        virtual ~Particles(){return;}

        Particles& operator*=(const std::valarray<double> &v);
        Particles& operator*=(const double &mul);
        virtual Particles& operator+=(const std::valarray<double> &v);

        static BoundingBox bounds(const std::valarray<double> &center,const double &r=0.0);

        virtual std::valarray<double> getDiff(const std::valarray<double> &from,const size_t &to) const;
        virtual std::valarray<double> getDiff(const size_t &from,const size_t &to) const;

        double getAngle(const size_t &origin,const size_t &a,const size_t &b) const;
        virtual std::deque< std::pair<size_t,size_t> > getBonds(const double &bondLength) const;

        void exportToFile(const std::string &filename) const;
        void exportToVTK(
			const std::string &filename,const std::deque< std::pair<size_t,size_t> > &bonds,
			const std::vector<scalarField> &scalars,	const std::vector<vectorField> &vectors,
			const std::string &dataName = "particles"
		) const;
        void exportToVTK(const std::string &filename,
			const std::vector<scalarField> &scalars,	const std::vector<vectorField> &vectors,
			const std::string &dataName = "particles"
		) const;
		void exportToVTK(const std::string &filename,
			const std::vector<scalarField> &scalars,
			const std::string &dataName = "particles"
		) const;

        //void exportToDb(const string &filename,const size_t &measurement_id,const size_t &time_step) const;
        double getMinDim() const;
        virtual double getNumberDensity() const;
        double getVF() const;

        void getBooFromFile(const std::string &filename,std::map< size_t,std::valarray<double> >&qw) const;

        //std::deque<unsigned char>& makeLabels(unsigned char &labelingFunction(const size_t&,Particles&));
        //std::deque<unsigned char>& makeLabels(unsigned char &labelingFunction(iterator));

        static bool areTooClose(const std::valarray<double> &c, const std::valarray<double> &d,const double &Sep);

};

std::valarray<double> cross_prod(const std::valarray<double> &u,const std::valarray<double> &v);
#endif

