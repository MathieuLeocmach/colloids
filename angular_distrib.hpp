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

 * \file boo.hpp
 * \brief Defines class for bond orientational order data
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 13 February 2009
 *
 *
 */
#ifndef angular_distrib_H
#define angular_distrib_H

#include <valarray>

/**
    \brief angular distribution between 0 and 180 degrees
    I chosse degree only because it gives an easy discreatization
*/
class AngularDistrib : public std::valarray<double>
{
    public:

        /** \brief default constructor */
        AngularDistrib() : std::valarray<double>(0.0,180){return;};
        AngularDistrib(const std::valarray<double> &va,const std::valarray<double> &vb);
} ;

#endif
