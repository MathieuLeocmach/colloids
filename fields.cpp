/**
    Copyright 2010 Mathieu Leocmach

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


 * \file fields.cpp
 * \brief Implement classes to manage fields of values : scalar, vectorial, static or dynamic
 * \author Mathieu Leocmach
 * \date 3 February 2010
 *
 */

#include "fields.hpp"

using namespace std;
using namespace Colloids;

/** @brief write as vtk legacy format  */
ostream & Colloids::operator<<(std::ostream &os, const ScalarField &s)
{
	os<<"SCALARS "<< s.name<<" double\n"
			"LOOKUP_TABLE default\n";
	copy(
		s.values.begin(), s.values.end(),
		ostream_iterator<double>(os,"\n")
		);
	return os;
}

/** @brief write as vtk legacy format  */
ostream & Colloids::operator<<(std::ostream &os, const VectorField &v)
{
	os<<"VECTORS "<<v.name<<" double\n";
	for(size_t p=0;p<v.values.size();++p)
	{
		for(size_t d=0;d<3;++d)
			os<<v.values[p][d]<<" ";
		os<<endl;
	}
	return os;
}
