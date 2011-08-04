/*
 * Center2D.h
 *
 *  Created on: 3 août 2011
 *      Author: mathieu
 */

#ifndef CENTER_H_
#define CENTER_H_

#include<boost/array.hpp>

namespace Colloids {

	template<int D>
	struct Center : boost::array<double, D>
	{
		double r, intensity;

		Center(const double &v=0.0, const double &r=0, const double &i=0) : r(r), intensity(i)
		{
			std::fill(this->begin(), this->end(), v);
		};
	};
	typedef Center<2> Center2D;
	typedef Center<3> Center3D;

}

#endif /* CENTER_H_ */
