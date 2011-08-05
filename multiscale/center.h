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

		explicit Center(const double &v=0.0, const double &r=0, const double &i=0) : r(r), intensity(i)
		{
			std::fill(this->begin(), this->end(), v);
		};
		explicit Center(const Center<D-1> &c, const double &additional_coord): r(c.r), intensity(c.intensity)
		{
			std::copy(c.begin(), c.end(), this->begin());
			this->back() = additional_coord;
		};
		inline double operator-(Center<D> other) const
		{
			double d = 0;
			for(size_t i=0; i<this->size(); ++i)
				d += pow((*this)[i]-other[i], 2);
			return d;
		}
	};
	typedef Center<2> Center2D;
	typedef Center<3> Center3D;

	template<int D>
	struct compare_radii : std::binary_function<bool, const Center<D>&, const Center<D>& >
	{
		bool operator()(const Center<D>& a, const Center<D> &b) const {return a.r < b.r;}
	};

}

#endif /* CENTER_H_ */
