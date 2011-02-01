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
 * \date 13 February 2009
 *
 *
 */
#ifndef boo_H
#define boo_H

#include "index.hpp"

#include <valarray>
#include <complex>
#include <string>
#include <boost/array.hpp>
//#include <tvmet/Vector.h>

namespace Colloids
{
    //typedef tvmet::Vector<double, 3>            Coord;

    /** \brief Bond-Orientational-Order data
     *
     *  Coordinates qlm of the local symmetry on the pair spherical harmonics base \f$Y_{lm}(\theta,\phi)\f$
     *
     *   0 <= l <=10 (pair)
     *
     *  -l <= m <=l
     *
     *  conjugate of \f$Y_{lm}\f$ is \f$(-1)^m Y_{l(-m)}\f$ so only positive m coefficients are kept in memory
    */
    class BooData : public std::valarray< std::complex<double> >
    {
        static size_t w3j_l_offset[6],w3j_m1_offset[11];
        public:
            /** the non redundant wigner 3j coefficients for l=0,2,4,6,8,10 */
            static double w3j[91];
            static double &getW3j(const size_t &l, const int &m1, const int &m2);
            static size_t i2l[36], i2m[36];

            /** \brief default constructor */
            BooData() : std::valarray< std::complex <double> >(std::complex <double>(0.0,0.0),36){return;};
            explicit BooData(const Coord &rij);
            explicit BooData(const std::string &str);
            explicit BooData(const double* buff);

            /** \brief access to members */
            //std::complex<double> &operator()(const size_t &l, const size_t &m){return (*this)[m + l*l/4];};
            const std::complex<double> operator()(const size_t &l, const int &m) const;
            double innerProduct(const BooData &boo, const size_t &l) const;
            double normedProduct(const BooData &boo, const size_t &l) const;
            double getSumNorm(const size_t &l) const;
            std::valarray<std::complex<double> > getL(const size_t &l) const
            {return std::valarray<std::complex<double> >::operator[](std::slice(l*l/4,l+1,1));}
            bool isnull() const {return std::norm((*this)[0])+1.0==1.0;};

            double getQl(const size_t &l) const;
            std::complex<double> getWl(const size_t &l) const;
            void getInvarients(const size_t &l, double &Q, std::complex<double> &W) const;
            void getInvarients(const size_t &l, double &Q, double &w) const
            {
                std::complex<double> W(0.0,0.0);
                getInvarients(l,Q,W);
                w=W.real();
            }

			BooData rotate_by_Pi(const Coord &axis) const;
			BooData reflect(const Coord &normal) const;

            std::string toString() const;
            char* toBinary(double *output) const;
    } ;

    std::ostream& operator<< (std::ostream& out, const BooData &boo );
    std::istream& operator>> (std::istream& in, BooData &boo );

    /**	\brief Wigner D matrix (large D) to rotate spherical harmonics by Euler angles (alpha, beta, gamma) in zyz convention */
    class Wigner_D
    {
    	/**	\brief The prefactor of Wigner d matrix (small d) is independant of the Euler angles */
    	static const boost::array<double, 6*11*11> prefactor;
    	/** Tables of powers of trigonometric functions depending on Euler Angles*/
    	boost::array<std::complex<double>, 11> e_a;
		boost::array<std::complex<double>, 21> e_g;
		boost::array<double, 21> c_b, s_b;


		static const double & getPrefactor(const size_t &l, const int &m1, const int &m2)
		{
			return prefactor[l/2 + 6*abs(m1) + 66*abs(m2)];
		};
		double small_d(const int &l, const int &m2, const int &m1) const;

		public:
			Wigner_D(const double &alpha, const double &beta, const double &gamma);
			std::complex<double> operator()(const size_t &l, const size_t &m2, const int &m1) const
			{
				return e_a[m2] * small_d(l, m2, m1) * e_g[10+m1];
			};
    };



    struct cloud_exporter : public std::unary_function<const BooData&, std::string>
	{
		std::string operator()(const BooData &boo)
		{
			boost::array<double, 8> qw;
			boo.getInvarients(4, qw[0], qw[4]);
			boo.getInvarients(6, qw[1], qw[5]);
			boo.getInvarients(8, qw[2], qw[6]);
			boo.getInvarients(10, qw[3], qw[7]);
			std::ostringstream os;
			std::copy(qw.begin(), qw.end(), std::ostream_iterator<double>(os, "\t"));
			return os.str();
		}
	};
};
#endif
