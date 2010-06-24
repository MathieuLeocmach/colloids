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

 * \file boo.cpp
 * \brief Implementation of the bond orientational order related functions
 * \author Mathieu Leocmach
 * \version 0.2
 * \date 17 December 2008
 * Needs boost library
 */

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/bind.hpp>
#include "boo_data.hpp"

//double wigner3j( int l, int m1, int m2, int m3);

using namespace std;
using namespace Colloids;
//using namespace tvmet;

double BooData::w3j[91] = {
    //l=0
    1,
    //l=2
    -sqrt(2/35.),				//(0, 0,0)
    sqrt(1/70.),				//(1,-1,0)
    sqrt(2/35.),-sqrt(3/35.),	//(2,-2,0),(2,-1,-1)
    //l=4
    3*sqrt(2/1001.),							        //(0, 0,0)
    -3*sqrt(1/2002.),							        //(1,-1,0)
    -sqrt(11/182.)/3, 2*sqrt(5/1001.),				    //(2,-2,0),(2,-1,-1)
    sqrt(7/286.), -sqrt(5/143.)/3,				        //(3,-3,0),(3,-2,-1)
    sqrt(14/143.)/3, -sqrt(35/143.)/3, sqrt(5/143.),	//(4,-4,0),(4,-3,-1),(4,-2,-2)
    //l=6
    -20*sqrt(1/46189.),									                //(0, 0,0)
    10*sqrt(1/46189.),									                //(1,-1,0)
    sqrt(11/4199.), -2*sqrt(105/46189.),					            //(2,-2,0),(2,-1,-1)
    -43*sqrt(1/46189.)/2, 3*sqrt(21/92378.),			                //(3,-3,0),(3,-2,-1)
    4*sqrt(1/46189.), 2.5*sqrt(35/46189.), -6*sqrt(14/46189.),          //(4,-4,0),(4,-3,-1),(4,-2,-2)
    2.5*sqrt(11/4199.), -3*sqrt(21/4199.)/2, sqrt(7/4199.) ,	        //(5,-5,0),(5,-4,-1),(5,-3,-2)
    sqrt(11/4199.), -sqrt(77/8398.), sqrt(70/4199.), -2*sqrt(21/4199.), //(6,-6,0),(6,-5,-1),(6,-4,-2),(6,-3,-3)
    //l=8
    7*sqrt(10/96577.),                                                                          //(0, 0,0)
    -7*sqrt(5/193154.),                                                                         //(1,-1,0)
    -37*sqrt(1/965770.), 6*sqrt(14/96577.),                                                     //(2,-2,0),(2,-1,-1)
    73*sqrt(1/965770.), -3*sqrt(66/482885.),                                                    //(3,-3,0),(3,-2,-1)
    -5*sqrt(5/193154.), -8*sqrt(3/96577.), 6*sqrt(77/482885.),                                  //(4,-4,0),(4,-3,-1),(4,-2,-2)
    -sqrt(65/14858.), 3*sqrt(5/7429.), -sqrt(42/37145.),                                        //(5,-5,0),(5,-4,-1),(5,-3,-2)
    sqrt(65/14858.), 0.0, -3*sqrt(3/7429.), 2*sqrt(66/37145.),                                  //(6,-6,0),(6,-5,-1),(6,-4,-2),(6,-3,-3)
    7*sqrt(13/74290.), -sqrt(78/7429.), 3*sqrt(26/37145.), -sqrt(33/37145.),                    //(7,-7,0),(7,-6,-1),(7,-5,-2),(7,-4,-3)
    sqrt(26/37145.), -3*sqrt(13/37145.), sqrt(273/37145.), -sqrt(429/37145.), 3*sqrt(11/7429.), //(8,-8,0),(8,-7,-1),(8,-6,-2),(8,-5,-3),(8,-4,-4)
    //l=10
    -126*sqrt(7/33393355.), //( 0,  0,0)
    63*sqrt(7/33393355.),   //( 1, -1,0)
    196*sqrt(7/33393355.)/3, -7*sqrt(462/6678671.),     //( 2, -2,0),( 2,-1,-1)
    -259*sqrt(7/33393355.)/2, 7*sqrt(1001/6678671.)/3,  //( 3, -3,0),( 3,-2,-1)
    1097*sqrt(1/233753485.)/3, 59*sqrt(77/6678671.)/6, -2*sqrt(6006/6678671.),      //( 4, -4,0),( 4,-3,-1),( 4,-2,-2)
    4021*sqrt(1/233753485.)/6, -113*sqrt(55/46750697.)/2, 3*sqrt(1155/13357342.),   //( 5, -5,0),( 5,-4,-1),( 5,-3,-2)
    -914*sqrt(1/233753485.), sqrt(2926/1757545.)/3, 48*sqrt(33/46750697.), -3*sqrt(3003/6678671.),  //( 6, -6,0),( 6,-5,-1),( 6,-4,-2),( 6,-3,-3)
    -7*sqrt(119/1964315.)/3, 65*sqrt(22/2750041.)/3, -sqrt(1914/474145.), 3*sqrt(429/5500082.),     //( 7, -7,0),( 7,-6,-1),( 7,-5,-2),( 7,-4,-3)
    214*sqrt(17/13750205.)/3, -3*sqrt(561/2750041.), -2*sqrt(77/392863.)/3, 71*sqrt(143/27500410.)/3, -sqrt(2002/392863.),  //( 8, -8,0),( 8,-7,-1),( 8,-6,-2),( 8,-5,-3),( 8,-4,-4)
    3*sqrt(323/723695.), -sqrt(1309/20677.)/3, 5*sqrt(374/144739.)/3, -2*sqrt(143/144739.), sqrt(1001/206770.)/3,           //( 9, -9,0),( 9,-8,-1),( 9,-7,-2),( 9,-6,-3),( 9,-5,-4)
    2*sqrt(323/723695.)/3, -sqrt(7106/723695.)/3, 2*sqrt(561/723695.), -4*sqrt(2431/723695.)/3, 2*sqrt(2002/103385.)/3, -sqrt(1001/103385.) //(10,-10,0),(10,-9,-1),(10,-8,-2),(10,-7,-3),(10,-6,-4),(10,-5,-5)
};

size_t BooData::w3j_l_offset[6] = {0,1,5,14,30,55};
size_t BooData::w3j_m1_offset[11] = {0,1,2,4,6,9,12,16,20,25,30};

size_t BooData::i2l[36] = {0,
                           2, 2, 2,
                           4, 4, 4, 4, 4,
                           6, 6, 6, 6, 6, 6, 6,
                           8, 8, 8, 8, 8, 8, 8, 8, 8,
                           10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
size_t BooData::i2m[36] = {0,
                           0, 1, 2,
                           0, 1, 2, 3, 4,
                           0, 1, 2, 3, 4, 5, 6,
                           0, 1, 2, 3, 4, 5, 6, 7, 8,
                           0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

/** @brief get the value of the Wigner 3j symbol ((l,l,l),(m1,m2,-m1-m2))
 *  l even
 *  -l<=m1,m2,m1+m2<=l
*/
double & BooData::getW3j(const size_t &l, const int &m1, const int &m2)
{
    boost::array<int,3> m = {{abs(m1), abs(m2), abs(m1+m2)}};
    sort(m.begin(),m.end());
    //comment the tests: never saw these exeptions.
    /*if(l/2>3) throw invalid_argument("l is too big");
    if(m.back()>6) throw invalid_argument("m.back() is too big");
    if(w3j_l_offset[l/2] + w3j_m1_offset[m.back()]+m.front()>29) throw invalid_argument("total too long");*/
    return w3j[w3j_l_offset[l/2] + w3j_m1_offset[m.back()]+m.front()];
}


/** \brief constructor from one bond */
BooData::BooData(const Coord &rij): valarray< complex <double> >(36)
{
    Coord diff(3);
	double theta, phi;

	// conversion to polar coordinates (physical convention, not mathematical)
    diff = normalize(rij);
    theta = acos(diff[2]);
	if(1.0+diff[0]*diff[0]==1.0)
	{
	    if(1.0+diff[1]*diff[1]==1.0) phi=0.0;
	    else phi = M_PI / 2.0 * (diff[1] > 0.0 ? 1.0 : -1.0);
	}
	else phi = atan( diff[1] / diff[0]) + (diff[0] > 0.0 ? 0.0 :M_PI);

	//fill in with spherical harmonics
	for(int i=0; i<36; ++i)
        (*this)[i] = boost::math::spherical_harmonic(i2l[i], i2m[i], theta, phi);
    return;
}

/** \brief return the qlm values */
const complex<double> BooData::operator()(const size_t &l, const int &m) const
{
    if(m>=0)
        return (*this)[m + l*l/4];
    if((-m)%2 == 0)
        return conj((*this)[l*l/4 -m]);
    else
        return -conj((*this)[l*l/4-m]);
}

/** \brief sum over m for a given l of the norms */
double BooData::getSumNorm(const size_t &l) const
{
    double sum = 0.0;
    for(size_t m = 1; m <= l; m++)
        sum += norm((*this)(l,m));
    sum *= 2.0;
    sum += norm((*this)(l, 0u));
    return sum;
}
/** \brief Steindhardt order parameter Ql */
double BooData::getQl(const size_t &l) const
{
    return sqrt( 4.0 * M_PI * getSumNorm(l) / ( 2 * l + 1));
}
/** \brief Steindhardt order parameter Wl */
complex<double> BooData::getWl(const size_t &l) const
{
    const double sumQl = getSumNorm(l);
    complex<double> sumWl(0.0,0.0);

    //m1,m2,m3 are in [-l,l] and m1+m2+m3=0
    /*for(int m1 = -(int)l; m1 <= (int)l; m1++)
        for(int m2 = -(int)l; m2 <= (int)l; m2++)
            for(int m3 = -(int)l; m3 <= (int)l; m3++)
                if(m1+m2+m3==0)
                    sumWl += wigner3j( l, m1, m2, m3) * (*this)(l,m1) * (*this)(l,m2) * (*this)(l,m3);*/

    //new method with no repeated calculation
    for(size_t m1 = 0; m1 <= l/2; m1++)
        for(size_t m2 = m1; m2 <= l && m2<=l-m1; m2++)
            sumWl += getW3j( l, m1, m2) * (*this)(l,m1) * (*this)(l,m2) * (*this)(l,-m1-m2);
    //There are 6 possible permutations in a triplet and 2 possible signs => 12
    sumWl *= 12.0;

    if( 1.0 + sumQl != 1.0) sumWl /= pow(sumQl, 1.5);
    return sumWl;
}

/** @brief get both Ql and Wl  */
void BooData::getInvarients(const size_t &l, double &Q, std::complex<double> &W) const
{
    const double sumQl = getSumNorm(l);
    Q = sqrt( 4.0 * M_PI * sumQl / ( 2 * l + 1));

    W = complex<double>(0.0,0.0);
    int m1,m2,m3;
    for( m1 = -(int)l; m1 <= (int)l; m1++)
		for( m2 = -(int)l; m2 <= (int)l; m2++)
		{
			m3 = -(m1+m2);
			if(-(int)l<=m3 && m3<=(int)l)
				W += getW3j( l, m1, m2) * (*this)(l,m1) * (*this)(l,m2) * (*this)(l,m3);
				//W += wigner3j( l, m1, m2, m3)* (*this)(l,m1) * (*this)(l,m2) * (*this)(l,m3);
		}
    /*for(size_t m1 = 0; m1 <= l/2; m1++)
        for(size_t m2 = m1; m2 <= l && m2<=l-m1; m2++)
            W += getW3j( l, m1, m2) * (*this)(l,m1) * (*this)(l,m2) * (*this)(l,-m1-m2);
    //There are 6 possible permutations in a triplet and 2 possible signs => 12
    W *= 12.0;*/

    if(1.0 + sumQl != 1.0) W /= pow(sumQl,1.5);
}

/** @brief Export the inner data to a String
  */
string BooData::toString() const
{
    ostringstream oss;
    for(size_t i=0;i<size();++i)
        oss<<(*this)[i]<<"\t";
    return oss.str();
}

/** @brief Constructor from string
  */
 BooData::BooData(const std::string &str): valarray< complex <double> >(complex <double>(0.0,0.0),16)
{
    istringstream iss(str);
    for(size_t i=0;i<36;++i)
        iss>>(*this)[i];
    return;
}

/** @brief Export to an array of 32 doubles and yeald a pointer to 32*sizeof(double) chars
  *
  */
char * BooData::toBinary(double* output) const
{
    for(size_t i=0;i<size();++i)
    {
        *(output+2*i) = (*this)[i].real();
        *(output+2*i+1) = (*this)[i].imag();
    }
    return (char*)output;
}

/** @brief constructor form a buffer of 72 doubles
  */
BooData::BooData(const double* buff) : std::valarray< std::complex <double> >(std::complex <double>(0.0,0.0),36)
{
	for(size_t i=0;i<36;++i)
		(*this)[i] = complex<double>(*(buff+2*i),*(buff+2*i+1));
	return;
}

/** \brief output to a stream */
ostream& Colloids::operator<< (ostream& out, const BooData &boo )
{
	for(size_t i=0;i<boo.size();++i)
        out << boo[i].real() <<"\t"<< boo[i].imag() <<"\t";

    return out;
}

/** \brief input from a stream */
istream& Colloids::operator>> (istream& in, BooData &boo )
{
	double re, im;
	for(size_t i=0;i<boo.size();++i)
	{
        in >> re >> im;
        boo[i] = complex<double>(re, im);
	}
	return in;
}
