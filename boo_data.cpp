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

#include "boo_data.hpp"
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/bind.hpp>

//double wigner3j( int l, int m1, int m2, int m3);

using namespace std;

double BooData::w3j[30] = {
    //l=0
    1,
    //l=2
    -0.239045722,				//(0, 0,0)
    0.119522861,				//(1,-1,0)
    0.239045722,-0.292770022,	//(2,-2,0),(2,-1,-1)
    //l=4
    0.134097047,							//(0, 0,0)
    -0.067048523,							//(1,-1,0)
    -0.081948195,0.141350699,				//(2,-2,0),(2,-1,-1)
    0.156446555,-0.062329799,				//(3,-3,0),(3,-2,-1)
    0.104297703,-0.164909148,0.186989398,	//(4,-4,0),(4,-3,-1),(4,-2,-2)
    //l=6
    -0.0930595,											//(0, 0,0)
    0.04652975,											//(1,-1,0)
    0.051182725,-0.095357612,							//(2,-2,0),(2,-1,-1)
    -0.100038963,0.045232087,							//(3,-3,0),(3,-2,-1)
    0.0186119,0.068818428,-0.10445903,					//(4,-4,0),(4,-3,-1),(4,-2,-2)
    0.127956813,-0.106078646,0.04082969,				//(5,-5,0),(5,-4,-1),(5,-3,-2)
    0.051182725,-0.095754111,0.129114817,-0.141438195	//(6,-6,0),(6,-5,-1),(6,-4,-2),(6,-3,-3)
};

size_t BooData::w3j_l_offset[4] = {0,1,5,14};
size_t BooData::w3j_m1_offset[7] = {0,1,2,4,6,9,12};

/** @brief get the value of the Wigner 3j symbol ((l,l,l),(m1,m2,-m1-m2))
 *  l even
 *  -l<=m1,m2,m1+m2<=l
*/
double & BooData::getW3j(const size_t &l, const int &m1, const int &m2)
{
    vector<int> m;
    m.push_back(abs(m1));
    m.push_back(abs(m2));
    m.push_back(abs(m1+m2));
    sort(m.begin(),m.end());
    if(l/2>3) throw invalid_argument("l is too big");
    if(m.back()>6) throw invalid_argument("m.back() is too big");
    if(w3j_l_offset[l/2] + w3j_m1_offset[m.back()]+m.front()>29) throw invalid_argument("total too long");
    return w3j[w3j_l_offset[l/2] + w3j_m1_offset[m.back()]+m.front()];

    /*sort(m.begin(),m.end(),boost::bind(less<int>,boost::bind(abs<int>,_1),boost::bind(abs<int>,_2)));
    return BooData::w3j[w3j_l[l/2] + m1*(l-m1+2)];*/
}


/** \brief constructor from one bond */
BooData::BooData(const valarray<double> &rij): valarray< complex <double> >(complex <double>(0.0,0.0),16)
{
    valarray<double> diff(0.0,3);
	double theta, phi;

	// conversion to polar coordinates (physical convention, not mathematical)
    diff=rij/sqrt((rij*rij).sum());
    theta = acos(diff[2]);
	if(diff[0]==0.0)
	{
	    if(diff[1]==0.0) phi=0.0;
	    else phi = M_PI / 2.0 * (diff[1] > 0.0 ? 1.0 : -1.0);
	}
	else phi = atan( diff[1] / diff[0]) + (diff[0] > 0.0 ? 0.0 :M_PI);

	//fill in with spherical harmonics
	for(size_t l = 0; l <= 6; l+=2)
        for(size_t m = 0; m <= l; m++)
            (*this)(l,m) = boost::math::spherical_harmonic(l, m, theta, phi);
    return;
}

/** \brief access to members */
complex<double> & BooData::operator()(const size_t &l, const size_t &m)
{
    return (*this)[m + l*l/4];
}

/** \brief return the qlm values */
const complex<double> BooData::operator()(const size_t &l, const int &m) const
{
    if(m>=0)
        return (*this)[m + l*l/4];
    if((-m)%2 == 0)
        return conj((*this)[-m + l*l/4]);
    else
        return -conj((*this)[-m + l*l/4]);
    //return (*this)[m+l + (l-1)*l/2];
}

/** \brief sum over m for a given l of the norms */
double BooData::getSumNorm(const size_t &l) const
{
    double sum=0.0;
    for(size_t m = 1; m <= l; m++)
        sum += norm((*this)(l,m));
    sum *= 2.0;
    sum += norm((*this)(l,0));
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
    double sumQl = getSumNorm(l);
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

    if( sumQl != 0) sumWl /= sqrt( sumQl * sumQl * sumQl);
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

    if( sumQl != 0) W /= pow(sumQl,1.5);
}

/** @brief get both Ql and Wl(real part)  */
void BooData::getInvarients(const size_t &l, double &Q, double &w) const
{
    complex<double> W(0.0,0.0);
    getInvarients(l,Q,W);
    w=W.real();
}

/** @brief Export the inner data to a String
  */
string BooData::toString() const
{
    ostringstream oss;
    for(size_t i=0;i<16;++i)
        oss<<(*this)[i]<<"\t";
    return oss.str();
}

/** @brief Constructor from string
  */
 BooData::BooData(const std::string &str): valarray< complex <double> >(complex <double>(0.0,0.0),16)
{
    istringstream iss(str);
    for(size_t i=0;i<16;++i)
        iss>>(*this)[i];
    return;
}

/** @brief Export to an array of 32 doubles and yeald a pointer to 32*sizeof(double) chars
  *
  */
char * BooData::toBinary(double* output) const
{
    for(size_t i=0;i<16;++i)
    {
        *(output+2*i) = (*this)[i].real();
        *(output+2*i+1) = (*this)[i].imag();
    }
    return (char*)output;
}

/** @brief constructor form a buffer of 32 doubles
  */
BooData::BooData(const double* buff) : std::valarray< std::complex <double> >(std::complex <double>(0.0,0.0),16)
{
	for(size_t i=0;i<16;++i)
		(*this)[i] = complex<double>(*(buff+2*i),*(buff+2*i+1));
	return;
}

/** \brief output to a stream */
ostream& operator<< (ostream& out, const BooData &boo )
{
    for(size_t l = 0; l <= 6; l+=2)
    {
        out<<"l="<<l;
        for(int m = -(int)l; m <= (int)l; m++)
            out << "\t" << boo(l,m);
        out <<endl;
    }
    return out;
}

/** \brief input from a stream */
istream& operator>> (istream& in, BooData &boo )
{
    string trash;
    for(size_t l = 0; l <= 6; l+=2)
    {
        in>>trash;
        for(size_t m=1; m <= l; m++)
            in >> trash;
        for(size_t m = 0; m <= l; m++)
            in >> boo(l,m);
    }
    return in;
}

/** \brief constructor from one bond */
BooFive::BooFive(const valarray<double> &rij): valarray< complex <double> >(complex <double>(0.0,0.0),6)
{
    valarray<double> diff(0.0,3);
	double theta, phi;

	// conversion to polar coordinates (physical convention, not mathematical)
    diff=rij/sqrt((rij*rij).sum());
    theta = acos(diff[2]);
	if(diff[0]==0.0)
	{
	    if(diff[1]==0.0) phi=0.0;
	    else phi = M_PI / 2.0 * (diff[1] > 0.0 ? 1.0 : -1.0);
	}
	else phi = atan( diff[1] / diff[0]) + (diff[0] > 0.0 ? 0.0 :M_PI);

	//fill in with spherical harmonics
    for(size_t m = 0; m <= 5; m++)
            (*this)[m] = boost::math::spherical_harmonic(5, m, theta, phi);
    return;
}

/** \brief sum over m of the norms */
double BooFive::getSumNorm() const
{
    double sum=0.0;
    for(size_t m = 1; m <= 5; m++)
        sum += norm((*this)[m]);
    sum *= 2.0;
    sum += norm((*this)[0]);
    return sum;
}

/** \brief Steindhardt order parameter Q5 */
double BooFive::getQ5() const
{
    return sqrt( 4.0 * M_PI * getSumNorm() / ( 2 * 5 + 1));
}
