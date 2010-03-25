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

#include "files_series.hpp"
#include <stdexcept>

using namespace std;
using namespace Colloids;

/** @brief Constructor  */
FileSerie::FileSerie(const std::string &namePattern, const std::string &token, size_t size, size_t offset) :
    length(size), offset(offset), token(token)
{
    string head = namePattern, tail;
	size_t digits=0, pos = namePattern.rfind(token);

    if(pos==string::npos)
        throw invalid_argument(("Name pattern \""+namePattern+"\" doesn't contain token \""+token+"\"").c_str());

    head.resize(pos);
    digits = namePattern.find_first_not_of("0123456789", pos+token.size());
    tail = namePattern.substr(digits);
    digits -= pos+token.size();

    ostringstream os;
    os<<head<<token<<"%|0"<<digits<<"|"<<tail;
    try
    {
		this->pattern.parse(os.str());
    }
    catch(exception &e)
    {
    	cerr<<"head: "<<head<<endl;
    	cerr<<"digits: "<<digits<<endl;
    	cerr<<"tail: "<<tail<<endl;
    	cerr<<"all: "<<os.str()<<endl;
    	throw;
    }
}

/** @brief return a new FileSerie with different extension  */
FileSerie FileSerie::changeExt(const std::string &ext) const
{
	boost::format fmt = this->pattern;
	fmt.clear();
	string pat = (fmt % 0).str();
	pat = pat.substr(0, pat.rfind(".")) + ext;
	return FileSerie(pat, this->token, this->size(), this->offset);
}

/** @brief return a new FileSerie with a postfix added just before the token  */
FileSerie FileSerie::addPostfix(const std::string &postfix) const
{
	boost::format fmt = this->pattern;
	fmt.clear();
	string pat = (fmt % 0).str();
	pat.insert(pat.rfind(this->token), postfix);
	return FileSerie(pat, this->token, this->size(), this->offset);
}

/** @brief return a new FileSerie with a postfix added just before the token and a different extension  */
FileSerie FileSerie::addPostfix(const std::string &postfix, const std::string &ext) const
{
	boost::format fmt = this->pattern;
	fmt.clear();
	string pat = (fmt % 0).str();
	pat.insert(pat.rfind(this->token), postfix);
	pat = pat.substr(0, pat.rfind(".")) + ext;
	return FileSerie(pat, this->token, this->size(), this->offset);
}

/** @brief return the head of the serie, without time dependence, with the given extension  */
string FileSerie::head() const
{
	boost::format fmt = this->pattern;
	fmt.clear();
	string pat = (fmt % 0).str();
	pat.erase(pat.rfind(this->token));
	return pat;
}



