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

/**	\brief constructor breaking the string recursively using the tokens */
TokenTree::TokenTree(const std::vector<string> &toks,const std::string &pattern,const size_t &lev)
{
    nbdigit=0;
    level = lev;
    tokens = toks;
    if(level<toks.size())
    {
        string token = toks[level];
        const size_t end_prefix = pattern.rfind(token);
        if(end_prefix!=string::npos)
        {
            //const int end_prefix = end_pre_z+token.size();
            prefix = new TokenTree(toks,pattern.substr(0,end_prefix),level+1);
            while(isdigit(pattern[end_prefix+token.size()+nbdigit]))
                nbdigit++;

            postfix = new TokenTree(toks,pattern.substr(end_prefix+token.size()+nbdigit),level+1);
        }
        else
        {
            prefix = new TokenTree(toks,pattern,level+1);
            postfix = new TokenTree(toks,"",level+1);
            //value = filename;
        }
    }
    else
        value = pattern;

    formatPattern = boost::format(getFormatPattern());
    return;
}

/**	\brief copy constructor	*/
TokenTree::TokenTree(const TokenTree &source)
{
    nbdigit = source.nbdigit;
    level = source.level;
    value = source.value;
    tokens = source.tokens;
    if(level<tokens.size())
    {
        prefix = new TokenTree(*prefix);
        postfix = new TokenTree(*postfix);
    }
}

/**	\brief get the digits string at the current level, with token and preceding zeros
	Example : the pattern is foo_t000_z00.tif and v=[22,47]
	At level 0 getDigits(v) => "_t022"
	At level 1 getDigits(v) => "_z47"
*/
string TokenTree::getDigits(std::vector<size_t> &v) const
{
    if(level>=tokens.size() || nbdigit==0) return "";
    ostringstream os;
    os << v[level];
    string str(os.str());
    for(size_t i = str.length();i<nbdigit;++i)
        str = "0"+ str;
    return token()+str;
}

/**	\brief Get recursively the pattern filled with the good numbers
	Example : the pattern is foo_t000_z00.tif, tokens=["_t","_z"] and v=[22,47]
	At level 2 => "foo_t" and "_z" and ".tif"
	At level 1 => "foo_t" and "_z47.tif"
	At level 0 => "foo_t022_z47.tif"
*/
string TokenTree::operator()(std::vector<size_t> &v) const
{
    if(level>=tokens.size()) return value;
    return (*prefix)(v)+ getDigits(v) + (*postfix)(v);
}

/** \brief get recursively the absolute prefix without the token	*/
string TokenTree::getPrefix() const
{
    if(level<tokens.size())	return prefix->getPrefix();
    return value;
}

/**	\brief	get the concatenation of everything but tokens and digits
	Example : The pattern foo_t000_foo_z00_bar.tif gives foo_foo_bar.tif
*/
string TokenTree::getNoIndexName() const
{
    if(level<tokens.size())
    {
        const string pref = prefix->getNoIndexName();
        return pref.substr(0,pref.size()-token().size())+postfix->getNoIndexName();
    }
    return value;
}

/** \brief get the concatenation of everything but tokens, digits and the extension (what is after the last ".") */
string TokenTree::getNoIndexNameNoExt() const
{
    const string noIndex = getNoIndexName();
    return noIndex.substr(0,noIndex.find_last_of("."));
}

/** @brief get the path out of the prefix (useful only at level 0)
  */
string TokenTree::getPath() const
{
	const string pref = getPrefix();
	const size_t endfolder = pref.find_last_of("/\\");
	if(endfolder==pref.npos)
		return "./";
	return pref.substr(0,endfolder+1);
}

/** @brief get the pattern with all digits at 0  */
string TokenTree::getPattern(const string &s) const
{
	vector<size_t> v(tokens.size()-level,0);
	return (*this)(v).insert(getPrefix().size(),s);
}

/** @brief get the pattern with all digits at 0 without the path  */
string TokenTree::getPatternNoPath() const
{
	const string pattern = getPattern();
	const size_t endfolder = pattern.find_last_of("/\\");
	return pattern.substr(endfolder+1);
}


/** \brief output recursively the structure of the token tree for debug */
std::ostream& operator<< (std::ostream& out, const TokenTree &t )
{
    string tab(t.level,'-');
    if(t.nbdigit!=0 || t.level<t.tokens.size())
    {
        //cout << t.level <<endl;
        out<<tab<< "prefix" << t.token() << " : " << *t.prefix << endl;
        out<<tab << "nbdigit" << t.token() << " : " << t.nbdigit <<endl;
        out<<tab << "postfix" << t.token() << " : " << *t.postfix << endl;
    }
    else
        out<< t.value;

    return out;
}

/**	\brief get the format for the digits at the current level, with token and preceding zeros
	Example : the pattern is foo_t000_z00.tif
	At level 0 getDigitsFormat() => "_t%|03d|"
	At level 1 getDigitsFormat() => "_z%|02d|"
*/
string TokenTree::getDigitsFormat() const
{
    if(level>=tokens.size() || nbdigit==0) return "";
    ostringstream os;
    os << token()<<"%|0"<<nbdigit<<"d|";
    return os.str();
}

/** @brief get the boost::format pattern of this branch  */
string TokenTree::getFormatPattern() const
{
    if(level>=tokens.size()) return value;
    return prefix->getFormatPattern()+ getDigitsFormat() + postfix->getFormatPattern();
}

/** @brief access and first argument assignment of the internat boost::format pattern
    cout << tt(foo_t000_z00.bar) % 42 % 30 => print "foo_t042_z30.bar"
    (tt(foo_t000_z00.bar) % 42 % 30).str() => yealds the string "foo_t042_z30.bar"
  */
boost::format & TokenTree::operator%(size_t n)
{
    return formatPattern % n;
}


