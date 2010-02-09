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


 * \file files_series.hpp
 * \brief Defines classes for time series and z series of files
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 15 December 2008
 *
 * Define functions and classes relative to the files series for the particle tracking code
 *
 */

#ifndef file_series_H
#define file_series_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/format.hpp>

namespace Colloids
{
    class FileSerie
    {
        boost::format pattern;
        size_t length, offset;

        public:
            FileSerie(const std::string &filesPattern, const std::string &token, size_t size, size_t offset=0);
            std::string operator%(const size_t &step){return (pattern%(step+offset)).str();}
            size_t size(){return length;}
    };
}

class TokenTree
{
    public:
        size_t nbdigit,level;
        std::vector<std::string> tokens;
        TokenTree *prefix, *postfix;
        std::string value;
        boost::format formatPattern;

        TokenTree(){return;};
        TokenTree(const TokenTree &source);
        TokenTree(const std::vector<std::string> &toks,const std::string &pattern,const size_t &lev=0);
        //~TokenTree();

        const std::string& token() const;
        std::string operator()(std::vector<size_t> &v) const;
        std::string getDigits(std::vector<size_t> &v) const;
        std::string getPrefix() const;
        std::string getNoIndexName() const;
        std::string getNoIndexNameNoExt() const;
        std::string getPath() const;
        std::string getPattern(const std::string &s="") const;
        std::string getPatternNoPath() const;

        std::string getDigitsFormat() const;
        std::string getFormatPattern() const;
        boost::format &operator%(size_t n);
};

std::ostream& operator<< (std::ostream& out, const TokenTree &t );

inline const std::string& TokenTree::token() const
{
    return tokens[level];
}

#endif
