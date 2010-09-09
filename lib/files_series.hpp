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

#include <boost/format.hpp>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

namespace Colloids
{
    class FileSerie
    {
        boost::format pattern;
        std::string token;
        size_t length, offset;

        public:
            FileSerie(const std::string &filesPattern, const std::string &token, size_t size, size_t offset=0);
            std::string operator%(const size_t &step){this->pattern.clear(); return (pattern%(step+offset)).str();}
            size_t size() const {return length;}
            size_t get_offset() const {return offset;}
            FileSerie changeExt(const std::string &ext) const;
            FileSerie addPostfix(const std::string &postfix) const;
            FileSerie addPostfix(const std::string &postfix, const std::string &ext) const;
            std::string head() const;

            static std::string get0th(const std::string &prefix, const size_t &size, const std::string &token="_t");
    };
}

#endif
