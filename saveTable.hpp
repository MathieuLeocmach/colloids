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


 * \file saveTable.hpp
 * \brief Defines utility function to save table
 * \author Mathieu Leocmach
 * \version 0.1
 * \date 16 April 2009
 *
 */

#ifndef save_table_H
#define save_table_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
/** \brief saving a range of colums into a single file */
template <class InputIterator>
void saveTable(InputIterator first, InputIterator last,const std::string &filename,const std::string &header,const double &scale=1.0)
{
     //std::cout << "export to " << filename << std::endl;
     std::ofstream output(filename.c_str(), std::ios::out | std::ios::trunc);
     if(output)
     {
       //header line
       int slash = filename.find_last_of("/");
       if(slash == (int)std::string::npos) slash=-1;
       output << "#"<<header << std::endl;

       for(size_t r=0; r<(*first).size();++r)
       {
         InputIterator it = first;
         output << r*scale;
         while ( it!=last )
         {
			output << "\t" << (*it)[r];
			++it;
         }
         output << std::endl;
       }
       output.close();
     }
     else
        throw std::invalid_argument("No such file as "+filename);
}
#endif
