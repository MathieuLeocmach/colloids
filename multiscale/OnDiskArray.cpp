/*
 * OnDiskArray.cpp
 *
 *  Created on: 27 sept. 2011
 *      Author: mathieu
 */

#include "OnDiskArray.h"
#include <numeric>
#include <iostream>
#include <stdexcept>

namespace Colloids {

OnDiskArray::OnDiskArray(const std::string &path, const std::vector<size_t> & dims, size_t element_size) :
		path(path), dims(dims)
{
	const size_t total_size = std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>());
	std::cout<<path<<std::endl;
	boost::iostreams::mapped_file_params params(path);
	params.new_file_size = total_size * element_size;
	params.flags = boost::iostreams::mapped_file::readwrite;
	this->file.open(params);
	this->steps.resize(dims.size(), element_size);
	//C order: larger dimension index => smaller step
	for(int d=dims.size()-1; d>0; --d)
		this->steps[d-1] = this->steps[d] * this->dims[d];
}

/**
 * \brief The file is deleted from disk when the object is destroyed
 */
OnDiskArray::~OnDiskArray() {
	remove(this->path.c_str());
}

const size_t OnDiskArray::get_offset(const std::vector<size_t> & pos) const
{
#ifndef NDEBUG
	if(this->steps.size() != pos.size())
		throw std::invalid_argument("OnDiskArray::get_offset: Too many or too few coordinates");
#endif
	return std::inner_product(this->steps.begin(), this->steps.end(), pos.begin(), 0);
}

const size_t OnDiskArray::get_line_offset(const std::vector<size_t> & pos) const
{
#ifndef NDEBUG
	if(this->steps.size() != pos.size()+1)
		throw std::invalid_argument("OnDiskArray::get_line_offset: Too many or too few coordinates");
#endif
	return std::inner_product(this->steps.begin(), this->steps.end()-1, pos.begin(), 0);
}

void OnDiskArray::readline(const std::vector<size_t> & pos, char *dest)
{
	const char * s = this->file.data() + this->get_line_offset(pos);
	for(size_t i=0; i<this->dims.back()*this->steps.back(); ++i)
		*dest++ = *s++;
}

void OnDiskArray::writeline(const std::vector<size_t> & pos, const char *src)
{
	char * d = this->file.data() + this->get_line_offset(pos);
	for(size_t i=0; i<this->dims.back()*this->steps.back(); ++i)
		*d++ = *src++;
}

}
