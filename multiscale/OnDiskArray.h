/*
 * OnDiskArray.h
 *
 *  Created on: 27 sept. 2011
 *      Author: mathieu
 */

#ifndef ONDISKARRAY_H_
#define ONDISKARRAY_H_

#include <boost/iostreams/device/mapped_file.hpp>
#include <vector>
#include <string>

namespace Colloids {

class OnDiskArray : boost::noncopyable{
public:
	explicit OnDiskArray(const std::string &path, const std::vector<size_t> & dims, size_t element_size);
	virtual ~OnDiskArray();

	//accessors
	const boost::iostreams::mapped_file & get_file() const {return file;}
	const std::string & get_path() const {return path;};
	const std::vector<size_t> & get_dims() const {return dims;}

	//low level pixel-wise I/O
	const size_t get_offset(const std::vector<size_t> & pos) const;
	template<typename T> T read(const std::vector<size_t> & pos) const;
	template<typename T> void write(const std::vector<size_t> & pos, const T& value) const;

	//low level line-wise I/O
	const size_t get_line_offset(const std::vector<size_t> & pos) const;
	void readline(const std::vector<size_t> & pos, char* dest);
	void writeline(const std::vector<size_t> & pos, const char* src);


protected:
	boost::iostreams::mapped_file file;
	std::string path;
	std::vector<size_t> dims;
	std::vector<int> steps;
};

template<typename T>
T OnDiskArray::read(const std::vector<size_t> & pos) const
{
	return *(T*)(this->file.const_data()+this->get_offset(pos));
}

template<typename T>
void OnDiskArray::write(const std::vector<size_t> & pos, const T& value) const
{
	*(T*)(this->file.data()+this->get_offset(pos)) = value;
}

}

#endif /* ONDISKARRAY_H_ */
