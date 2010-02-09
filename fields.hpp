/**
    Copyright 2010 Mathieu Leocmach

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


 * \file fields.hpp
 * \brief Defines classes to manage fields of values : scalar, vectorial, static or dynamic
 * \author Mathieu Leocmach
 * \date 3 February 2010
 *
 */

#ifndef fields_H
#define fields_H

#include "index.hpp"
#include "traj.hpp"

namespace Colloids
{
	/**	\brief	A view of existing data organised as a field : one value per position */
    template<typename V>
    struct Field
    {
    	//typedef boost::ptr_vector<const V, boost::view_clone_allocator>::const_iterator const_iterator;
    	std::string name;
    	boost::ptr_vector<V, boost::view_clone_allocator> values;

    	Field(const std::string &name, const size_t size) : name(name), values(size){};
    	template<class ForwardInterator>
    	Field(ForwardInterator first, ForwardInterator last,const std::string &name="");
    	template<class tableInterator>
    	Field(tableInterator first, tableInterator last, const std::string &name, const size_t &column);

    	const double& operator[](const size_t &p) const {return this->values[p];}
    };
    typedef Field<double>							ScalarField;
    typedef Field<Coord>	                		VectorField;

    std::ostream& operator<< (std::ostream& os, const ScalarField &s);
	std::ostream& operator<< (std::ostream& os, const VectorField &s);

    /**	\brief	Contains data averaged over a given interval

		Averaging is made by forward scheme.
		Values are defined for all the range of the input,
		but the averaging interval vanishes when approaching the upper bound.
    */
    template<typename V>
    class DynamicField
    {
    	boost::ptr_vector< std::vector<V> > values;
    	std::deque< std::vector<size_t> > divisors;
    	const TrajIndex trajectories;
    	size_t averaging, actual_time;
    	const size_t size;

    	public:
			std::string name;

			DynamicField(const TrajIndex &trajectories, const size_t &averaging, const V& defaultValue, const std::string &name="") :
				name(name), trajectories(trajectories), divisors(averaging/2),
				values(trajectories.inverse.size()), size(trajectories.inverse.size()),
				averaging(averaging) {};
			DynamicField(const TrajIndex &trajectories, boost::ptr_vector< std::vector<V> > &values, const std::string &name="")
				: trajectories(trajectories), name(name), size(values.size()) {this->values.swap(values);};
			DynamicField(const DynamicField &d)
				: name(d.name), trajectories(d.trajectories), divisors(d.divisors),
				values(d.values), size(d.size),	averaging(d.averaging), actual_time(d.actual_time) {};

		void push_back(const Field<V> &frame);
		void assign(boost::ptr_vector< std::vector<V> > &values){this->values.swap(values);};
		Field<V> operator[](const size_t &t);
    };
    typedef DynamicField<double>	ScalarDynamicField;
    typedef DynamicField<Coord>		VectorDynamicField;

    /**	\brief	create a view of the content of an existing container */
    template<class V> template<class ForwardInterator>
	Field<V>::Field<V>(ForwardInterator first, ForwardInterator last,const std::string &name) : name(name), values(distance(first, last))
	{
		while(first!=last)
			values.push_back(&((*first++)));
	}

    /**	\brief	constructor from a table containing some columns and one field per column */
    template<class V> template<class tableInterator>
	Field<V>::Field<V>(tableInterator first, tableInterator last, const std::string &name, const size_t &column) : name(name), values(distance(first, last))
	{
		while(first!=last)
			values.push_back(&((*first++)[column]));
	}

	/** @brief bin a new frame into the average	  */
	template<class V>
	void DynamicField<V>::push_back(const Field<V> &frame)
	{
		//averaging on less and less time steps on the upper boundary
		if(actual_time + averaging/2+1 < values.size())
			divisors.push_back(std::vector<size_t>(trajectories.inverse[actual_time+averaging/2+1].size()));
		//bin the values at actual time in the next "average interval" time steps (including actual)
		for(size_t t=actual_time; t<min(actual_time+averaging/2+1, values.size()); ++t)
			for(size_t p=0; p<trajectories.inverse[t].size();++p)
			{
				const Traj& tr = trajectories[trajectories.inverse[t][p]];
				if(tr.exist(actual_time))
				{
					values[t][tr[actual_time]] = frame[p];
					divisors[t-actual_time][tr[actual_time]]++;
				}
			}
		//the time step "average interval" ago is totally binned. We divide it unless lower boundary
		if(actual_time>=averaging/2)
		{
			std::transform(
				values[actual_time-divisors.size()+1].begin(),
				values[actual_time-divisors.size()+1].end(),
				divisors.front().begin(),
				values[actual_time-divisors.size()+1].begin(),
				std::multiplies<double>()
				);
			divisors.pop_front(); //discard the divisors of the time step we just output
		}
		actual_time++;
	}
	/** @brief return a view of the field at frame t in terms of positions, not trajectories  */
	template<class V>
	Field<V> DynamicField<V>::operator[](const size_t &t)
	{
		return Field<V>(values[t].begin(), values[t].end(), name);
	}

};

#endif
