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
    	boost::ptr_vector<const V, boost::view_clone_allocator> values;
    	Field(const std::string &name, const size_t size) : name(name), values(size){};

    	template<class tableInterator>
    	Field(tableInterator first, tableInterator last, const std::string &name, const size_t &column);
    };
    typedef Field<double>							ScalarField;
    typedef Field<Coord>	                		VectorField;

    std::ostream& operator<< (std::ostream& os, const ScalarField &s);
	std::ostream& operator<< (std::ostream& os, const VectorField &s);

    /**	\brief	Contains data averaged over a given interval */
    template<typename V>
    class DynamicField
    {
    	std::vector< std::vector<V> > values;
    	std::deque< std::vector<size_t> > divisors;
    	const TrajIndex &trajectories;
    	size_t averaging, actual_time;

    	public:
			std::string name;

			DynamicField(const std::string &name, const TrajIndex &trajectories, const size_t &averaging, const V& defaultValue) :
				name(name), trajectories(trajectories),
				values(trajectories.getMaxTime()-averaging, std::vector<V>(trajectories.size(), defaultValue)),
				averaging(averaging) {};

		void push_back(const Field<V> &frame);
		Field<V> operator[](const size_t &t) const;
    };
    typedef DynamicField<double>	ScalarDynamicField;
    typedef DynamicField<Coord>		VectorDynamicField;

    /**	\brief	constructor from a table containing some columns and one field per column */
    template<class V> template<class tableInterator>
	Field<V>::Field<V>(tableInterator first, tableInterator last, const std::string &name, const size_t &column) : name(name), values(distance(first, last))
	{
		while(first!=last)
			values.push_back(&((*first++)[column]));
	};

	/** @brief bin a new frame into the average	  */
	template<class V>
	void DynamicField<V>::push_back(const Field<V> &frame)
	{
		divisors.push_back(std::vector<size_t>(trajectories.size()));
		for(size_t t=actual_time; t<min(actual_time+averaging, values.size()); ++t)
			for(size_t p=0; p<trajectories.inverse[t].size();++p)
			{
				values[t][trajectories.inverse[t][p]]=frame[p];
				divisors[t-actual_time][trajectories.inverse[t][p]]++;
			}

		if(actual_time>=averaging)
		{
			for(std::vector<size_t>::const_iterator tr=trajectories.inverse[actual_time-averaging].begin();
				tr!=trajectories.inverse[actual_time-averaging].end();++tr)
			{
				values[actual_time-averaging][*tr] /= (double)divisors.front()[*tr];
			}
			divisors.pop_front();
		}
		actual_time++;
	}
	/** @brief return a view of the field at frame t in terms of positions, not trajectories  */
	template<class V>
	Field<V> DynamicField<V>::operator[](const size_t &t) const
	{
		Field<V> view(name, trajectories.inverse[t].size());
		for(size_t p=0; p<trajectories.inverse[t].size();++p)
			view.values.push_back(&values[t][trajectories.inverse[t][p]]);
		return view;
	}

};

#endif
