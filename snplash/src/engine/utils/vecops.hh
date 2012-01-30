/*
 *      vecops.hh
 *      
 *      Copyright 2010 Richard T. Guy <guy@cs.toronto.edu>
 *      
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 2 of the License, or
 *      (at your option) any later version.
 *      
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *      
 *      You should have received a copy of the GNU General Public License
 *      along with this program; if not, write to the Free Software
 *      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *      MA 02110-1301, USA.
 */

#ifndef VECOPS_H
#define VECOPS_H

/**
 * Provide common vector operations.
 * 
 * @author Richard T. Guy
 */

#include <vector>

namespace vecops{
	
	/**
	 * vecAdd.  Add a double to a vector of doubles.
	 * 
	 * @param v A vector of type T.  The return is v .+ a.
	 * @param a A T.
	 */
	template<class T> 
	void vecAdd(std::vector<T> &v, T a){
		for(unsigned int i=0;i < v.size(); i++)
			v.at(i) += a;
	}
	
	/**
	 * vecDiv.  Divide a vector by a double.
	 * 
	 * @param v A vector of type T.  The return is v ./ a.
	 * @param a A T.
	 */
	template<class T> 
	void vecDiv(std::vector<T> &v, T a){
	for(unsigned int i=0;i < v.size(); i++)
		v.at(i) /= a;
	}
	
	/**
	 * Sum a vector of doubles.
	 * 
	 * @param v Sum the components of this.
	 * @return The sum.
	 */
	template<class T>
	T vecCumSum(const std::vector<T> &v){
		T d=0;
		for(unsigned int i=0;i < v.size(); i++)
			d += v.at(i);
		return d;
	}
	
	/**
	 * vecVariance.  Compute variance using E(X^2) - E(X)^2
	 */
	template<class T>
	T vecVariance(const std::vector<T> &v){
		T a = 0.0;
		for (unsigned int i=0; i < v.size(); i++)
			a += v[i]*v[i];
		a /= v.size();
		T b = vecCumSum(v);
		b /= v.size();
		
		return a - b*b;
	}
	
	// This one is still defined in .cpp
	std::vector<std::vector<double> > getDblVec(int size1, int size2);	
	
}

#endif
