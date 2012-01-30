/*
 *      condition.h
 *      
 *      Copyright 2009 Richard T. Guy <guyrt7@wfu.edu>
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
#ifndef CONDITION_H
#define CONDITION_H

/**
 * A rule for an ADTree node.
 * 
 * Contains
 * 	pointer to precondition (another rule)
 * 	condition 
 * 	scores for true and false 
 * 
 */

#include "../data_plugin.h"
#include <string>
#include <vector>
#include <sstream>

using namespace std;

class Condition {
	
	public :
	
		Condition();
		
		enum comparison {GE, GT, LE, LT, EQ, NE};
	
		long attribute_index;		// This holds position in the data array.
		comparison genotype_conditional; 
		short genotype_reference;  // will be either 1,2,3, or 4

		// Key:
		//   1 : >=
		//	 2 : >
		//	 3 : <=
		//	 4 : <
		//	 5 : ==
		//	 6 : !=
		
		bool evaluate(vector<short> *);
		
		void inverse();
		string print(DataAccess *data);
		string print_reverse(DataAccess *data);
		string hash();
		
		int operator==(const Condition& right);
};

#endif

