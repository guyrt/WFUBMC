/*
 *      precondition.h
 *      
 *      Copyright 2009 Richard T. Guy <guy@cs.toronto.edu>
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
#ifndef PRECONDITION_H
#define PRECONDITION_H

#include "condition.h"
#include <vector>

class Precondition {
	
	public : 
		Precondition();
		~Precondition();
		vector<Condition> conditions;
		
		Condition last_condition();
		bool evaluate_truth(vector<short> *);	
		string to_string(DataAccess *data);
		string hash();
		void used(vector<long> &);
	
};

#endif
