/*
 *      precondition.cpp
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
#include "precondition.h"

/**
 * Evaluate truth at each condition in the vector.
 * 
 * @param *v address to instance.  No checking of sizes because this should be fast.
 * @return boolean truth value
 */
bool Precondition::evaluate_truth(vector<short> *v){
	for(unsigned int i=0; i < conditions.size(); ++i){
		if(! conditions.at(i).evaluate(v) ){return false;}
	}
	return true;
}

/**
 * Push a condition that always evaluates to true onto the vector.
 */
Precondition::Precondition(){
	Condition c;
	c.attribute_index = -1;
	c.genotype_reference = -1; // fake
	c.genotype_conditional = Condition::GE;
	// This is always true.
	conditions.push_back(c);	
}
Precondition::~Precondition(){
	conditions.clear();
}

/**
 * Return the last condition on the vector.
 * 
 * @return Condition in last position on vector.
 */
Condition Precondition::last_condition(){
	return conditions.at(conditions.size()-1);
}

/**
 * Produce string representation of all conditions in vector.
 * 
 * @return string representation.
 */
string Precondition::to_string(DataAccess *data){
	string h = "";
	for(unsigned int i=0; i < conditions.size(); ++i){
		h += conditions.at(i).print(data);
		h += " " ;
	}
	return h;
}

string Precondition::hash(){
	string r = "";
	for(unsigned int i=0; i < conditions.size(); ++i){
		r += conditions.at(i).hash() + "_";
	}
	return r;
}

void Precondition::used(vector<long> &v){
	for(unsigned int i=0;i < conditions.size();i++){
		v.push_back(conditions.at(i).attribute_index);
	}
}
	
