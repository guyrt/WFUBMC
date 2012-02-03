/*
 *      ad_rule.cpp
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
#include "ad_rule.h"

AD_Rule::AD_Rule(){}

/**
 *	Set up the rule with the correct instance information.
 */
AD_Rule::AD_Rule(Precondition p, Condition c, double a1, double a2){

	this->precon = p;
	this->con = c;
	this->score_true = a1;
	this->score_false = a2;
}


/**
 * Perform a test of this rule on the input data passed in.
 */
double AD_Rule::evaluate(vector<short> *vec){
	if(precon.evaluate_truth(vec)){
		if(con.evaluate(vec)){
			return score_true;
		}else{return score_false;}
	}else{return 0;}
}

/**
 * Evaluate truth value of this rule on the instance vector.
 * 
 * name: evaluate_truth
 * @param vector<short> address representing instance.  No checks on correctness of length.
 * @return boolean truth value.
 */
bool AD_Rule::evaluate_truth(vector<short> * vec){
	return con.evaluate(vec) && precon.evaluate_truth(vec);
}

bool AD_Rule::operator< (const AD_Rule &compare) const{
	
	Precondition compare_p = compare.getPrecondition();
	
	bool localIsLarger = false;
	unsigned int sz = precon.conditions.size();
	if (sz > compare_p.conditions.size()){
		localIsLarger = true;
		sz = compare_p.conditions.size();
	}
	
	for (unsigned int i=0; i < sz; i++){
		if (precon.conditions.at(i).attribute_index != compare_p.conditions.at(i).attribute_index){
			return precon.conditions.at(i).attribute_index < compare_p.conditions.at(i).attribute_index;
		}
	}
	
	return !localIsLarger; // reverse this, so if local is the shorter one, it is returned first.
}

/**
 * Returns a string representation of this node based on the precondition
 * condition and two scores.
 * 
 * @return string representation. 
 */
string AD_Rule::to_string(DataAccess *data){

	string pr = "";
	for(unsigned int i=1;i < precon.conditions.size();i++){
		pr += precon.conditions.at(i).print(data);
	}

	stringstream s;
	
	if (pr.length() > 0){
		s << pr << ": " ;
	}
	
	s << con.print(data) << " " << score_true;
	if (con.attribute_index >= 0 ){
		s << endl;
		if (pr.length() > 0)
			s << pr << ": " ;
		s << con.print_reverse(data) << " " << score_false;
	}
	return s.str();

}

/**
 * Create a string object representing this rule as a hash (so unique)
 * 
 * @return string object containting the hash
 */
string AD_Rule::hash(){
	stringstream s;
	s << precon.hash() << con.hash();
	return s.str();
}

/**
 * Create an array with list of the nodes used.
 * 
 * @param vector<long> reference filled with positions in array. 
 */
void AD_Rule::used(vector<long> &v){
	precon.used(v);
	v.push_back(con.attribute_index);
}

/**
 * Check for equality between two rules.  Equality is based on equality
 * of connection only, NOT equality of preconditions.
 * 
 * name: ==
 * @param const AD_Rule& rule to compare to
 * @return int: 0 means not equal, 1 means equal.
 */
/**int AD_Rule::operator==(const AD_Rule& right){

	if( con == right.con ){
		cout << con.attribute_index << endl;
		cout << right.con.attribute_index << endl;
		return 1;
	}
	else{
		cout << con.attribute_index << endl;
		cout << right.con.attribute_index << endl;
		return 0;
	}
}*/
