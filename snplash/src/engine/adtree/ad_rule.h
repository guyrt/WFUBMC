/*
 *      ad_rule.h
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

#ifndef AD_RULE_H
#define AD_RULE_H

/**
 * A rule for an ADTree node.
 *
 * Contains
 * 	pointer to precondition (another rule)
 * 	condition
 * 	scores for true and false
 * 
 * Functions:
 * 	Evaluate the rule on a single individual and return either numeric
 * score or boolean truth value.
 *
 */

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include "condition.h"
#include "precondition.h"

using namespace std;

class AD_Rule {

	public :

		AD_Rule();
		AD_Rule(Precondition, Condition , double a1, double a2);

		double evaluate(vector<short> *);
		bool evaluate_truth(vector<short> * vec);
		long attribute_value(){return con.attribute_index;}

		string to_string(DataAccess *data);
		string hash();
		void used(vector<long> &);

		bool operator< (const AD_Rule &compare) const;
		
		Precondition const& getPrecondition() const {return precon;}

	private :
		Condition con;
		Precondition precon;
		double score_true, score_false;

};

#endif
