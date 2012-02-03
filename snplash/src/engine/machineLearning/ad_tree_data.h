/*
 *      ad_tree_data.h
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
#ifndef AD_DATA_H
#define AD_DATA_H
/**
 * Holds a single tree.  Contains a vector of nodes and methods to access information about the tree.
 */

#include "ad_rule.h"
#include "precondition.h"
#include <map>		// Gives us multimap.
#include <stdlib.h> // so we can exit!

#include <algorithm> // so we can sort

using namespace std;

class AD_Data {

	public :

		~AD_Data();

		/* Functions */
		// Evaluation methods.
		double evaluate_instance(vector<short> *); // Perform evaluation of this instance.
		double evaluate(vector<short> *, int ); // holds start index and vector to evaluate.
		double evaluate_on_last(vector<short> * v){ return evaluate(v,nodes.size()-1); }
		bool evaluate_truth(vector<short> *, int) ; // return only truth value.
		
		void push_node(AD_Rule , unsigned int);
		void push_precondition(Precondition);

		AD_Rule* node_at(int i){return &nodes.at(i);}
		Precondition* precon_at(int i){return &preconditions.at(i);}
		string print(DataAccess *data);
		
		int node_size(){return nodes.size();}
		int precondition_size(){return preconditions.size();}

		//int merged(AD_Rule *, unsigned int);
		void report(map<string, int> &hash , map<string, AD_Rule> &key);
		void report_leaves(map<string, int> &hash , map<string, AD_Rule> &key);
		void reset();

	private :
		vector<AD_Rule> nodes;
		vector<Precondition> preconditions; // These hold conjunctions of rules.
		vector<unsigned int> parent; // Each spot holds a location of the parent of that node.

};

#endif
