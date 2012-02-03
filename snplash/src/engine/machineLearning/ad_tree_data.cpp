/*
 *      ad_tree_data.cpp
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
#include "ad_tree_data.h"

AD_Data::~AD_Data(){
	nodes.clear(); // we made these...
}

/**
 * Return score for a given node start on vector v.
 * 
 * @param vector<short> pointer to instance
 * @param int denoting position of the node to be evaluated.
 * @return double value of node i on instance v 
 */
double AD_Data::evaluate(vector<short> *v, int start){

	return nodes.at(start).evaluate(v);

}

/**
 * Evaluate the truth value (boolean) of an instance on a single rule.
 * 
 * @param vector<short> pointer to instance
 * @param int denoting position of the node to be evaluated.
 * @return boolean value of rule i evaluated on instance v.
 */
bool AD_Data::evaluate_truth(vector<short> *v, int i){
	return nodes.at(i).evaluate_truth(v);
}

/**
 * Push a single node onto the nodes vector and record its parent.
 *
 * @param AD_Rule to be added to the tree.
 * @param int that denotes position of the parent.
 * @return none 
 */
void AD_Data::push_node(AD_Rule r, unsigned int i){
	if(i > parent.size() && nodes.size() > 0){
		cerr << "Error pushing new rule: asked for out of range: " << i << " " << parent.size() << endl;
		exit(0);
	}
	nodes.push_back(r);
	parent.push_back(i);
}

/**
 * Adds a precondition to the preconditions vector, which is private.
 *
 * @param Precondition object to be added.
 * @return none
 */
void AD_Data::push_precondition(Precondition r){
	preconditions.push_back(r);
}

/**
 * Create a formatted string ready for printing.
 *
 * @return Return formatted string.
 */
string AD_Data::print(DataAccess *data){

	string tree;

	// Sort by first attribute.
	sort(nodes.begin()+1, nodes.end());

	for (unsigned int i=0;i < nodes.size();++i){
		tree += nodes.at(i).to_string(data);
		tree += "\n";
	}

	return tree;
}

/**
 * Evaluate the instance (given in the vec) on this tree.
 *
 * @param vector pointer to instance data
 * @return double score for this individual.
 */
double AD_Data::evaluate_instance(vector<short> *vec){
	double return_value = 0;
	for(unsigned int i=0;i < nodes.size(); ++i){
		return_value += nodes.at(i).evaluate(vec);
	}
	return return_value;
}

/**
 * @param map<string, int> reference that is updated with information about the SNPs in this map.
 * @param map<string, AD_Rule> reference tying hashs to Rules.
 */
void AD_Data::report(map<string, int> &hash , map<string, AD_Rule> &key){
	
	string h;
	map<string, int> seenIt;
	for(unsigned int i=0;i < nodes.size(); ++i){
		h = nodes.at(i).hash();
		if(h.compare("-1>=-1_-1>=-1")!=0){ // Ignore the occurence of the trivial rule.
			if(seenIt.count(h) < 1){
				hash[h]++;
				key[h] = nodes.at(i);
				seenIt[h]++;
			}
		}
	}
	
}

/**
 * Report only rules that are leaf rules.
 * 
 * @param map<string, int> reference that is updated with information about the SNPs in this map.
 * @param map<string, AD_Rule> reference tying hashs to Rules.
 */
void AD_Data::report_leaves(map<string, int> &hash , map<string, AD_Rule> &key){
	
	vector<bool> leaf_node;
	leaf_node.resize(nodes.size(), true);
	for(unsigned int i=0;i < nodes.size(); ++i){
		leaf_node.at(parent.at(i)) = false;
	}
	
	string h;
	for(unsigned int i=0;i < nodes.size(); ++i){
		h = nodes.at(i).hash();
		if(h.compare("-1>=-1_-1>=-1")!=0 && leaf_node.at(i)){ // Ignore the occurence of the trivial rule.
			hash[h]++;
			key[h] = nodes.at(i);
		}
	}
}

/**
 * Clear out all containers.
 * 
 * @param none
 * @return none
 */
void AD_Data::reset(){
	
	nodes.clear();
	preconditions.clear();
	parent.clear();
	
}
