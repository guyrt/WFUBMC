/*
 *      bagging.cpp
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

#include "bagging.h"

Bagging::Bagging(){
	param_reader = ParamReader::Instance();
	data = new DataAccess;
	data->init(NULL);
	bag_param = new EngineParamReader;
	haveOwner = false;
	order_in_bag = 0;
}

Bagging::Bagging(DataAccess *d){
	this->param_reader = ParamReader::Instance();
	this->data = d;
	this->bag_param = new EngineParamReader;
	this->haveOwner = false;
	order_in_bag = 0;
}

/*
 * If we have ownership of our own properties, then delete them.
 * Ownership flows up: if this object is owned by another, then do
 * not worry about garbage collection.
 */
Bagging::~Bagging(){
	if(!haveOwner) delete_my_innards();
}

/**
 * Perform garbage collection.
 */
void Bagging::delete_my_innards(){
	this->engines.clear();
	delete param_reader;
	param_reader = NULL;
	delete data;
	data = NULL;
	delete bag_param;
	bag_param = NULL;
	delete reader;
	reader = NULL;
}

/**
 *	Retrieve data from files, prep data.
 * 
 * Only called if this is not a slave.
 */
void Bagging::init(){

	bag_param->read_parameters(param_reader->get_engine_specific_params());

	initializeReader(); // defined in engine.h
	
	this->outstream.open(this->param_reader->get_out_file().c_str() );
	if(!this->outstream.is_open()){
		cerr << "Unable to open file " << this->param_reader->get_out_file() << endl;
		exit(1);
	}

	reader->process(data->getDataObject(), param_reader);
	data->getDataObject()->prep_data(param_reader);
	data->make_categorical(-1,1);
	data->getDataObject()->remove_haplotype();
}

/**
 * Delete current adparam reader and replace it with the referant.
 */
void Bagging::enslave(EngineParamReader *b){
	bag_param = b;
	haveOwner = true;
}

/**
 * Create bootstrap samples for every set by calling the snp_data set.
 * Right now we are always using the same size.
 *
 */
void Bagging::preProcess(){
	
//	int bag_num;
	/*engines.reserve(bag_param->get_number_of_bags());

	ADTree *a;
	DataAccess *d;
	d = new DataAccess();
	d->init(data->getDataObject());
	a = 0;
	a = new ADTree(param_reader,d);

	if(bag_param->get_engine_type() == ParamReader::ADTREE){
		for(int i=0; i < bag_param->get_number_of_bags(); i++){
			
			engines.push_back(*a);

		}
		for(int i=0;i < bag_param->get_number_of_bags(); ++i){
			engines.at(i).enslave(bag_param);
		}
		
	}else{
		cerr << "You must use ADTREE reader." << endl;
		return;
	}*/
	
	engines.reserve(bag_param->get_number_of_bags());
	
	if(bag_param->get_engine_type() == ParamReader::ADTREE){
		
		for(int i=0; i < bag_param->get_number_of_bags(); i++){
			
			ADTree *a;
			DataAccess *d;
			d = new DataAccess();
			d->init(data->getDataObject());
			a = new ADTree(d);
			engines.push_back(*a);

		}
		for(int i=0;i < bag_param->get_number_of_bags(); ++i){
			engines.at(i).enslave(bag_param);
		}
		
	}else{
		cerr << "You must use ADTREE reader." << endl;
		return;
	}
}

/**
 *
 */
void Bagging::process(){
	int num_bags = bag_param->get_number_of_bags();
	int iteration = 1;
	do{
		//if(iteration++ > 3) break; // experimental
		#pragma omp parallel
		{
			#pragma omp for schedule(static) nowait
			for(int i=0;i < num_bags;i++){
				this->engines.at(i).preProcess();
				this->engines.at(i).process();
			}
		} // end parallel
		// Examine then repeat.
	}
	while(evaluate_trees2());
	
	outstream.close();
}

void Bagging::test(){

}

/*
 * Evaluate the set of trees for structural motifs that we want to keep.
 *
 * Right now, the only available option is to evaluate rules, which include
 * the precondition.
 *
 * @return boolean true if something was removed.
 */
bool Bagging::evaluate_trees(){

	map<string, int> hash_count; // Hold number of occurences of each rule.
	map<string, AD_Rule> rule_map; // Hold key of hash to rule.
	map<string, int>::iterator it;
	
	bool ret_val = false;
	
	if(this->bag_param->get_bag_type() == EngineParamReader::ALL){
		for(int i=0;i < this->bag_param->get_number_of_bags();i++){
			this->engines.at(i).report(hash_count, rule_map);
		}
	}else if(this->bag_param->get_bag_type() == EngineParamReader::LEAVES){
		for(int i=0;i < this->bag_param->get_number_of_bags();i++){
			this->engines.at(i).report_leaves(hash_count, rule_map);
		}
	}
	
	outstream << "Run completed: " << endl;
	
	// Get maximum number of times an element occurs:
	int max_occurs = 0;
	for(it = hash_count.begin();it != hash_count.end(); it++){
		if((*it).second > max_occurs){
			max_occurs = (*it).second;
		}
		if((*it).second >= this->bag_param->get_threshold_of_bags()){
			print_and_process((*it).second, rule_map[(*it).first]);
			ret_val = true;
		}
	}
	cout << "Maximum occurences: " << max_occurs << " out of " << this->bag_param->get_number_of_bags() << " bags." << endl;
	if(param_reader->get_verbosity() >= 2){
		for(it = hash_count.begin();it != hash_count.end(); it++){
			if((*it).second == max_occurs){
				cout << "Rule: " << (*it).first << endl;
			}
		}
	}
	
	data->getDataObject()->snp_flush();
	hash_count.clear();
	rule_map.clear();
	return ret_val;
}

/**
 * Evaluate with a more accurate method that reorders SNPs.
 */
bool Bagging::evaluate_trees2(){
	map<string, int> hash_count; // Hold number of occurences of each rule.
	map<string, AD_Rule> rule_map; // Hold key of hash to rule.
	map<string, int>::iterator it;
	map<string, AD_Rule>::iterator rit;
	
	bool ret_val = false;
	
	if(this->bag_param->get_bag_type() == EngineParamReader::ALL){
		for(int i=0;i < this->bag_param->get_number_of_bags();i++){
			this->engines.at(i).report(hash_count, rule_map);
		}
	}else if(this->bag_param->get_bag_type() == EngineParamReader::LEAVES){
		for(int i=0;i < this->bag_param->get_number_of_bags();i++){
			this->engines.at(i).report_leaves(hash_count, rule_map);
		}
	}
	
	map<vector<long>, int> counts;
	
	/* Bug here: make sure doesn't happen more than once per tree! */
	for(rit = rule_map.begin(); rit != rule_map.end(); rit++){
		vector<long> u;
		(*rit).second.used(u);
		int num = hash_count[(*rit).first];
		sort(u.begin(), u.end());
		u.erase(u.begin(), u.begin()+1);
		counts[u]+=num;
	}
	
	// Print.
	map<vector<long>, int>::iterator vit;
	outstream << "New run: " << endl;
	for(vit = counts.begin(); vit != counts.end(); vit++){
		if((*vit).second >= this->bag_param->get_threshold_of_bags()){
			print_and_process((*vit).second, (*vit).first);
			ret_val = true;
		}
	}
	
	/// This path will remove one at a time.
	/*vector<long> maxPath; int max = 0;
	for(vit = counts.begin(); vit != counts.end(); vit++){
		if((*vit).second > max){
			max = (*vit).second;
			maxPath = (*vit).first;
		}
	}
	if(max >= this->bag_param->get_threshold_of_bags()){
		print_and_process(max, maxPath);
		ret_val = true;
	}*/
	
	
	data->getDataObject()->snp_flush();
	hash_count.clear();
	rule_map.clear();
	return ret_val;
}

/**
 * Given a rule and number of hits, print output.
 */
void Bagging::print_and_process(int bags, vector<long> SNPs){
	
	outstream << "Feature: " << bags << " hits out of " << this->bag_param->get_number_of_bags() << " possible ";
	for(unsigned int ii=0;ii <SNPs.size();ii++){
		outstream << " " << data->snp_name(SNPs.at(ii));
		this->data->getDataObject()->delete_snp( SNPs.at(ii));
	}
	outstream << endl;
	
}

/**
 * For a given rule, print it using the map data and add to delete list.
 */
void Bagging::print_and_process(int bags, AD_Rule &a){
	
	vector<long> rules;
	a.used(rules);
	unsigned int i;
	outstream << "In " << bags << " out of " << this->bag_param->get_number_of_bags() << " bags:  ";
	for(i=0;i<rules.size()-1;i++){ // NOTE: ignore last one.
		outstream << data->snp_name(rules.at(i)) << " ";
	}
	outstream << data->snp_name(rules.at(i)) << endl;
	this->data->getDataObject()->delete_snp( rules.at(i));
}


//////////////////////////////////////////////////////////////////////
// Classification.
/////////////////////////////////////////////////////////////////////
void Bagging::classify(int a[4]){} // Run test on the actual data.
void Bagging::classify(int a[4], SnpData *){} // Run test on passed data.
int Bagging::classify(vector<short> &){return 0;}// Run test on single individual.
