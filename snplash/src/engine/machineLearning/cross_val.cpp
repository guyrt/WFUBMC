/*
 *      cross_val.cpp
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
#include "cross_val.h"

CrossValidation::CrossValidation(){
	param_reader = ParamReader::Instance();
	data = new DataAccess;
	cross_param = new EngineParamReader;
	haveOwner = false;
	order_in_bag = 0;
}

/**
 * Open the data and read it into the data class. * 
 */
void CrossValidation::init(){
	
	cross_param->read_parameters(param_reader->get_engine_specific_params());

	initializeReader(); // defined in engine.h

	reader->process(data->getDataObject(), param_reader);

	data->getDataObject()->prep_data(param_reader);
	data->make_categorical(-1,1);
	data->getDataObject()->remove_haplotype();
	delete reader; // No longer needed.	
	
}
/**
 * Create all of the engines and set up the master-slave system.
 * Create data first.
 */
void CrossValidation::preProcess(){
	
	// Create data.
	int start_pos = 0;//data->create_cv_data(cross_param->get_number_of_bags(), this->order_in_bag, this->out_of_bag);
	
	if(cross_param->get_engine_type() == ParamReader::ADTREE){
		ADTree *e ;//= new ADTree(param_reader, data);
		for(int i=0; i < cross_param->get_number_of_bags(); i++){
			e = new ADTree(data);
			engines.push_back(e);
		}
		for(int i=0; i < cross_param->get_number_of_bags(); i++){
//			engines.at(i)->enslave(cross_param, i+start_pos);
		}
	}else if(cross_param->get_engine_type() == ParamReader::BAGGING){
		Bagging *e ;//= new ADTree(data);
		for(int i=0; i < cross_param->get_number_of_bags(); i++){
			e = new Bagging(data);
			engines.push_back(e);
		}
		for(int i=0; i < cross_param->get_number_of_bags(); i++){
	//		engines.at(i)->enslave(cross_param, i+start_pos);
		}
	}else{
		cerr << "Warning: no engines assigned." << endl;
	}
}
/**
 * Make this engine controllable by another.
 */
void CrossValidation::enslave(EngineParamReader *b){
	cross_param = b;
	haveOwner = true;
}


/**
 * Build a model on each of the cross-validation parts.
 * Run the model on the left out individuals.
 * Report individuals.
 */
void CrossValidation::process(){
	
	int num_bags = cross_param->get_number_of_bags();
	
	#pragma omp parallel
	{
		#pragma omp for schedule(static) nowait
		for(int i=0;i < num_bags;i++){
			this->engines.at(i)->preProcess();
			this->engines.at(i)->process();
		}
	} // end pragma.
	
	// Classify each.
	for(unsigned int i=0;i < data->pheno_size();++i){
		
		int oob = out_of_bag.at(i); // Get correct OOB tree.
		int score = engines.at(oob)->classify(*(data->get_data(i)));
		if(equal(data->get_phenotype(i) , -1)){
			if(score < 0){
				results[3]++;
			}else{results[2]++;}
		}else if(equal(data->get_phenotype(i) , 1)){
			results[ score >= 0 ? 0 : 1 ]++;
		}
	}	
}
void CrossValidation::test(){}

//////////////////////////////////////////////////////////////////////
// Classification.
/////////////////////////////////////////////////////////////////////
/**
 * Report classification, which was precomputed.
 * 
 */
void CrossValidation::classify(int a[4]){
	for(int i=0;i<4;++i)
		a[i] = results[i];
}
void CrossValidation::classify(int a[4], SnpData *){} // Run test on passed data.
int CrossValidation::classify(vector<short> &){return 0;} // Run test on single individual.
