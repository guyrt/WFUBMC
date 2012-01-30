/*
 *      engine_param_reader.cpp
 *      
 *      Copyright 2009 Richard T. Guy <guyr7t@wfu.edu>
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
#include "engine_param_reader.h"

/**
 * 
 * @todo Make all of these defaults for speed.
 */
EngineParamReader::EngineParamReader(){
	
	number_of_crossvalidations = 10;
	number_of_bags = 10;
	threshold_of_bags = 3;
	number_of_nodes = 5;
	division_threshold = 0.005;
	bag_type = EngineParamReader::ALL;
	engine_type = ParamReader::ADTREE;
	
	dprime_smartpairs = false;
	dprime_fmt = 3;
	dprime_window = -1;
	
	snpgwa_dohaptest = true;
	
	output_haplo = output_geno = output_hwe = false;
	output_val = false;
	
	dandelion_pprob = false;
	haplo_thresh = -1;
	dandelion_window = -1;
	
	regression_condition_threshold = 1e12;
}

void EngineParamReader::read_parameters(vector<string> *params){
	
	string token;
	bool bad_start = false;
	
	for(unsigned int i = 0 ; i < params->size(); i++){
		token = params->at(i);
		if(token.compare("--nodes") == 0){
			i++;
			if(i >= params->size()){
				cerr << "Expected --node <int>" << endl;
				bad_start = true;
			}else{	// Get an integer
				token = params->at(i);
				int j = atoi(token.c_str());
				if(j <= 0){
					bad_start = true;
					cerr << "Node number must be positive.  Received " << token << endl;
				}else{
					number_of_nodes = j;
				}
			}
		}else if(token.compare("--bags") == 0){
			i++;
			if(i >= params->size()){
				cerr << "Expected --bags <int>" << endl;
				bad_start = true;
			}else{	// Get an integer
				token = params->at(i);
				int j = atoi(token.c_str());
				if(j <= 0){
					bad_start = true;
					cerr << "Bag number must be positive.  Received " << token << endl;
				}else{
					number_of_bags = j;
				}
			}
		}else if(token.compare("--bagthresh") == 0){
			i++;
			if(i >= params->size()){
				cerr << "Expected --bagthresh <int>" << endl;
				bad_start = true;
			}else{	// Get an integer
				token = params->at(i);
				int j = atoi(token.c_str());
				if(j <= 0){
					bad_start = true;
					cerr << "Bag threshold must be positive.  Received " << token << endl;
				}else{
					this->threshold_of_bags = j;
				}
			}
		}else if(token.compare("--method") == 0){
			i++;
			if(i >= params->size()){
				cerr << "Expected --method <int>" << endl;
				bad_start = true;
			}else{	// Get an integer
				token = params->at(i);
				int j = atoi(token.c_str());
				if(j < 1 || j > 3){
					bad_start = true;
					cerr << "Method number must be between 1 and 3.  Received " << token << endl;
				}else{
					if(j == 1){
						bag_type = EngineParamReader::ALL;
					}else if(j == 2){
						bag_type = EngineParamReader::LEAVES;
					}else{
						bag_type = EngineParamReader::SUBTREE;
					}
				}
			}
		}else if(token.compare("--engine") == 0){
			i++;
			if(i >= params->size()){
				cerr << "Expected --engine [adtree|bagging]" << endl;
				bad_start = true;
			}else{	// Get an integer
				token = params->at(i);
				if(token.compare("adtree") != 0){
					bad_start = true;
					cerr << "Only adtree engine available right now." << endl;
				}else{
					engine_type = ParamReader::ADTREE;
				}
			}
		}else if(token.compare("--dprime_smartpairs") == 0){
			dprime_smartpairs = true;
		}else if(token.compare("--snpgwa_nohap") == 0){
			snpgwa_dohaptest = false;
		}else if(token.compare("--val") == 0){
			output_val = true;
		}else if(token.compare("--geno_file") == 0){
			output_geno = true;
		}else if(token.compare("--haplo_file") == 0){
			output_haplo = true;
		}else if(token.compare("--hwe_file") == 0){
			output_haplo = true;
		}else if(token.compare("--condition_number") == 0){
			i++;
			if(i >= params->size()){
				cerr << "Expected --condition_number <number>" << endl;
				bad_start = true;
			}else{	// Get an integer
				token = params->at(i);
				int j = atoi(token.c_str());
				regression_condition_threshold = j;
			}
		}else if(token.compare("--dprime_window") == 0){
			i++;
			if(i >= params->size()){
				cerr << "Expected --dprime_window <number>" << endl;
				bad_start = true;
			}else{	// Get an integer
				token = params->at(i);
				int j = atoi(token.c_str());
				dprime_window = j;
			}
		}else if(token.compare("--dprime_fmt") == 0){
			i++;
			if(i >= params->size()){
				cerr << "Expected --dprime_fmt <number>" << endl;
				bad_start = true;
			}else{	// Get an integer
				token = params->at(i);
				int j = atoi(token.c_str());
				dprime_fmt = j;
			}
		}else if(token.compare("--dandelion_window") == 0){
			i++;
			if(i >= params->size()){
				cerr << "Expected --dandelion_window <number>" << endl;
				bad_start = true;
			}else{	// Get an integer
				token = params->at(i);
				int j = atoi(token.c_str());
				dandelion_window = j;
			}
		}else if(token.compare("--haplo_thresh") == 0){
			i++;
			if(i >= params->size()){
				cerr << "Expected --haplo_thresh <number>" << endl;
				bad_start = true;
			}else{	// Get an integer
				token = params->at(i);
				int j = atoi(token.c_str());
				haplo_thresh = j;
			}
		}else if(token.compare("--dandelion_pprob") == 0){
			dandelion_pprob = true;
		}else if(token.compare("--threshold") == 0){
			i++;
			if(i >= params->size()){
				cerr << "Expected --threshold number" << endl;
				bad_start = true;
			}else{	// Get an integer
				token = params->at(i);
				double j = atof(token.c_str());
				if(j < 0){
					bad_start = true;
					cerr << "Only positive thresholds allowed." << endl;
				}else{
					this->division_threshold = j;
				}
			}
		}else if(token.compare("--partition") == 0){
			i++;
			if(i >= params->size()){
				cerr << "Expected --partition <number>" << endl;
				bad_start = true;
			}else{	// Get an integer
				token = params->at(i);
				double j = atof(token.c_str());
				if(j < 0){
					bad_start = true;
					cerr << "Only positive partitions allowed." << endl;
				}else{
					this->division_threshold = j;
					if(j > 50){
						cerr << "Caution!  Using too many cross-validation partitions might yield unexpected results." << endl;
					}
				}
			}
		}else{
			cerr << "Unexpected argument: " << token << endl;
			bad_start = true;
		}
		
	}
		
	if(bad_start){exit(0);}
	
	// Verify that necessary parameters are here.
	bad_start = check_necessary();
	if(bad_start){exit(1);}
	bad_start = check_conflicting();
	if(bad_start){exit(1);}
}

bool EngineParamReader::check_necessary(){
	return false;
}

bool EngineParamReader::check_conflicting(){
	if(this->number_of_bags < this->threshold_of_bags){
		cerr << "Parameter error: Number of bags must be greater than or equal to threshold of bags." << endl;
		return true;
	}
	return false;
}

bool EngineParamReader::check_engine_param_conflict(){
	if(dprime_smartpairs && engine_type != ParamReader::DPRIME){
		cerr << "Parameter error: --dprime_smartpairs flag but engine is not DPRIME." << endl;
		return true;
	}
	if(dprime_fmt == 2 && dprime_window){
		cerr << "Parameter error: --dprime_window flag and --dprime_fmt 2 flag are in conflict.  Do not use this " ;
		cerr << "setting as the matrix output will not make sense." << endl;
		return true;
	}
	return false;
}
