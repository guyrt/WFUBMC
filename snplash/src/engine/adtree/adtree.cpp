/*
 *      adtree.cpp
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
#include "adtree.h"

/**
 * Set up the param_reader, data classes.
 * Leave the reader for init.
 */
ADTree::ADTree(){
	this->param_reader = ParamReader::Instance();
	this->data = new DataAccess;
	this->data->init(NULL);
	this->ad_param = new EngineParamReader;
	this->order_in_bag = 0;
	this->fudge = 0.0001;
	this->haveOwner = false;
}

ADTree::ADTree(DataAccess *d){
	this->param_reader = ParamReader::Instance();
	this->data = d;
	this->ad_param = new EngineParamReader;
	this->order_in_bag = 0;
	this->fudge = 0.0001;
	this->haveOwner = false;
}
ADTree::~ADTree(){
	if(!haveOwner) delete_my_innards();
}

void ADTree::delete_my_innards(){
	
	this->param_reader = NULL;
	delete this->data;
	data = NULL;
	delete this->ad_param;
	ad_param = NULL;
}

/** init()
 *
 * Set up the reader and input all data.
 * Call the data checker in the snp_data method.
 *
 */
void ADTree::init(){

	
	ad_param->read_parameters(param_reader->get_engine_specific_params());

	if(param_reader->get_input_type() == ParamReader::LINKAGE){
		reader = new LinkageReader;
	}else{
		cerr << "Other types don't work yet." << endl;
	}

	reader->process(data->getDataObject(), param_reader);
	
	data->getDataObject()->prep_data(param_reader);
	
	try{
		data->make_categorical(-1,1);
	}catch(DataException d){
		cerr << "Data Exception with message: " << endl << d.message << endl;
		exit(0);
	}catch(...){
		cerr << "Exception making categorical." << endl;
	}
	data->getDataObject()->remove_haplotype();
	//data->getDataObject()->pass_through_bootstrap();
	delete reader; // No longer needed.

}

/**
 * Delete current adparam reader and replace it with the referant.
 */
void ADTree::enslave(EngineParamReader *b){
	ad_param = b;
	//order_in_bag = order;
	haveOwner = true;
	data->setRedirectBootstrap(); // Set up a random redirect (bootstrap)
}

/** preProcess()
 *
 * Push the head node onto the stack.
 *
 */
void ADTree::preProcess(){
	
	this->weight_vec.clear();
	this->tree.reset();

	// Make weights.
	unsigned int size_of_pheno = this->data->pheno_size();
	for(unsigned int i=0;i < size_of_pheno ; i++){
		this->weight_vec.push_back(1.0 / size_of_pheno);
	}

	// Get simple weights
	double *simple_weights = new double[2];
	double weight_ratio = 0;
	weights(simple_weights);
	if(equal(simple_weights[1] , 0)){
		cerr << "All instances are in the same class!" << endl;
		cerr << "Bag: " << this->order_in_bag << endl;
		cerr << "Num instances: " << this->data->pheno_size() << endl;
		exit(0);
	}
	weight_ratio = .5 *log(simple_weights[0] / simple_weights[1]);

	// Update weights with correct values
	for(int i = 0; i < this->data->pheno_size() ; i++){
		weight_vec.at(i) *= (exp (-1.0*weight_ratio * data->get_phenotype(i)));
	}
	delete[] (simple_weights);

	Precondition *p = new Precondition;
	AD_Rule a(*p, p->last_condition(), weight_ratio,0);
	tree.push_node(a,0);
	tree.push_precondition(*p);
	delete p;
}

/** process()
 *
 * This calls the guts.  Start the algorithm and run it.
 * 
 * If there are covariates, run them first.
 *
 */
void ADTree::process(){

	if(param_reader->get_covariates().size() > 0){
		// First, build a tree on covariates.
		cout << "No covariates at this time." << endl;
	}else{
		processNoCov();
	}

	if(param_reader->get_verbosity() > 1){
		
		cout << print_tree();
		cout << "Along the way:" << endl;
		for(unsigned int i=0;i < percentTrue.size(); ++i){
			cout << percentTrue.at(i) << " ";
		}
		cout << endl;
	}


}

/**
 * Build the ADTree, assuming no covariates.
 * 
 */
void ADTree::processNoCov(){
	int temp_precon=0; long temp_attr=0; short temp_val=0;
	Condition::comparison temp_test = Condition::GE;
	double new_score_t, new_score_f;
	double *temp_weights = new double[5];

	AD_Rule *temp_r;
	Condition *c;

	//
	// READY TO RUN
	//
	while(keep_going()){
			
		double z = 0;
		// Get minimum score
		//z = minimize(&temp_precon, &temp_attr, &temp_test, &temp_val);
		z = minimize_ordinal(temp_precon, temp_attr, temp_test, temp_val);
		// Calc new scores
	
		weights(temp_weights, temp_precon, temp_attr, temp_test, temp_val);
		
		if(param_reader->get_verbosity() > 2){
			cout << "Num nodes: " << tree.node_size() << " score " << z << endl;
			cout << "Attribute: " << temp_attr << endl;
		}
		
		// Using 0.0005 for now, will change later to be 1/size
		new_score_t = .5 * log((temp_weights[0]+this->fudge) / (temp_weights[1]+this->fudge));
		new_score_f = .5 * log((temp_weights[2]+this->fudge) / (temp_weights[3]+this->fudge));
		
		// Make new rule and make two new preconditions.
		c = new Condition;
		c->attribute_index = temp_attr;
		c->genotype_conditional = temp_test;
		c->genotype_reference = temp_val;

		vector<Condition> condi_v = tree.precon_at(temp_precon)->conditions;
		Precondition *p1 = new Precondition;
		p1->conditions = condi_v;

		temp_r = new AD_Rule(*p1,*c, new_score_t, new_score_f);
		(*p1).conditions.push_back(*c);

		tree.push_node(*temp_r, (temp_precon+1)/2);

		tree.push_precondition(*p1);
		// Make a precondition with this set.
		c->inverse();
		condi_v = tree.precon_at(temp_precon)->conditions;
		Precondition *p2 = new Precondition;
		p2->conditions = condi_v;
		p2->conditions.push_back(*c);
		tree.push_precondition(*p2);
		delete c;
		delete p1;
		delete p2;
		delete temp_r;

		// Reweight data.
		for(unsigned int i=0; i < weight_vec.size(); i++){
			weight_vec.at(i) *= exp(-1.0 * tree.evaluate_on_last(data->get_data(i)) * data->get_phenotype(i));
		}
	}


	delete[] (temp_weights);	
}

/**
 * This is more extensible, but if data is known, use weights2 for speed.
 * Calc total weights of all individuals that satisfy each of the conditions.
 *
 * W_2 means are a 2 in phenotype.
 *
 * Passed the weight vector, as well as an index to the rule under the precondition, the attribute to test against
 * test type (see condition.h) and the value we are comparing to.  This could also be a Condition object but I
 * don't want to build one every time we use this method (waste of cycles...)
 *
 * Returned in ret_weights:
 * 		0 -> W_1 for p and c (W_+(p and c)) where W_+ is phenotype 1.
 * 		1 -> W_2 for p and c (W_-(p and c))
 * 		2 -> W_1 for p and not c (W_+(p and c))
 * 		3 -> W_2 for p and not c
 * 		4 -> W_1 or W_2 and not p
 * 
 * A better strategy would be to return the count for each value of the categorical variable.
 */
void ADTree::weights( double *ret_weights , int pre_cond, long attribute, Condition::comparison test, short val){

	bool val_at_p = false;
	bool val_at_c = false;

	// Initialize ret_weights to 0.
	for(int i=0;i < 5;i++){ret_weights[i] = 0;}

	Precondition *pre = tree.precon_at(pre_cond);
	for(int i=0; i < this->data->pheno_size() ; i++){
		if(pre->evaluate_truth(data->get_data(i))){
			val_at_p = true;
			switch (test){
				case Condition::GE :
					val_at_c = ((data->get_data(i))->at(attribute) >= val);	
				break;
				case Condition::GT :
					val_at_c = ((data->get_data(i))->at(attribute) > val);
				break;
				case Condition::LT :
					val_at_c = ((data->get_data(i))->at(attribute) < val);
				break;
				case Condition::LE :
					val_at_c = ((data->get_data(i))->at(attribute) <= val);
				break;
				case Condition::EQ :
					val_at_c = ((data->get_data(i))->at(attribute) == val);
				break;
				case Condition::NE :
					val_at_c = ((data->get_data(i))->at(attribute) != val);
				break;
			}
			
		}else{
			val_at_p = false;
			val_at_c = false;
		}
		if(val_at_p){

			if(equal(data->get_phenotype(i) , 1)){
				if(val_at_c){
					ret_weights[0] = ret_weights[0] + weight_vec.at(i);
				}else{
					ret_weights[2] = ret_weights[2] + weight_vec.at(i);
				}
			}else{
				if(val_at_c){
					ret_weights[1] = ret_weights[1] + weight_vec.at(i);
				}else{
					ret_weights[3] = ret_weights[3] + weight_vec.at(i);
				}
			}

		}else{
			ret_weights[4] = ret_weights[4] + weight_vec.at(i);
		}
	}
}


/**
 * Calc W_+ and W_- for each possible value of genotype and each value of 
 * p && c, p && !c.  Also, calculate for !p, which does not differentiate phenotypes.
 * 
 * @param Vector of weights to be filled.
 * @param long attribute on which we fill.
 */
void ADTree::weights2(double *ret_weights, int pre_cond, long attribute){
	
	bool val_at_p = false;
	bool val_at_c1 = false;
	bool val_at_c2 = false;
	vector<short>* vec;

	// Initialize ret_weights to 0.
	for(int i=0;i < 9;i++){ret_weights[i] = 0;}

	Precondition *pre = tree.precon_at(pre_cond);
	for(int i=0; i < this->data->pheno_size() ; i++){
		if(pre->evaluate_truth(data->get_data(i))){
			val_at_p = true;
			vec = data->get_data(i);
			val_at_c1 = (vec->at(attribute) >= 2);	
			val_at_c2 = (vec->at(attribute) >= 4); // so really just 4	
			
		}else{
			val_at_p = false;
			val_at_c2 = val_at_c1 = false;
		}
		
		if(val_at_p){

			if(equal(data->get_phenotype(i) , 1)){
				if(val_at_c1){
					ret_weights[0] = ret_weights[0] + weight_vec.at(i);
				}else{
					ret_weights[2] = ret_weights[2] + weight_vec.at(i);
				}
				if(val_at_c2){
					ret_weights[4] = ret_weights[4] + weight_vec.at(i);
				}else{
					ret_weights[6] = ret_weights[6] + weight_vec.at(i);
				}
			}else{
				if(val_at_c1){
					ret_weights[1] = ret_weights[1] + weight_vec.at(i);
				}else{
					ret_weights[3] = ret_weights[3] + weight_vec.at(i);
				}
				if(val_at_c2){
					ret_weights[5] = ret_weights[5] + weight_vec.at(i);
				}else{
					ret_weights[7] = ret_weights[7] + weight_vec.at(i);
				}
			}

		}else{
			ret_weights[8] = ret_weights[8] + weight_vec.at(i);
		}
	}
	
}

/**
 * Simply calculate weights for pos or neg class (1 or 2 in my coding)
 * Completely skip missings.
 *
 * Returns [phen==2,phen==1] for W_+, W_-
 */
void ADTree::weights(double *ret_weight){
	ret_weight[0] = 0; ret_weight[1] = 0;
	for(int i=0;i < data->pheno_size();i++){
		if(equal(data->get_phenotype(i) , -1)){
			ret_weight[1]++;
		}else if(equal(data->get_phenotype(i) , 1)){
			ret_weight[0]++;
		}
	}
}

/**
 * calculate stopping criteria for the main loop.  The simplest scheme is
 * to go based only on size.  A more nuanced scheme will rely on sufficient 
 * decrease or (even better) on a mix of the two.
 * @return boolean true if perform another round.
 *
 */
bool ADTree::keep_going(){
	
		if(param_reader->get_verbosity() > 2){
			cout << "Num nodes: " << tree.node_size() << " out of " << ad_param->get_number_of_nodes() << endl;
		}
		int *temp = new int[4];
		classify(temp);
		double d1 = static_cast<double>(temp[0]+temp[3]);
		double d2 = static_cast<double>(temp[0]+temp[1]+temp[2]+temp[3]);
		percentTrue.push_back(d1/d2);
		//percentTrue.push_back(temp[3]);
		delete[] temp;
		
		return (tree.node_size() < ad_param->get_number_of_nodes());
}

/**
 * Minimize the value Z(p,c) given in white paper.
 *
 * Note: this will need to be done in parallel soon...
 *
 * For each snp for each locus for each condition (roll this together) for each precondition
 * get score and compare.
 */
double ADTree::minimize(int *precon, long *attribute, Condition::comparison *test, short *val){

	double min_score = LONG_MAX;

	long final_attribute = -1;
	Condition::comparison final_test = Condition::GE;
	short final_val = -1;
	int final_precon = -1;

	int r;
	long i;
	long stop_size = static_cast<long> (data->snp_size(order_in_bag));
	// r is a precon and i is a snp we care about.

	for(r = 0; r < tree.precondition_size() ; r++){
		// Worth putting optimization right here?
		//#pragma omp parallel shared(final_attribute,final_precon,final_test,final_val)
		//{
			//#pragma omp for schedule(static)
			for(i=0; i < stop_size ; i++){
				double z;
				Condition::comparison c = Condition::GE;
				short s = 4;
				z = score_categorical(r,i,c,s);
				#pragma omp critical
				{
					if(z < min_score){
						final_attribute = i;
						final_precon = r;
						final_test = c;
						final_val = s;
						min_score = z;
					}
				}
				
			} // end inner for
	//	} // end pragma
	}

	// Assign the final values.
	*precon = final_precon;
	*attribute = final_attribute;
	*test = final_test;
	*val = final_val;

	return min_score;
}

/**
 * Perform a minimization using each possible score from the set {1,2,4}
 * 
 */
double ADTree::minimize_ordinal(int &precon, long &attribute, Condition::comparison &test, short &val){
	
	double min_score = LONG_MAX;
	long final_attribute = -1;
	
	int final_precon = -1;
	vector<double> weights;

	int r;
	long i;
	long stop_size = static_cast<long> (data->geno_size());
	
		for(r = 0; r < tree.precondition_size() ; r++){
		// Worth putting optimization right here?
		#pragma omp parallel shared(final_attribute,final_precon)
		{
			#pragma omp for schedule(static)
			for(i=0; i < stop_size ; i++){
				
				double z;
				
				short s = 4;
				vector<double> temp_w;
				
				z = score_ordinal(r,i,s,temp_w);
				
				#pragma omp critical
				{
					if(z < min_score){
						final_attribute = i;
						final_precon = r;
						min_score = z;
						weights = temp_w;
					}	
				}
				
			} // end inner for
		} // end pragma
	}
	
	/*
	 * Now, given the correct SNP to use, which is in final_attribute and final_precon,
	 * Calculate the three scores and group them.
	 */
	
	double w1 = .5 * log((weights.at(1) + this->fudge) / (weights.at(0) + this->fudge));
	double w2 = .5 * log((weights.at(3) + this->fudge) / (weights.at(2) + this->fudge));
	//double w4 = .5 * log((weights.at(5) + this->fudge) / (weights.at(4) + this->fudge));
	
	bool b1 = w1 < 0;
	bool b2 = w2 < 0;
	
	if(b1 == b2){
		val = 4;
	}else{
		val = 2;
	}
	precon = final_precon;
	attribute = final_attribute;
	test = Condition::GE;
	return min_score;
}

/**
 * Score for ordinal, which treats each type separately.
 * 
 * This is a customized version of the minimization in ADTrees that 
 * was created by Richard T. Guy and is under current investigation.
 */
double ADTree::score_ordinal(int pre_cond, long snp, short &s, vector<double> &temp_weights){
	
	double w1p = 0, w2p = 0, w4p = 0;
	double w1m = 0, w2m = 0, w4m = 0;
	double notp = 0;
	
	Precondition *pre = tree.precon_at(pre_cond);
	
	double qq = 0, pq = 0, pp = 0; // counts.
	
	for(int ii=0; ii < this->data->pheno_size() ; ii++){
		if(pre->evaluate_truth(data->get_data(ii))){
			double phen = data->get_phenotype(ii);

			switch (data->get_data(ii)->at(snp)){
				case 1:
					if(equal(phen,-1)){
						w1m += weight_vec.at(ii);
					}else{
						w1p += weight_vec.at(ii);
					}
					qq++;
					//qq += weight_vec.at(ii);
				break;
				case 2:
					if(equal(phen,-1)){
						w2m += weight_vec.at(ii);
					}else{
						w2p += weight_vec.at(ii);
					}
					pq++;			
					//pq += weight_vec.at(ii);
				break;
				case 4:
					if(equal(phen,-1)){
						w4m += weight_vec.at(ii);
					}else{
						w4p += weight_vec.at(ii);
					}
					pp++;
					//pp += weight_vec.at(ii);
				break;
				default:
				
				break;
			}			
		}else{
			notp += weight_vec.at(ii);
		}
	}
	
	temp_weights.clear();
	temp_weights.push_back(w1m);
	temp_weights.push_back(w1p);
	temp_weights.push_back(w2m);
	temp_weights.push_back(w2p);
	temp_weights.push_back(w4m);
	temp_weights.push_back(w4p);
	temp_weights.push_back(notp);
	
	//if(w1p > w2p || w2p > w4p) return 10000;
	double inner_sum = 0;
	if(qq>0) inner_sum += w1p*w1m * (pp+pq+qq) / qq;
	if(pq>0) inner_sum += w2p*w2m * (pp+pq+qq) / pq;
	if(pp>0) inner_sum += w4p*w4m * (pp+pq+qq) / pp;
	//cout << snp << " " << 2*sqrt(inner_sum )+notp << endl;
	return 2*sqrt(inner_sum )+notp;
}


/**
 * score_categorical
 * 
 * Perform scoring for a given precondition and variable.
 * 
 * @param int r: Precondition in use
 * @param long i: SNP in use
 * @param Condition::comparison &: Type of comparison returned.
 * @param short s: Optimal split point.
 */
double ADTree::score_categorical(int r, long i, Condition::comparison &c, short &s){
	
	double *ret_score = new double[9];
	double z1, z2, ret;
	weights2(ret_score,r,i); // run both tests.
	z2 = 2.0*(sqrt(ret_score[0]*ret_score[1]) + sqrt(ret_score[2]*ret_score[3])) + ret_score[8];
	z1 = 2.0*(sqrt(ret_score[4]*ret_score[5]) + sqrt(ret_score[6]*ret_score[7])) + ret_score[8];
	
	if(z1 < z2){
		ret = z1;
		c = Condition::GE;
		s = 4;
	}else{ // so default to z2.
		ret = z2;
		c = Condition::GE;
		s = 2;
	}
	
	delete[] ret_score;
	return ret;
}

///////////////////////////////////////////////////////////////////////////////////////
/*
 * Below this line is primarily reporting routines.
 * 
 */

/* Fill a hash with all rules in this tree.
 * 
 * name: ADTree::report
 * @param map<string, int> reference that is updated with information about the SNPs in this map.
 * @param map<string, AD_Rule> reference tying hashs to Rules.
 * @return none
 */
void ADTree::report(map<string, int> &hash , map<string, AD_Rule> &key){
	tree.report(hash, key);
}
void ADTree::report_leaves(map<string, int> &hash , map<string, AD_Rule> &key){
	tree.report_leaves(hash, key);
}

/**
 * 
 * name: ADTree::classify
 * @param int a[4] will contain TP,FP,FN,TN counts for the individuals in this test set.
 * Controls classified as controls would be in a[0]
 * P and N are defined as >=, < 0.
 */
void ADTree::classify(int a[4]){
	double score = 0;
	a[0] = a[1] = a[2] = a[3] = 0;
	for(int i=0;i < data->pheno_size();i++){
		score = tree.evaluate_instance(data->get_data(i));
		if(equal(data->get_phenotype(i) , -1)){
			if(score < 0){
				a[3]++;
			}else{a[2]++;}
		}else if(equal(data->get_phenotype(i) , 1)){
			a[ score >= 0 ? 0 : 1 ]++;
		}
	}
}
void ADTree::classify(int a[4], SnpData * data){
	
}

/**
 * Return the classification for this individual.
 * @return int 1 or -1 based on sign.
 */
int ADTree::classify(vector<short> &data){
	double score = tree.evaluate_instance(&data);
	return score > 0 ? 1 : -1;
}

/**
 * Print the tree to a file.
 * 
 * 
 */
void ADTree::print_tree_to_file(){
	
	string fileName = param_reader->get_out_file();
	ofstream outstream;
	outstream.open(fileName.c_str());
	
	outstream << print_tree();
	
	int *results = new int[4];
	classify(results);
	
	outstream << endl << endl << "Results: " << endl;
	outstream << "      Classified as" << endl;
	outstream << "Control   Cases" << endl;
	outstream << results[0] << "        " << results[1] << " | Controls" << endl;
	outstream << results[2] << "        " << results[3] << " | Cases " << endl;
	delete[] results;
	
	outstream.close();
	
}

/**
 *Perform a few tests, see internal comments.
 */
void ADTree::test(){

/*	for(int i=0;i < tree.precondition_size();i++){
		cout << i << " " << tree.precon_at(i)->to_string() << endl;
	}*/
	
	cout << "Size: " << tree.node_size() << endl;
	

}
