//      intertwolog.cpp
//      
//      Copyright 2010 Richard T. Guy <richardtguy84@gmail.com>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.

#include "intertwolog.hh"

#include "../utils/lr.hh"
#include "../../logger/log.hh"


InterTwoLog::InterTwoLog(){
	this->param_reader = ParamReader::Instance();
	this->data = new DataAccess;
	this->data->init(NULL);
	this->itl_param = new EngineParamReader;
	haveOwner = false;
}

InterTwoLog::~InterTwoLog(){
	if(!haveOwner) delete_my_innards();
}

void InterTwoLog::delete_my_innards(){
	
	this->param_reader = NULL;
	delete this->data;
	data = NULL;
	delete this->itl_param;
	itl_param = NULL;
}

/*
 * Read engine specific parameters and process the data.
 */
void InterTwoLog::init(){
	itl_param->read_parameters(param_reader->get_engine_specific_params());

	initializeReader(); // defined in engine.h
	
	if(param_reader->get_linkage_map_file().compare("none") == 0){
		cerr << "Intertwolog requires that you enter a map file.  Aborting." << endl;
		//exit(0);
	}

	reader->process(data->getDataObject(), param_reader);
	
	numInitSNPs = data->geno_size();
	numInitPhen = data->pheno_size();
	
	Logger::Instance()->init(param_reader->get_out_file() + ".log");
	
	data->getDataObject()->remove_phenotype(numeric_limits<double>::max());
	data->getDataObject()->remove_covariate(numeric_limits<double>::max());
	data->getDataObject()->remove_phenotype(0);
	data->getDataObject()->prep_data(param_reader);
	try{
		data->make_categorical(1,2);
	}catch(DataException d){
		cerr << "Data Exception with message: " << endl << d.message << endl;
		exit(0);
	}catch(...){
		cerr << "Exception making categorical." << endl;
	}
	delete reader; // No longer needed.
	
	if(data->pheno_size() == 0){
		cerr << "Error: there are no individuals.  Program halting." << endl;
		exit(0);
	}
	
	numFinalPhen = data->pheno_size();
	
}

/*
 * Prep the output engine.
 */
void InterTwoLog::preProcess(){
	
	int numCase = 0;
	for(int i=0;i < data->pheno_size(); ++i)
		if(equal(data->get_phenotype(i), 2)) numCase++;

	stringstream ss;
	ss << "INDIVIDUALS READ FROM THE INPUT FILE:            " << numInitPhen << endl;
	ss << "INDIVIDUALS DELETED                              " << numInitPhen - numFinalPhen << endl;
	ss << "INDIVIDUALS LEFT FROM THE INPUT FILE:            " << numFinalPhen << "  (" << numCase << " cases and " << numFinalPhen - numCase << " controls)" << endl;

	cout << ss.str();

	if(!out.init(param_reader, itl_param, data->max_map_size(), ss.str())){
		cerr << "Error opening output files.  Aborting." <<
		endl << "If you did not expect this error, and you are on a Linux/Unix machine, please verify that you have permission to open the file requested." << endl;
		exit(0);
	}
	
}

/*
 * Process the entire data set using a series of parallel OMP for loops.
 */
void InterTwoLog::process(){
	
	int sz = data->geno_size();
	
	int idx = 0;
	
	for(int i=0;i < sz; i++){
		
		#if RUN_IN_PARALLEL
		#pragma omp parallel shared(sz, idx) 
		{
			#pragma omp for schedule(static)
		#endif
		for(int j=i+1;j < sz; j++){
			
			#if INTERTWOLOG_TESTLOOP
			cout << "Running " << i << " " << j << endl;
			#endif
			
			InterTwoLogMeasures itlm;
			processPair(i,j,itlm);
			itlm.index1 = i+1;
			itlm.index2 = j+1;
			
			string tempS;
			int tempP;
			data->get_map_info(i, tempS, itlm.name1, tempP);
			data->get_map_info(j, tempS, itlm.name2, tempP);
			
			#if INTERTWOLOG_TESTLOOP
				cout << "Running " << i << " " << j << " printing." << endl;
			#endif
			out.printLine(itlm, idx+j-i-1);
			#if INTERTWOLOG_TESTLOOP
				cout << "Running " << i << " " << j << " done." << endl;
			#endif
		}
		#if RUN_IN_PARALLEL
		}
		#endif
		idx += sz-i-1;
	}
	out.close();
	
}

/**
 * Create data and run test for a single pair of SNPs given by i and j.
 * 
 * @param i First SNP
 * @param j Second SNP
 * @return itlm Structure filled with p, beta, se.
 */
void InterTwoLog::processPair(int i, int j, InterTwoLogMeasures &itlm){
	
	vector<double> snp1, snp2, snpInt, ones;
	vector<double> phen_vec;
	vector<vector<double> > cov;
	
	#if INTERTWOLOG_TESTLOOP
		cout << "Running " << i << " " << j << " start." << endl;
	#endif
	
	/* Prep covariate matrix */
	vector<double> *tmp = data->get_covariates(0);
	if(tmp != NULL){
		for(unsigned int ii=0;ii<tmp->size();ii++){
			vector<double> a;
			cov.push_back(a);
		}
	}
	
	#if INTERTWOLOG_TESTLOOP
		cout << "Running " << i << " " << j << " cov filled." << endl;
	#endif
	
	// First pass: Get individual SNP averages. 
	// Do not build the interaction.
	double mean1 = 0, mean2 = 0;
	
	for(int people = 0;people < data->pheno_size();++people){
		
		double ph = data->get_phenotype(people)-1;
		short s1, s2;
		s1 = data->get_data(people)->at(i);
		s2 = data->get_data(people)->at(j);
		if(s1 * s2 > 0){  // if neither is 0.
			ones.push_back(1.0);
			phen_vec.push_back(ph);
			if(s1 == 1){
				snp1.push_back(-1);
				mean1--;
			}else if(s1 == 2){
				snp1.push_back(0);
				
			}else if(s1 == 3){
				snp1.push_back(0);
				
			}else if(s1 == 4){
				snp1.push_back(1);
				mean1++;
			}
			if(s2 == 1){
				snp2.push_back(-1);
				mean2--;
			}else if(s2 == 2){
				snp2.push_back(0);
				
			}else if(s2 == 3){
				snp2.push_back(0);
				
			}else if(s2 == 4){
				snp2.push_back(1);
				mean2++;
			}
		}
	}
	mean1 /= snp1.size(); mean2 /= snp2.size();
	
	// Second pass: do mean adjusted interactions.
	for(int people = 0;people < data->pheno_size();++people){
		
		short s1, s2;
		s1 = data->get_data(people)->at(i);
		s2 = data->get_data(people)->at(j);
		if(s1 * s2 > 0){
			double interaction;
			if(s1 == 1){
				interaction = -1 - mean1;
			}else if(s1 == 2){
				interaction = -mean1;
			}else if(s1 == 3){
				interaction = -mean1;
			}else if(s1 == 4){
				interaction = 1 - mean1;
			}
			if(s2 == 1){
				interaction *= (-1 - mean2);
			}else if(s2 == 2){
				interaction *= (0 - mean2);
			}else if(s2 == 3){
				interaction *= (0 - mean2);
			}else if(s2 == 4){
				interaction *= (1 - mean2);
			}
			snpInt.push_back(interaction);
			
			vector<double> *t = data->get_covariates(people);
			for(unsigned int jj=0;jj<t->size();jj++){
				cov.at(jj).push_back(t->at(jj));
			}
		}
	}
		
	#if INTERTWOLOG_TESTLOOP
		cout << "Running " << i << " " << j << " data filled." << endl;
	#endif
	// Use the covariate matrix as the actual test matrix.
	cov.push_back(ones);
	cov.push_back(snp1);
	cov.push_back(snp2);
	cov.push_back(snpInt);
	
	vector<vector<double> > inv_infmatrix = vecops::getDblVec(cov.size(), cov.size());
	vector<double> betas;
	LRStats l;
	LogisticRegression lr(itl_param->getRegressionConditionNumberThreshold());

	int retry = 0;
	double startVal = 0;  // value to start betas with.

	while(retry < 3){
		try{
			betas = lr.newtonRaphson(cov, phen_vec, inv_infmatrix, startVal);
			l = lr.getSingleStats(betas, inv_infmatrix, betas.size()-1);
			
			itlm.beta = betas.at(betas.size() - 1);
			itlm.pVal = l.pVal;
			itlm.SE = l.invInf;
			
			#if INTERTWOLOG_TESTLOOP
				cout << "Running " << i << " " << j << " try " << retry << " worked." << endl;
			#endif
			
			return;
		}catch (ConditionNumberEx ex){
			// The condition number of the information matrix is large. Check for separation.
			int separableVariable = lr.dataIsSeparable(cov, phen_vec);
			
			string tempS;
			int tempP;
			string name1;
			string name2;
			data->get_map_info(i, tempS, name1, tempP);
			data->get_map_info(j, tempS, name2, tempP);
			
			if (separableVariable < 0){
				// Error: poor conditioning.
				stringstream ss;
				ss << "Intertwolog SNPs " << name1 << " and " << name2 << ": Poor conditioning in information matrix." << endl;
				Logger::Instance()->writeLine(ss.str());
			}else{
				// Error: the poor conditioning is due to separable variable.
				stringstream ss;
				ss << "Intertwolog SNPs " << name1 << " and " << name2 << ": ";
				if (static_cast<unsigned int>(separableVariable) == cov.size() - 1){
					ss << "Interaction term is separable variable. ";
				}else if (static_cast<unsigned int>(separableVariable) == cov.size() - 2){
					ss << "Separable by " << name2;
				}else if (static_cast<unsigned int>(separableVariable) == cov.size() - 3){
					ss << "Separable by " << name1;
				}else{
					ss << "Separable by a covariate.";
				}
				ss << endl;
				Logger::Instance()->writeLine(ss.str());
			}
			retry = 100;
		}catch(...){

			#if INTERTWOLOG_TESTLOOP
				cout << "Running " << i << " " << j << " try " << retry << " failed." << endl;
			#endif

			retry++;
			if(retry == 1) {startVal = 0.5;}
			else if(retry == 2){ startVal = -0.5;}
			else{
				stringstream ss;
				ss << "Error on snp set " << i << " " << j << ".  Run error in Logistic Regression" << endl;
				Logger::Instance()->writeLine(ss.str());
				itlm.beta = -1;
				itlm.pVal = 2.000;
				itlm.SE = -1.000;
			}
		}
	}
	
}


void InterTwoLog::enslave(EngineParamReader *){
	cerr << "Enslave not supported for INTERTWOLOG." << endl;
}

void InterTwoLog::test(){}
