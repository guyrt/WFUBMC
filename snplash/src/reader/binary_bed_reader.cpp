#include "binary_bed_reader.h"
#include <stdio.h>
#include <stdlib.h>

#define DBG_PROGRESS 0

BinaryBedReader::BinaryBedReader(){
	
}

BinaryBedReader::~BinaryBedReader(){}

/** process()
 *
 * 
 * 
 * This is the entry point to any reader class.  
 *
 */
void BinaryBedReader::process(SnpData *data, ParamReader *params){
		
	this->getPhenotype(data, params);
	
	// Create a vector per person.
	long l = data->numIndividuals();
	for(long i=0;i<l;i++){
		vector<short> vec;
		data->snp_data.push_back(vec);
	}
	
#if DBG_PROGRESS
	cout << "Binary reader has " << l << " individuals." << endl;
#endif
	
	this->getBedFile(data, params);
	// don't need to verify this forever but for now, yes.
	// Take out if we implement map based readin.
	data->verify_data_size_match();


	if( 0 != ( params->get_linkage_map_file().compare("none"))){
		// If we have a map file, process it.
		this->getMapFile(data, params);
	}else{
		// Build a dummy map file.
		for(unsigned long i=0; i < data->snp_data.at(0).size(); i++){
			stringstream ss;
			ss << i+params->get_begin();
			data->push_map("0",ss.str(),0);
		}
	}

	data->verify_map_length();
	data->cleanIndividualNames();
}

/** getGenotype()
 *
 */
void BinaryBedReader::getBedFile(SnpData *data, ParamReader *params){

	char file_buffer[4];

	FILE *in;
	
	long minSnp = params->get_begin();
	long maxSnp = params->get_end();
	
	
	// Open file (rb for read binary)
	in = fopen(params->get_binary_geno_file().c_str(), "rb");
	if(in == NULL){
		cerr << "Error opening file " << params->get_linkage_geno_file() << endl;
		exit(0);
	}

	// Check the 3-byte header.
	int numBytes = fread(file_buffer, sizeof(char), 3, in);
	if (numBytes < 3){
		cerr << "Error reading file" << params->get_linkage_geno_file() << ": please check path and file size." << endl;
		exit(0);
	}
	if (file_buffer[0] != 108 || file_buffer[1] != 27){
		cerr << "File is not binary .bed format" << params->get_linkage_geno_file() << endl;
		exit(0);
	}
	if(file_buffer[2] != 1){
		cerr << "File is not SNP-major format" << params->get_linkage_geno_file() << endl;
		exit(0);
	}

#if DBG_PROGRESS
	cout << "Binary reader ready to start." << endl;
#endif

	// Read SNP 1 for all people, SNP2 for all people, etc (SNP-major mode)
	int personCntr = 0;
	short t; // hold the new bits.
	char c;
	long rowCntr=1;
	while (fread(file_buffer, sizeof(char), 1, in) > 0){
	
		c = file_buffer[0];
		t = c & 0x3;
		if (params->in_window(rowCntr))
			data->get_column(personCntr++)->push_back(binToCode(t));
		// These lines ensure that the rest of the byte is discarded if the
		// byte includes a row end. Check for all zeros for sanity.
		if (personCntr == data->numIndividuals()){
			
			if ((c & 0xFC) != 0){
				cerr << "Error in row " << rowCntr << " col 1 end of row expected but not found." << endl << "Aborting." << endl;
				exit(0);
			}
			personCntr = 0;
			rowCntr++;
			continue;
		}
		
		t = (c & 0xc) >> 2;
		if (params->in_window(rowCntr))
			data->get_column(personCntr++)->push_back(binToCode(t));
		if (personCntr == data->numIndividuals()){
			
			if ((c & 0xf0) != 0){
				cerr << "Error in row " << rowCntr << " col 2 end of row expected but not found." << endl << "Aborting." << endl;
				exit(0);
			}
			personCntr = 0;
			rowCntr++;
			continue;
		}
		
		
		t = (c & 0x30) >> 4;
		if (params->in_window(rowCntr))
			data->get_column(personCntr++)->push_back(binToCode(t));
		if (personCntr == data->numIndividuals()){
			if ((c & 0xCF) != 0){
				cerr << "Error in row " << rowCntr << " col 3 end of row expected but not found." << endl << "Aborting." << endl;
				exit(0);
			}
			rowCntr++;
			personCntr = 0;
			continue;
		}
		
		
		t = (c & 0xc0) >> 6;
		if (params->in_window(rowCntr))
			data->get_column(personCntr++)->push_back(binToCode(t));
		if (personCntr == data->numIndividuals()){
			personCntr = 0;
			rowCntr++;
			continue;

		}
	}

}

/** getPhenotype()
 *
 * Expects a dedicated phenotype file.  Read it.
 *
 * Header line has trait names.
 * Format:
 * header
 * personID <trait> <trait> ...
 *
 * Updates the classwide map variable.
 *
 */
void BinaryBedReader::getPhenotype(SnpData *data, ParamReader *params){

	ifstream infile;
	vector<string> line;

	unsigned int line_length = 0; // keeps track of length.
	int line_no = 0;			// Keeps track of individual order.
	int trait_pos = -1;			// This and next used to hold pull positions.
	vector<int> cov_pos;
	vector<int>::iterator it;

	vector<double> *cov_vec = new vector<double>;

	infile.open(params->get_linkage_pheno_file().c_str(), ifstream::in);
	if(! infile.is_open()){
		cerr << "Unable to open file " << params->get_linkage_pheno_file() << endl;
		exit(0);
	}

	// Get the header.
	if(this->getLine(&infile, &line, params, false)){
		cerr << "Only a header line in the phenotype file.  Did you mean to include any samples?" << endl;
		cerr << line.at(0) << endl;
	}

	line_length = line.size();

	if(params->get_trait().compare("none") == 0){
		trait_pos = 1;
	}

	// Identify trait and cov position.  Still run for cov even if trait
	// was defaulted.
	if(line.size() < 2){
		cerr << "Phenotype file has only one column (space delimited).  Aborting." << endl;
		exit(0);
	}

	for(unsigned int i=0; i < line.size();i++){
		if(line.at(i).compare(params->get_trait()) == 0){
			trait_pos = i;
		}
		for(unsigned int j=0; j < params->get_covariates().size() ; j++){
			if(line.at(i).compare(params->get_covariates().at(j)) == 0){
				cov_pos.push_back(i);
				break;
			}
		}
	}

	cout << "Param - cov: " << trait_pos << " - ";
	for(unsigned int i=0;i<cov_pos.size();i++)
		cout << cov_pos.at(i) << "  " ;
	cout << endl;

	// Verify that all covariates found and trait established.
	verify_phen_header(trait_pos, &cov_pos, &line, params);

	while(!this->getLine(&infile, &line, params, false)){
		if(line.size() != line_length){
			cerr << "Phenotype line " << line_no << " is incomplete.  Aborting." << endl;
			exit(0);
		}

		// Pull line apart and put it in data.
		// Push trait spot onto list.
		double t = atof(line.at(trait_pos).c_str());
		string s = line.at(trait_pos);
		if( s.compare(".") == 0 ){
			t = numeric_limits<double>::max();
		}
		data->phenotypes.push_back(t);
		// Add covariace
		data->covariance.push_back(*(cov_vec));
		for(it = cov_pos.begin(); it != cov_pos.end(); it++){
			t = atof(line.at(*(it)).c_str());
			s = line.at(*(it));
			if( s.compare(".") == 0 ){
				t = numeric_limits<double>::max();
			}
			data->covariance.at(data->covariance.size() - 1).push_back(t);
		}

		if(order_in_file.count(line.at(0)) > 0){
			cerr << "Phenotype " << line.at(0) << " was repeated on line " << line_no + 2<< " and " << order_in_file[line.at(0)] + 2 << " in the phenotype file.  Aborting." << endl;
			exit(0);
		}

		// Add position.
		order_in_file[line.at(0)] = line_no++;
	}
	infile.close();
	delete cov_vec;

}

/** getMapFile()
 *
 * Retrieve and process all information from the map file.
 *
 * Expects four columns plus two optional:
 * chr     snp             zero    pos     <a1>      <a2>
 * 1       SNP_A-8575125   0       554484  <T>       <C>
 *
 */
void BinaryBedReader::getMapFile(SnpData *data, ParamReader *params){

	ifstream infile;
	vector<string> line;

	infile.open(params->get_linkage_map_file().c_str(), ifstream::in);
	if(! infile.is_open()){
		cerr << "Unable to open file " << params->get_linkage_map_file() << endl;
		exit(0);
	}

	long l = 1;
	while(!this->getLine(&infile, &line, params, false)){
		// Pull line apart and put it in data.
		if(line.size() < 4){
			cerr << "Map line " << l << " has poor format.  Use format " << endl;
			cerr << "chr name 0 pos <opt ref> <opt minor>" << endl << "example: 22 rs1000 0 32432 [A] [T]" << endl;
			exit(0);
		}
		if( params->in_window(l)){
			if(line.size() == 4){
				data->push_map(line.at(0), line.at(1), atol(line.at(3).c_str()));
				data->character_list.push_back(' ');
				data->character_list.push_back(' ');
			}else if (line.size() == 5){
				data->push_map(line.at(0), line.at(1), atol(line.at(3).c_str()), line.at(4).at(0));
				// push min then maj allele onto character list.
				data->character_list.push_back(' ');
				data->character_list.push_back(line.at(4).at(0));
			}else if (line.size() > 5){
				data->push_map(line.at(0), line.at(1), atol(line.at(3).c_str()), line.at(4).at(0));
				// push min then maj allele onto character list.
				data->character_list.push_back(line.at(5).at(0));
				data->character_list.push_back(line.at(4).at(0));
			}
		}
		l++;
	}

	infile.close();
}

/**
 * Verify that all traits and covariates were found.
 */
void BinaryBedReader::verify_phen_header(int trait_pos, vector<int> *cov_pos, vector <string> *line, ParamReader *params){

	bool quit = false;

	if(trait_pos == -1){
		cerr << "Trait was not found." << endl;
		quit = true;
	}

	if(cov_pos->size() != params->get_covariates().size()){
		quit = true;
		cerr << "Missing covariate: " << endl;
		for(unsigned int i=0; i < params->get_covariates().size(); i++){
			string s = params->get_covariates().at(i);
			bool f = false;
			for(unsigned int j=0; j < line->size(); j++){
				if(s.compare(line->at(j)) == 0){f = true;}
			}
			if(!f){cerr << "   " << s << endl;}
		}
	}

	if(quit){exit(0);}
}
