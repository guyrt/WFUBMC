#include "binary_bed_reader.h"


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
		


}

/** getGenotype()
 *
 */
void BinaryBedReader::getBedFile(SnpData *data, ParamReader *params){

	bool warned = false; // If there are any warnings, we want to print the first one we see.
	char file_buffer[4];
	
	// This will let us know if we should include extra checks.
	bool highVerbosity = params->get_verbosity(); 

	FILE *in;
	unsigned long fileLen;
	
	// Open file (rb for read binary)
	in = fopen(params->get_linkage_geno_file().c_str(), "rb");
	if(in == NULL){printf("Error opening file\n.");}

	//Get file length
	fseek(in, 0, SEEK_END);
	fileLen=ftell(in);
	fseek(in, 0, SEEK_SET);	
	
	// Check the 3-byte header.
	fread(file_buffer, sizeof(char), 3, in);
	if (file_buffer[0] != 108 || file_buffer[1] != 27){
		cerr << "File is not binary .bed format" << params->get_linkage_geno_file() << endl;
		exit(0);
	}
	if(file_buffer[2] != 1){
		cerr << "File is not SNP-major format" << params->get_linkage_geno_file() << endl;
		exit(0);
	}

	// Read SNP 1 for all people, SNP2 for all people, etc (SNP-major mode)
	int personCntr = 0;
	short t; // hold the new bits.
	char c;
	int rowCntr=0;
	while (fread(file_buffer, sizeof(char), 1, in) > 0){
	
		c = file_buffer[0];
		t = c & 0x3;
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
		data->get_column(personCntr++)->push_back(binToCode(t));
		if (personCntr == data->numIndividuals()){
			personCntr = 0;
			rowCntr++;
			continue;
		}
		
		personCntr %= data->numIndividuals();
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

		
}

/**
 * Verify that all traits and covariates were found.
 */
void BinaryBedReader::verify_phen_header(int trait_pos, vector<int> *cov_pos, vector <string> *line, ParamReader *params){
	
}
