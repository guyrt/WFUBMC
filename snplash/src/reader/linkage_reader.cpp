#include "linkage_reader.h"

LinkageReader::LinkageReader(){

}

LinkageReader::~LinkageReader(){}

/**
 * This is the entry point to any reader class.  Will perform
 * any necessary work to actually capture data.
 *
 * @param data A SnpData object that will hold the data we are collecting
 * @param params A parameter object holding information about what to collect.
 */
void LinkageReader::process(SnpData *data, ParamReader *params){

	this->getPhenotype(data, params);

	#if DBG_PROGRESS
		cout << "Entering genotype." << endl;
	#endif

	this->getGenotype(data, params);

	#if DBG_PROGRESS
		cout << "Leaving genotype." << endl;
		cout << "Entering verify_snp_lengths" << endl;
	#endif

	data->verify_snp_lengths();

	#if DBG_PROGRESS
		cout << "Leaving verify_snp_lengths" << endl;
		cout << "Entering verify_data_size" << endl;
	#endif

	// don't need to verify this forever but for now, yes.
	// Take out if we implement map based readin.
	data->verify_data_size_match();

	#if DBG_PROGRESS
		cout << "Leaving verify_data_size" << endl;
	#endif

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
 * Retrieve and process all information from the genotype file.
 *
 * A line has the following format after skipped columns are removed:
 *
 * identifier snp1.1 snp1.2 snp2.1 snp2.2 snp3.1 snp3.2 ect.
 *
 * We have to store two symbols per locus so we know how to assign initial
 * case/con.  Choose arbitrary order.
 *
 * Future coders: don't be afraid to bust this 100 line method up if so inclined.
 *
 */
void LinkageReader::getGenotype(SnpData *data, ParamReader *params){

	bool warned = false; // If there are any warnings, we want to print the first one we see.
	bool highVerbosity = params->get_verbosity(); // This will let us quickly know if we should include extra checks.

	ifstream infile;
	char s1, s2;
	vector<string> line;
	string id; // Will hold first elt of the line.

	int line_count = 1;
	
	// Push all vectors onto stack. Could get better performance here by
	// reserving space.
	
	long l=0;
	map<string, int>::iterator it;
	for(it = order_in_file.begin(); it != order_in_file.end(); it++){
		if(l < (*it).second) l = (*it).second;
	}
	for(long i=0;i<=l;i++){
		vector<short> vec;
		data->snp_data.push_back(vec);
	}

	infile.open(params->get_linkage_geno_file().c_str(), ifstream::in);
	if(!infile.is_open()){
		cerr << "Unable to open file " << params->get_linkage_geno_file() << endl;
		exit(0);
	}

	unsigned long start_spot = 1;  // Change one for skip.  Are we sure?
	while( ! this->getLine(&infile, &line, params, true)){

		// Data verification:
		if((line.size()) % 2 != 1){
			cerr << "Uneven number of columns, line " << line_count << endl << "Aborting." << endl;
			exit(0);
		}

		id = line.at(0);

		// Code the char sets.
		if(line_count == 1){
			for(unsigned long i=start_spot;i < line.size(); i+=2){

				s1 = *(line.at(i).c_str());
				s2 = *(line.at(i+1).c_str());

				// EDIT 3-1-2010
				if(s1 == '0' || s2 == '0'){
					data->character_list.push_back(' ');
					data->character_list.push_back(' ');
				}else if(s1 == s2){
					data->character_list.push_back(line.at(i).at(0));
					data->character_list.push_back(' ');
				}else{
					data->character_list.push_back(line.at(i).at(0));
					data->character_list.push_back(line.at(i+1).at(0));
				}

			}
		}else if (line_count == 2){
			// Now, go through and reserve space for each vector.
			size_t number_of_snps = data->character_list.size() / 2;
			
			for(size_t ii=0;ii<data->snp_data.size();ii++){
				data->snp_data.at(ii).reserve(number_of_snps);
			}
		}

		if(order_in_file.count(id) == 0){
			if(params->get_verbosity() > 1) cout << "Individual " << id << " line " << line_count << " unused." << endl;
			line_count++;
			continue;
		}else{
			data->setIndividualName(static_cast<unsigned int>(order_in_file[id]), id);
		}

		// Pull line apart and put it in data.
		for(unsigned long i=start_spot;i < data->character_list.size(); i+=2){

			// THIS IS A BIT OF A HACK
			// I am checking here for windowing and omitting the element if it is outside
			// the window.
			if(!params->in_window( (i-start_spot)/2+1 )){
				continue;
			}

			// Get elements
			s1 = *(line.at(i).c_str());
			s2 = *(line.at(i+1).c_str());

			if(highVerbosity){
				if(!warned && !((s1 >= '0' && s1 <= '9' ) || s1 == 'A'|| s1 == 'C'|| s1 == 'T'|| s1 == 'G'|| s1 == 'a'|| s1 == 'c'|| s1 == 'g'|| s1 == 't'|| s1 == 'B'|| s1 == 'b'|| s1 == '.')){
					warned = true;
					cerr << "SNP: " << (i-start_spot)/2 + 1 << "  Locus 1" << endl << "Line " << line_count << " had coding '" << s1 << " " << s2 <<"'.  If unexpected, please check data." << endl;
					cerr << "Further warnings suppressed." << endl;
				}
				if(!warned && !((s2 >= '0' && s2 <= '9') || s2 == 'A'|| s2 == 'C'|| s2 == 'T'|| s2 == 'G'|| s2 == 'a'|| s2 == 'c'|| s2 == 'g'|| s2 == 't'|| s2 == 'B'|| s2 == 'b'|| s2 == '.')){
					warned = true;
					cerr << "SNP: " << (i-start_spot)/2 + 1 << "  Locus 2" << endl << "Line " << line_count << " had coding '" << s1 << " " << s2 <<"'.  If unexpected, please check data." << endl;
					cerr << "Further warnings suppressed." << endl;
				}
			}


			if(s1 == '0' || s2 == '0' || s1 == '.' || s2 == '.'){
				// A missing in either means push missing onto stack.
				data->snp_data.at(order_in_file[id]).push_back(0);
			}else{
				// First, check if we have stored ' ' , ' ', which means we've always seen blanks here.

				if(' ' == data->character_list.at(i-start_spot) && ' ' == data->character_list.at((i-start_spot)+1)){

					if(s1 == s2){
						data->character_list.at(i-start_spot) = s1;

					}else{
						data->character_list.at(i-start_spot) = s1;
						data->character_list.at(i-start_spot+1) = s2;
					}

				}
				else{
					// Have we seen these elements before?
					// Check against both spots.  If not match first, check second is available
					// or matching.  Note that zeros (missing) are skipped.
					if(s1 != data->character_list.at((i-start_spot))){
						if(data->character_list.at((i-start_spot)+1) != s1){
							if(data->character_list.at((i-start_spot)+1) != ' '){
								if(!warned){
									cerr << "Message: " << data->character_list.at((i-start_spot)+1)  << " " << data->character_list.at((i-start_spot)) << endl;
									cerr << "Have: " << s1 << " " << s2 << endl;
									cerr << "Non-biallelic SNP: " << (i-start_spot)/2 + 1 << "  Locus 1" << endl << "Line " << line_count << endl;
									cerr << "Individuals that do not have characters " << data->character_list.at((i-start_spot)+1)  << " or " << data->character_list.at((i-start_spot)) << " will be treated as missing." << endl;
								}
								warned = true;
							}else{
								data->character_list.at((i-start_spot)+1) = s1;
							}
						}
					}
					if(s2 != data->character_list.at((i-start_spot))){
						if(data->character_list.at((i-start_spot)+1) != s2){
							if(data->character_list.at((i-start_spot)+1) != ' '){
								if(!warned){
									cerr << "Message: " << data->character_list.at((i-start_spot)+1)  << " " << data->character_list.at((i-start_spot)) << endl;
									cerr << "Have: " << s1 << " " << s2 << endl;
									cerr << "Non-biallelic SNP: " << (i-start_spot)/2 + 1 << "  Locus 1" << endl << "Line " << line_count << endl;
									cerr << "Individuals that do not have characters " << data->character_list.at((i-start_spot)+1)  << " or " << data->character_list.at((i-start_spot)) << " will be treated as missing." << endl;
								}
								warned = true;
							}else{
								data->character_list.at((i-start_spot)+1) = s2;
							}
						}
					}
				}

				// Push this element onto the end of the last genotype.
				// Use coding described in snp_data.h.
				short push_val = 1;
				char d1, d2;
				d1 = data->character_list.at(i-start_spot);
				d2 = data->character_list.at(i-start_spot+1);

				if(s1 == d1 && s2 == d1){
					push_val = 1;
				}else if(s1 == d1 && s2 == d2){
					push_val = 2;
				}else if(s1 == d2 && s2 == d1){
					push_val = 3;
				}else if(s1 == d2 && s2 == d2){
					push_val = 4;
				}else{
					push_val = 0;
				}

				data->snp_data.at( order_in_file[id] ).push_back( push_val );
			}

		}

		line_count++;
	}

	#if DBG_PROGRESS
		cout << "First readthrough done. " << endl;
	#endif

	if(line_count != static_cast<int>(order_in_file.size())+1){ // note: line_count starts at 1 for the user.
		cout << "Warning: Uneven number of individuals in phen and gen.  If this message was unexpected, please verify data integrity." << endl;
	}

	infile.close();

	// Fix char_list to window size.
	vector<char> t;
	for(unsigned int i=0; i < data->character_list.size(); i+=2){
		if(params->in_window(i/2+1)){
			t.push_back(data->character_list.at(i));
			t.push_back(data->character_list.at(i+1));
		}
	}
	data->character_list = t;
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
void LinkageReader::getPhenotype(SnpData *data, ParamReader *params){

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
 * Expects five columns:
 *
 */
void LinkageReader::getMapFile(SnpData *data, ParamReader *params){

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
			cerr << "chr name 0 pos <opt ref>" << endl << "example: 22 rs1000 0 32432 [A]" << endl;
			exit(0);
		}
		if( params->in_window(l)){
			if(line.size() == 4){
				data->push_map(line.at(0), line.at(1), atol(line.at(3).c_str()));
			}else{
				data->push_map(line.at(0), line.at(1), atol(line.at(3).c_str()), line.at(4).at(0));
			}
		}
		l++;
	}

	infile.close();
}

/**
 * Verify that all traits and covariates were found.
 */
void LinkageReader::verify_phen_header(int trait_pos, vector<int> *cov_pos, vector <string> *line, ParamReader *params){

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
