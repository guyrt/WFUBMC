/*
 *      param_reader.cpp
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
#include "param_reader.h"

// Global static pointer used to ensure a single instance of the class.
ParamReader* ParamReader::m_readerInstance = NULL;  

ParamReader* ParamReader::Instance()
{
   if (!m_readerInstance)   // Only allow one instance of class to be generated.
      m_readerInstance = new ParamReader();

   return m_readerInstance;
}

ParamReader::ParamReader(){

	binary_geno_file = "none";
	linkage_map_file = "none";
	linkage_geno_file = "none";
	linkage_pheno_file = "none";

	arff_file = "none";

	output_file = "none";
	log_file = "none";
	verbosity = 1;

	trait = "none";

	begin = 1;
	end = INT_MAX;
	window = -1;  // What does window mean?

	engine = UNASSIGNED;
	file_type = ARFF;
	
	missing_ignore = false;
	
}

/**
 *
 * Read the input parameters and record them.  Call methods to check
 * for mistakes.
 *
 * See method verify() header for description of several bad cases.
 *
 */
bool ParamReader::process_parameters(int argc, char * argv[]){

	if(argc == 1){
		print_usage();
		return false;
	}

	string token;
	bool bad_start = false;

	// Loop through parameters.  Assign each
	for(int i=1; i < argc ; i++){

		token = argv[i];

		if(token.find("--") == 0 ){
			resolve_engine_params(argc, i, token, argv);
		}else if(token.compare("-bed") == 0){
			bad_start = bad_start || resolve_single_string(argc, i, this->binary_geno_file, token, argv);
		}else if(token.compare("-geno") == 0){
			bad_start = bad_start || resolve_single_string(argc, i, this->linkage_geno_file, token, argv);
		}else if(token.compare("-pheno") == 0 || token.compare("-phen") == 0){
			bad_start = bad_start || resolve_single_string(argc, i, this->linkage_pheno_file, token, argv);
		}else if(token.compare("-trait") == 0){
			bad_start = bad_start || resolve_single_string(argc, i, this->trait, token, argv);
		}else if(token.compare("-cov") == 0){
			// Pull apart the string that contains covariates.
			bad_start = bad_start || split_multiple_strings(argc, i, this->covariates, token, argv);
		}else if(token.compare("-map") == 0){
			bad_start = bad_start || resolve_single_string(argc, i, this->linkage_map_file, token, argv);
		}else if(token.compare("-arff") == 0){
			bad_start = bad_start || resolve_single_string(argc, i, this->arff_file, token, argv);
		}else if(token.compare("-out") == 0){
			bad_start = bad_start || resolve_single_string(argc, i, this->output_file, token, argv);
			this->log_file = this->output_file + ".log";
		}else if(token.compare("-beg") == 0){
			bad_start = bad_start || resolve_single_int(argc, i, this->begin, token, argv);
		}else if(token.compare("-end") == 0){
			bad_start = bad_start || resolve_single_int(argc, i, this->end, token, argv);
		}else if(token.compare("-win") == 0){
			bad_start = bad_start || resolve_single_int(argc, i, this->window, token, argv);
		}else if(token.compare("-v") == 0){
			bad_start = bad_start || resolve_single_int(argc, i, this->verbosity, token, argv);
		}else if(token.compare("-ign") == 0){
			missing_ignore = true;
		}else if(token.compare("-engine") == 0){
			i++;
			if(i >= argc){
				cerr << "Expected -engine [adtree | bagtree | snpgwa | dprime | dandelion | intertwolog]." << endl;
				bad_start = true;
			}else{
				token = argv[i];
				if(token.find('-') == 0){
					cerr << "Expected -engine [adtree | bagtree | snpgwa | dprime | dandelion | intertwolog] but got -engine " << token << endl;
					bad_start = true;
				}else{
					if(token.compare("adtree") == 0){
						engine = ADTREE;
					}else if(token.compare("bagging") == 0){
						engine = BAGGING;
					}else if(token.compare("cv") == 0){
						engine = CROSSVAL;
					}else if(token.compare("dprime") == 0){
						engine = DPRIME;
					}else if(token.compare("snpgwa") == 0){
						engine = SNPGWA;
					}else if(token.compare("qsnpgwa") == 0){
						engine = QSNPGWA;
					}else if(token.compare("dandelion") == 0){
						engine = DANDELION;
					}else if(token.compare("intertwolog") == 0){
						engine = INTERTWOLOG;
					}else{
						cerr << "The only acceptable engines are adtree, bagging, snpgwa, qsnpgwa, intertwolog, dandelion, and dprime." << endl;
						bad_start = true;
					}
				}
			}
		}else if(token.compare("-skipcols") == 0){
			i++;
			if(i >= argc){
				cerr << "Expected -skipcols <skipped columns>." << endl;
				bad_start = true;
			}else{
				token = argv[i];
				bool b = process_skip(token);
				if(!b){
					cerr << "Error in skip column formation: " << token << endl;
					bad_start = true;
				}
			}
		}else{
			cerr << "Token \"" << argv[i] << "\" unrecognized.  " << endl;
			bad_start = true;
		}
	}

	// Verify all parameters.
	if(bad_start){return false;}
	bad_start = this->resolve_input_type();
	if(bad_start){return false;}
	bad_start = this->check_necessary();
	if(bad_start){return false;}
	bad_start = this->check_conflicting();
	if(bad_start){return false;}
	return true;
}

/**
 * Resolves input type according to Enum InputTypes

 * binary_geno_file implies BINARY
 * linkage_geno_file implies LINKAGE
 * arff_file implies ARFF.
 *
 * Return: bool true if we found a problem.
 *
 */
bool ParamReader::resolve_input_type(){

	bool un_resolved = true;

	int num_hits = 0; // should be one at end.

	if (binary_geno_file.compare("none") != 0){
		this->file_type = BINARY;
		num_hits++;
		cout << "Reading linkage file." << endl;
	}
	if(arff_file.compare("none") != 0 ){
		this->file_type = ARFF;
		num_hits++;
		cout << "Arff file: " << arff_file.length() << endl;
	}
	if(linkage_geno_file.compare("none") != 0)
	{
		this->file_type = LINKAGE;
		num_hits++;
		cout << "Reading linkage file." << endl;
	}

	if(num_hits > 1){
		// There were too many hits (not enough covered later)
		cerr << "Conflicting file types discovered.  You must enter variables for one and only one input type." << endl;
	}else{
		un_resolved = false;
	}
	return un_resolved;
}

/**
 * Checks necessary options are there for each kind of input type.
 *
 * Assumes that input_type has been resolved.
 *
 * Special note on trait flag.  The default is last attribute
 * for arff files and SECOND column for linkage (after the head)
 *
 * Return type: true means that we identified a problem.
 */
bool ParamReader::check_necessary(){

	bool bad_start = false;

	// General to all input types.
	if(this->output_file.compare("none")==0){
		cerr << "You must include an output file." << endl;
		bad_start = true;
	}
	if(this->trait.compare("none") == 0){
		string str = "Trait flag not set.  Default behavior applies.";
		cerr << str << endl;
	}
	if(this->engine == UNASSIGNED){
		cerr << "No engine assigned.  Nothing to do!" << endl;
		bad_start = true;
	}

	// Specific to LINKAGE type
	if(this->file_type == LINKAGE){

		if(this->linkage_geno_file.compare("none") == 0){
			cerr << "Linkage format requires a genotype file." << endl;
			bad_start = true;
		}
		if(this->linkage_pheno_file.compare("none") == 0){
			cerr << "Linkage format requires a phenotype file." << endl;
			bad_start = true;
		}

	}else if(this->file_type == ARFF){
		cerr << "Arff file type not yet supported." << endl;
		bad_start = true;
	}

	if(begin < 1){
		cerr << "-beg cannot have value less than 1.  Input was " << begin << endl;
		bad_start = true;
	}

	return bad_start;
}

/**
 * Checks for lists of conflicting options.  Verifies both that paired
 * items are both present and that mutually exclusive not both there.
 *
 * Assumes that necessary items are all there and type is resolved.
 *
 * Return type: true means that we identified a problem.
 */
bool ParamReader::check_conflicting(){
	bool ret_type = false;

	return ret_type;
}

/**
 * Print usage string from a file.
 */
void ParamReader::print_usage(){
	
	stringstream ss;
	ss << "Snplash! by Richard T. Guy, Joshua D. Grab, Matt L. Stiegert, and Carl D. Langefeld." << endl;
	ss << "Command line usage:" << endl;
	ss << endl;
	ss << "For text input documents:" << endl;
	ss << "    -geno <geno file>      Input file for genetic data." << endl;
	ss << "    -phen <phenotype file> The phenotype input file." << endl;
	ss << "    -map <map file>        Contains information about the SNPs." << endl;
	ss << endl;
	ss << "For binary input documents:" << endl;
	ss << "    -bed <geno file>      Input file for genetic data in Plink format v0.991 or later." << endl;
	ss << "    -phen <phenotype file> The phenotype input file." << endl;
	ss << "    -map <map file>        Contains information about the SNPs." << endl;
	ss << endl;
	ss << endl;
	ss << "    -out <output file>     The primary output file.  Some engines may create several files by appending extra information." << endl;
	ss << endl;
	ss << "    -trait <string>    An element of the header in the phenotype file.  Optional, with the default being the second column in the phenotype file." << endl;
	ss << "    -cov <string,string,...,string>    An arbitrary number of covariates from the phenotype file.  Note that they should be comma separated."  << endl;
	ss << endl << "Data manipulation options" << endl;
	ss << "    -skipcols <int,int,int>   A list of columns to ignore in the input file.  Note that these columns are treated as if they do not exist, and this command is treated before begin, end.  Example: 1,3-6" << endl;
	ss << "    -beg <int>    The first SNP (1-based) to include " << endl;
	ss << "    -end <int>      The last SNP to include "<< endl;
	
	ss << endl;
	ss << "    -v <1,2, or 3>    The amount of printing to perform.  Not supported by all engines.  Primarily intended for use with machine learning engines." << endl;
	ss << "    -ign              If passed, any individuals with missing data in any SNP are excluded from the data set." << endl;
	
	ss << endl << endl << "    -engine <adtree | bagging | snpgwa | qsnpgwa | dprime | dandelion > " << endl;
	
	ss << endl << endl;
	ss << "Machine specific parameters:" << endl;
	ss << "GLOBAL" << endl;
	ss << endl;
	ss << "     --condition_number <number>    Maximum allowable condition number in logistic and linear regression." << endl;
	ss << "                               Condition number is measured in the 1-norm. Default is 1e12" << endl;
	ss << endl;
	ss << "ADTree" << endl;
	ss << "     --nodes <int>     The number of nodes to include in the tree" << endl;
	ss << endl;
	ss << "Bagging (also, see options for the engine being bagged)" << endl;
	ss << "     --bags <int>      The number of bootstraps to build" << endl;
	ss << endl;
	ss << "SNPGWA " << endl;
	ss << "     --val           If present, print the statistics value to <outfile>.statvals" << endl;
	ss << "     --geno_file	    If present, print extra genotypic files to <outfile>.geno[1,2,3]" << endl;
	ss << "     --haplo_file    If present, print extra haplotype files to <outfile>.haplo[1,2,3]" << endl;
	ss << "     --hwe_file      If present, print extra Hardy-Weinberg files to <outfile>.hwe[ | cntrl | case]" << endl;
	ss << "     --snpgwa_nohap  If present, haplotype tests are not computed.  This saves run time." << endl;
	ss << "     --haplo_thresh <int> If present, sets the threshold number of chromosomes on which a haplotype" << endl;
	ss << "                          must exist to be used for global haplotype association testing." << endl;
	ss << endl;
	ss << "DANDELION " << endl;
	ss << "     --dandelion_pprob  If present, create a file <outfile>.pprob and list each individual's personal probability of having each possible haplotype.  " << endl;
	ss << "     --dandelion_window <int> If present, perform dandelion on each set of <int> SNPs contiguously through the file." << endl;
	ss << "Send bug reports, including your computer's operating system, the full command line, and any additional information to dmcwilli@wfubmc.edu" << endl;
	cout << ss.str();
}

/**
 * Perform a simple parameter check and update on a single string.
 * 
 * @param argc Integer number of command line parameters
 * @param i Integer index in the parameter array.  Updated.
 * @param place String that will be updated.
 * @param token Token passed in
 * @param argv[] Pointer to char* of input parameters.
 */
bool ParamReader::resolve_single_string(int argc, int &i, string &place, string token, char * argv[]){
	
	i++;
	string next_token;
	if(i >= argc){
		cerr << "Expected another parameter after " << token << endl;
		return true;
	}else{
		next_token = argv[i];
		if(next_token.find('-') == 0){
			cerr << "Expected " << token << " <file> but got " << token << " " << next_token << endl;
			cerr << "Note: you can not start a file name with - " << endl;
			return true;
		}
		place = next_token;
	}
	return false;
}

/**
 * Split a string on comma
 * 
 * @param argc
 * @param i   The counter, updated.
 * @param &vec  Vector we put things in.
 * @param token  Token passed in.
 * @param *argv[] The parameter list.
 */
bool ParamReader::split_multiple_strings(int argc, int &i, vector<string> &vec, string token, char *argv[]){
	char c = ','; // Should parameterize this...
	
	i++;
	string next_token;
	if(i >= argc){
		cerr << "Expected string after " << token << "." << endl;
		return true;
	}else{
		next_token = argv[i];
		if(next_token.find('-') == 0){
			cerr << "Expected " << token << " <file> but got " << token << " " << next_token << endl;
			cerr << "Note: you can not start the cov list with - " << endl;
			return true;
		}
		// Split the string.
		string::size_type si = 0;
		string::size_type sj = next_token.find(c);
		while (sj != string::npos) {
			vec.push_back(next_token.substr(si, sj-si));
			si = ++sj;
			sj = next_token.find(c, sj);
		}
		vec.push_back(next_token.substr(si, next_token.length( )));
	}
	return false;
}

bool ParamReader::resolve_single_int(int argc, int &i, int &place, string ptoken, char *argv[]){
	i++;
	string token;
	if(i >= argc){
		cerr << "Expected integer after " << token << "." << endl;
		return true;
	}else{
		token = argv[i];
		place = atoi(token.c_str());
		if(this->window == 0){
			cerr << "Unable to resolve integer type after " << ptoken << ": " << place << endl;
			return true;
		}
	}
	token.clear();
	return false;
}

bool ParamReader::resolve_engine_params(int argc, int &i, string token, char *argv[]){
			
	if(token.compare("--nodes") == 0 || token.compare("--bags") == 0 || token.compare("--threshold")  == 0 ||
		token.compare("--bagthresh") == 0 || token.compare("--method") == 0 || token.compare("--partition") == 0 ||
		token.compare("--dprime_fmt") == 0 || token.compare("--dprime_window") == 0 || 
		token.compare("--haplo_thresh") == 0 || token.compare("--dandelion_window") == 0
		|| token.compare("--condition_number") == 0){
		engine_specific_params.push_back(token);
		i++;
		if(i >= argc){
			cerr << "Expected another token past " << token << endl;
			return true;
		}else{
			engine_specific_params.push_back(argv[i]);
		}
	}else if(token.compare("--dprime_smartpairs") == 0 || token.compare("--snpgwa_nohap") == 0
				|| token.compare("--val") == 0
				|| token.compare("--dandelion_pprob") == 0 || token.compare("--geno_file") == 0
				|| token.compare("--haplo_file") == 0 || token.compare("--hwe_file") == 0){
		engine_specific_params.push_back(argv[i]);
	}
	return false;
}

/**
 * Process the skip option, which would look something like
 * 1,2,3-8
 * 
 * Updates global variable skip_cols.
 * 
 * @param token representing the skipped columns in above format.
 * @return true if passed.
 */
bool ParamReader::process_skip(string token){
	
	/// Split on ,
	vector<string> toks;
	int start = 0;
	int end = 0;
	unsigned int i = 0;
	while(i < token.size()){
		
		if(token.at(i) == ','){
			if(start >= end)
				return false;
			toks.push_back(token.substr(start, end-start));
			end++;
			start = end; // skip over the ,
		}
		else if((token.at(i) >= '0' && token.at(i) <= '9' ) || token.at(i) == '-')
			end++;
		else
			return false;
		i++;
	}
	// push last one.
	
	if(start >= end){cout << start << " " << end << endl;
		return false;}
	toks.push_back(token.substr(start, end-start));
	end++;
	start = end; // skip over the ,

	/// Search for '-'.  
	for(unsigned int i=0;i < toks.size();i++){
		if(toks.at(i).find("-") != string::npos){
			string s = toks.at(i);
			
			unsigned int j = 0;
			while(j < s.size()){
				if(s.at(j) == '-')
					break;
				j++;
			}
			if(j == s.size() || j == 0) return false; // Misformed.
			int start = atoi(s.substr(0,j).c_str());
			int end = atoi(s.substr(j+1).c_str());
			if(end < start) return false; // Misformed: '6-5'
			for(int j = start;j <= end;j++)
				skip_cols.push_back(j);
		}else{
			skip_cols.push_back(atoi(toks.at(i).c_str()));
		}
	}
	return true;
}

/**
 * Return true if this is a column we are supposed to use.
 * @param i column label
 * @return bool true if use.
 */
bool ParamReader::use_this_column(int i){
		
	vector<int>::iterator it;
	for(it = skip_cols.begin(); it != skip_cols.end(); it++){
		if((*it) == i) return false;
	}
	
	return true;
}

/**
 * Return true if this is a SNP in the accepted window.
 * @param i A SNP index.
 */
bool ParamReader::in_window(int i){
	if( i >= this->begin && i <= this->end ){
		return true;
	}
	return false;
}
