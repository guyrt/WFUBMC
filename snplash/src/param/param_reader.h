/*
 *      param_reader.h
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

#ifndef PARAM_READER_H
#define PARAM_READER_H

#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <limits.h>
#include <sstream>

using namespace std;

/**
 *
 * General parameter reader.  Two sets of parameters exist for this program.
 *
 * Engine specific would have a -- start.  General will have a - start.  This
 * class handles all - parameters and passes -- on.
 *
 *
 * Available parameters:
 * Linkage file format (mutually exclusive with arff format.)
 * 	-geno \<file\>
 * 	-phen \<file\>
 * 	-map \<file\>
 * Arff format
 * 	-arff <file>
 * Engine specification:
 * 	-engine <adtree | bagging>
 *
 * Windowing parameters:
 * 	-beg \<int\> -end \<int\> paired
 *  -win <int> alone.  (two are exclusive)
 *
 * 	-skip <cols>
 *
 */

class ParamReader {

	public:

		static ParamReader* Instance();

		bool process_parameters(int argc, char * argv[]); // processes '-' params and
										// records rest in engine_specific_params.


		// For now, only LINKAGE and BINARY work.  Will add more.
		enum InputTypes { ARFF , LINKAGE, BINARY };
		enum EngineTypes { UNASSIGNED , BAGGING , ADTREE , SNPGWA , FORMAT, CROSSVAL , DPRIME, QSNPGWA, DANDELION, INTERTWOLOG};

		/// Skip related function
		bool use_this_column(int i);
		/// Windowing function
		bool in_window(int i);

		/// Keep all parameters private and use getter/setters
		string get_linkage_geno_file() {return linkage_geno_file;}
		string get_binary_geno_file() {return binary_geno_file;}
		
		string get_linkage_pheno_file() {return linkage_pheno_file;}
		
		string get_linkage_map_file() {return linkage_map_file;}
		
		string get_arff_file() {return arff_file;}
		
		string get_out_file() {return output_file;}
		
		string get_log_file() {return log_file;}
		
		int get_verbosity() {return verbosity;}
		string get_trait(){return trait;}
		
		
		bool get_ign(){return missing_ignore;}

		vector<string> get_covariates(){return covariates;}

		int get_begin() {return begin;}
		
		int get_end() {return end;}
		
		int get_win() {return window;}
		
		vector<string> * get_engine_specific_params() {return &engine_specific_params;}
		

		InputTypes get_input_type(){return file_type;}
		
		EngineTypes get_engine_types(){return engine;}

	private:

		/// Singleton stuff:
		ParamReader();  // Private so that it can  not be called
		ParamReader(ParamReader const&){};             // copy constructor is private
		ParamReader& operator=(ParamReader const& ){return *this;}  // assignment operator is private
		
		static ParamReader* m_readerInstance;

		void print_usage();

		/// Methods to verify input flags.
		bool resolve_input_type();
		bool check_necessary();
		bool check_conflicting();

		// Pull out a single string and put it in &place
		bool resolve_single_string(int argc, int &i, string &place, string token, char * argv[]);
		bool split_multiple_strings(int argc, int &i, vector<string> &vec, string token, char *argv[]);

		bool resolve_single_int(int argc, int &i, int &place, string ptoken, char *argv[]);
		bool resolve_engine_params(int argc, int &i, string token, char *argv[]);

		vector<string> engine_specific_params; // Hold any parameters that start with '--'

		/// Input file types
		string binary_geno_file;
		string linkage_geno_file;
		string linkage_pheno_file;
		string linkage_map_file;   // init to "none"
		string arff_file;

		/// Output file types
		string output_file;
		string log_file;

		/// Output constant
		int verbosity; // Can be 1,2,3 with increasing levels at each.  Implementation engine specific.

		/// Phenotype information
		string trait;
		vector<string> covariates;
		
		/// Data information
		bool missing_ignore;

		/// Data localization parameters
		int begin, end;  // Start and end of the data we want to use.  (1 based)
		int window;		 // UNUSED.
		vector<int> skip_cols; // vector of skipped columns.  (1 based)

		/// Data localization helper fxns
		bool process_skip(string);

		// Engine type
		EngineTypes engine;
		// Type of expected input files.
		InputTypes file_type;
};

#endif
