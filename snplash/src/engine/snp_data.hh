/*
 *      snp_data.h
 *
 *      Copyright 2009 Richard T. Guy <richardtguy84@gmail.com>
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


#ifndef SNP_DATA_H
#define SNP_DATA_H

/**
 * Contains the primary data set for all SNP methods.
 *
 * The only things we might implement here (later) are methods
 * related to data cleaning, unless this is something that
 * we would want to do on a method specific type.
 *
 * In the interest of quick access, we are going to implement all
 * data in vectors, which are thread-safe for access.  Shorts for snps.
 *
 * 0: Missing
 * 1: 1 1
 * 2: 1 2
 * 3: 2 1
 * 4: 2 2
 */

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "../param/param_reader.h"
#include "utils/exceptions.h"
#include "utils/float_ops.hh"
#include "randwh.h"

class SnpData {


	friend class LinkageReader;
	friend class BinaryBedReader;
	friend class DataAccess;

	protected:
	struct MapData {
		long pos;
		string chr;
		string name;
		char refAllele;
		bool flipped; // true if ref is not minor.
	};

	public:
		SnpData();
		~SnpData();

		vector<char> character_list;  		// stores the characters that were used to represent SNPs.
									  		// 0,1 are snp1 and 2,3 are snp2, ect.
		
		// Input function(s)
		// Used to insert an element into the map array.
		void push_map(string, string, long, char ref = (char) 0);
		void setIndividualName(unsigned int i, string name);
		void cleanIndividualNames(); // Clean this up.
		
		/// Data cleaning functions
		void prep_data(ParamReader *);
		/* Call this to categoricalize the phenotype */
		void categorical_phenotype(int low, int hi);
									  
		void remove_haplotype();      // 3 -> 2 in data.

		int remove_phenotype(double); // remove all individuals with phenotype listed.
		int remove_covariate(double); // same, but for covariates.
		int remove_indiv_with_snp_value(int); // remove all individuals with a single SNP that has given value.

		int remove_missing_geno(); // remove if all individuals are missing.

		bool isUsable(int i){return !all_missing.at(i);}

		int numIndividuals(){return phenotypes.size();}

		/// Data manipulation functions
		bool delete_snp(long l); // Delete a SNP by inserting into buffer.
		void snp_flush(); // Actually perform delete.  Two part process: delete_snp(); flush(); flush calls perform_delete();

		bool delete_indiv(int i); // Delete an indiv by inserting into buffer.
		void indiv_flush(); // Actually perform the delete.

		void dump_data();

	protected:

		/* The following data access routines are allowed for friends only. */
		
		/// Access functions
		// Returns column using bootstrap in second param.
		vector<short>* get_column(unsigned int); 
		// Return a pointer to a column of covariates.
		vector<double>* covariate_column(int person);
		// Return a single phenotype.
		double phenotype_at(int i){	return phenotypes.at(i);}
		
		// Return the name of a SNP.
		string snp_name(long l);
		
		// Returns a single chromosome
		string snp_chr(int i){return map.at(i).chr;}
		// Returns the major and minor allele for an individual.
		void fill_allele_codes(int i, char &maj, char &min, char &ref);
		// Return number of snps in sample.
		long snp_size(unsigned int o){return snp_data.at(o).size();}
		
		vector<vector< short> > snp_data;  // Each inner vec stores an individual's genotypes.
		vector<vector< double > > covariance; 
		vector<double> phenotypes;			
		vector<MapData> map;
		vector<long> snp_delete_records; // Holds set of records that are to be deleted. (snps)
		vector<int> indiv_delete_records;
		vector<bool> all_missing; // Holds whether missing for all individuals or not.
		
		vector<string> individualName; // Hold original individual names.

		unsigned int maxMapLength; // holds maximum length of a SNP name.
		int maxPersonIDLength; // holds max length of person ID.

		void normalize();
		void recode_snp(long);
		void perform_delete_of_snp(long); // Actually performs a delete operation.
		void perform_delete_of_indiv(int); // Actually performs a delete operation.

		// Check fxns.
		void verify_snp_lengths();
		void verify_data_size_match();
		void verify_map_length();

		RandWH random;

};

#endif
