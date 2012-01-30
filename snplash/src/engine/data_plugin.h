/*
 *      data_plugin.h
 *      
 *      Copyright 2010 Richard T. Guy <guyrt@guyrt-lappy>
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

/**
 * Act as an intermediary between the data storage engine and data
 * accessors (generally, engines.) 
 * 
 * Utilize an optional permutation vector for redirection. 
 */

#ifndef DATA_PLUGIN_H
#define DATA_PLUGIN_H

#include "snp_data.hh"
using namespace std;

class DataAccess {
	
	public:
	
		DataAccess();
		void init(SnpData *);
		~DataAccess();
		
		/* Set up the process. */
		void setRedirectBootstrap();
		
		/* Return pointer to all SNPs for an individual. */
		vector<short> *get_data(int);
		/* Return phenotype value for individual. */
		double get_phenotype(int);
		/* Return number of phenotypes (so number of people) */
		int pheno_size();
		/* Return number of SNPs */
		int geno_size();
		/* Return covariates for an individual. */
		vector<double> *get_covariates(int);
		/* Return the chromosome for a given SNP */
		string get_chrom(int);
		/* Return several pieces of information about a SNP */
		void get_map_info(int, string &chr, string &name, int &position);
		/* Return the major and minor alleles */
		void get_allele_codes(int, char &maj, char &min, char &ref);
		/* Return the maximum size of any map */
		int max_map_size();
		/* Return a single person's ID */
		string get_person_ID(int);
		/* Return the name of a single SNP */
		string snp_name(long l){ return data->snp_name(l); };
		/* Return number of snps in a sample */
		long snp_size(unsigned int o){return data->snp_size(o);}

		
		int max_person_ID_size(){return data->maxPersonIDLength;}
		
		void make_categorical(int low, int hi);
		
		SnpData *getDataObject(){return data;} // This will need to probably go away...
		
	protected:
		SnpData *data;
		vector<unsigned int> redirect;
		bool uses_redirect;
		bool owns_data; // if it was created here, then kill it here.
	
};

#endif
