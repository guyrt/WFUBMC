/*
 *      data_plugin.cpp
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

#include "data_plugin.h"

using namespace std;

DataAccess::DataAccess(){
	uses_redirect = false;	
	data = NULL;
	owns_data = false;
}

DataAccess::~DataAccess(){
	if(owns_data) delete data;
}

void DataAccess::init(SnpData *d){
	delete data;
	if(d == NULL){
		data = new SnpData();
		owns_data = true;
	}else{
		data = d;
		owns_data = false;
	}
}
/* Get a redirect set using bootstrapping */
void DataAccess::setRedirectBootstrap(){
	uses_redirect = true;
	redirect.clear();
	for(unsigned int i=0;i < data->phenotypes.size();i++){
		redirect.push_back(static_cast<int>(data->random.get() * static_cast<double>(data->phenotypes.size())));
	}
	
}

/**
 * Return all SNPs for a single individual in form of a pointer.
 * If this DataAccess object uses redirection, then the method figures out
 * which individual to actually return.
 *
 * @param i  The individual whose data is returned
 * @return vector<short> *  A pointer to the vector of data for the individual queried. 
 */
vector<short> * DataAccess::get_data(int i){
	if(uses_redirect){
		return data->get_column(redirect.at(i));
	}else{
		return data->get_column(i);
	}
}

double DataAccess::get_phenotype(int i){
	if(uses_redirect){
		return data->phenotype_at(redirect.at(i));
	}else{
		return data->phenotype_at(i);
	}
}

/* Return number of phenotypes (so number of people) */
int DataAccess::pheno_size(){
	if(uses_redirect){
		return redirect.size();
	}else{
		return data->phenotypes.size();
	}
}

/**
 * Return number of SNPs in the data set
 * 
 * @return int Number of SNPs.
 */
int DataAccess::geno_size(){
	return data->snp_size(0);
}

/**
 * Return a person ID
 * @param int the person index.
 */
string DataAccess::get_person_ID(int i){
	return data->individualName[i];
}

/* Return covariates for an individual. */
vector<double> * DataAccess::get_covariates(int i){
	if(uses_redirect){
		return data->covariate_column(redirect.at(i));
	}else{
		return data->covariate_column(i);
	}
}

/* Return chromosome for a given SNP. */
string DataAccess::get_chrom(int i){
	return data->snp_chr(i);
}

/*
 * Return the physical information about a given marker.
 */
void DataAccess::get_map_info(int i, string &chr, string &name, int &position){
	SnpData::MapData m = data->map.at(i);
	chr = m.chr;
	name = m.name;
	position = m.pos;
}

/* Return major and minor alleles */
void DataAccess::get_allele_codes(int i, char &maj, char &min, char &ref){
	data->fill_allele_codes(i,maj,min,ref);
}

/* Return maximum size of the SNP names */
int DataAccess::max_map_size(){
	return data->maxMapLength;
}

void DataAccess::make_categorical(int low, int hi){
	data->categorical_phenotype(low, hi);
}
