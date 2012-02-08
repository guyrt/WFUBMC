/*
 *      snp_data.cpp
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
#include "snp_data.hh"

using namespace std;

SnpData::SnpData(){
	maxMapLength = 0;
	maxPersonIDLength = 8;
}

SnpData::~SnpData(){

}

/**
 * Remove anyone with phenotype listed.
 * @param d All individuals with given phenotype are removed.
 * @return int Number of individuals that are removed.
 */
int SnpData::remove_phenotype(double d){
	int deleted = 0;
	for(unsigned int i=0;i<phenotypes.size();i++){
		if(phenotypes.at(i) == d){
			if(delete_indiv(i)) deleted++;
		}

	}
	indiv_flush();

	return deleted;
}

/**
 * Remove anyone with any covariate matching these.
 *
 * @param d All individuals with given covariate are removed.
 */
int SnpData::remove_covariate(double d){
	int deleted = 0;
	for(unsigned int i=0; i < phenotypes.size(); i++){
		for(unsigned int j=0;j<covariance.at(i).size();j++){
			if(covariance.at(i).at(j) == d) {
				if(delete_indiv(i)) deleted++;
			}
		}
	}
	indiv_flush();
	return deleted;
}

int SnpData::remove_missing_geno(){

	int deleted = 0;
	for(unsigned long i=0;i<snp_data.at(0).size();i++){
		bool flag = true;
		for(unsigned long j=0;j<snp_data.size();j++){
			if(snp_data.at(j).at(i) != 0) {
				flag = false;
				break;
			}
		}
		if(flag){
			// all missing.  erase.
			all_missing.push_back(true);
			deleted++;
		}else{
			all_missing.push_back(false);
		}
	}
	cout << "Omitting " << deleted << " SNPs that are missing data on all individuals." << endl;
	return deleted;

}

/*
 * Remove any individual that has the given value anywhere.
 */
int SnpData::remove_indiv_with_snp_value(int d){

	int deleted = 0;
	for(unsigned long i=0;i<snp_data.size();i++){
		for(unsigned long j=0;j<snp_data.at(i).size();j++){
			if(snp_data.at(i).at(j) == d) {
				if(delete_indiv(i)) deleted++;
				break;
			}
		}
	}
	indiv_flush();
	return deleted;
}

/**
 * Perform initial analysis of data.  Jobs
 * 	1) Code data as major allele and minor allele (may have to switch variables)
 *
 */
void SnpData::prep_data(ParamReader * params){

	remove_missing_geno();
	if (params->get_input_type() == ParamReader::LINKAGE)
		normalize();

}

/**
 *
 * Check locus counts and flip if no minor allele present.
 * Flip score if necessary.  Flip should be 1<->4 ; 3<->2
 *
 * Also store the minor allele.
 */
void SnpData::normalize(){

	vector<long> locus_counts;
	for(unsigned long l = 0; l < 2*snp_data.at(0).size(); l++)
		locus_counts.push_back(0);


	//#pragma omp parallel
	//{
		long stop_crit = static_cast<long>(snp_data.size());
		short s;
		//#pragma omp for schedule(static)
		for(long i = 0; i < stop_crit; i++){
			for(unsigned long j = 0 ; j < snp_data.at(i).size() ; j++){

				// Build based on count.
				s = snp_data.at(i).at(j);
				switch (s){
					case 1 :
						locus_counts.at(j*2)+=2;
					break;
					case 2 :
						locus_counts.at(j*2)++;
						locus_counts.at(j*2+1)++;
					break;
					case 3 :
						locus_counts.at(j*2)++;
						locus_counts.at(j*2+1)++;
					break;
					case 4 :
						locus_counts.at(j*2+1)+=2;
					break;
					default:

					break;
				}
			}
		}
	//}

	// Now correct those SNPs that have a reference allele.
	// If ref allele preassigned, then we want to check that they
	// are in the correct order and flip if not.
	// If the ref allele was not assigned, then we want to get
	// the locus counts.
	for(unsigned long l = 0; l < snp_data.at(0).size(); l++){

		if(map.at(l).refAllele != (char) 0){

			// If the first individual's first location is equal to the ref allele,
			// then keep same.  Otherwise, force swap.
			if(character_list.at(l*2) == map.at(l).refAllele){
				locus_counts.at(l*2) = 0;
				locus_counts.at(l*2+1) = 1;
				map.at(l).flipped = true;
			}else if (character_list.at(l*2+1) == map.at(l).refAllele){
				locus_counts.at(l*2) = 1;
				locus_counts.at(l*2+1) = 0;
			}else{
				cerr << "Error in map file, SNP " << l << ": reference allele " << map.at(l).refAllele << " is not found in geno file." << endl;
				cerr << "Alleles in geno file are " << character_list.at(l*2) << " " << character_list.at(l*2+1) << endl;
			}
		}
	}

	stop_crit = static_cast<long>(locus_counts.size());

	for(long i = 0; i < stop_crit; i+=2){
		if(locus_counts.at(i) < locus_counts.at(i+1)){
			// flip so the minor allele is second.
			char c = character_list.at(i);
			character_list.at(i) = character_list.at(i+1);
			character_list.at(i+1) = c;
			// recode
			recode_snp(i/2);
		}
	}

	locus_counts.clear();
}

/**
 * 1<->4 ; 2<->3
 */
void SnpData::recode_snp(long l){
	short s;
	for(unsigned long i = 0; i < snp_data.size(); i++){
		s = snp_data.at(i).at(l);
		switch (s) {
			case 1:
				snp_data.at(i).at(l) = 4;
			break;
			case 2:
				snp_data.at(i).at(l) = 3;
			break;
			case 3:
				snp_data.at(i).at(l) = 2;
			break;
			case 4:
				snp_data.at(i).at(l) = 1;
			break;
		}
	}
}

/**
 *
 * Verify that all of the individuals have the same number of snps.
 *
 * Will call exit(0) if this fails.
 *
 */
void SnpData::verify_snp_lengths(){

	unsigned int i = 0;
	unsigned long length;
	// Get base length;
	length = snp_data.at(i).size();
	while(length == 0){
		i++;
		length = snp_data.at(i).size();
	}
	for(i=1; i < snp_data.size(); i++){
		if(snp_data.at(i).size() == 0){
			// This was unused in phenotype.
			delete_indiv(i);
		}else{
			if(length != snp_data.at(i).size()){
				cerr << "Line " << i+1 << " did not have same number of snps as line 1.  " << endl;
				cerr << "Size was: " << snp_data.at(i).size() << endl;
				cerr << "Aborting" << endl;
				exit(0);
			}
		}
	}
	indiv_flush();
	cout << "Read " << length << " SNPs" << endl;
}

void SnpData::verify_map_length(){
	if(map.size() != snp_data.at(0).size()){
		cerr << "Map has incorrect number of SNPs.  Aborting." << endl;
		exit(0);
	}
}

/**
 * Verify that num phenotypes == num genotypes.  Prints a warning if mismatch,
 * but does not exit because the program is designed to match differently sized
 * files by individual.
 *
 * name: SnpData::verify_data_size_match()
 */
void SnpData::verify_data_size_match(){

	unsigned long length;
	// Get base length;
	length = snp_data.size();
	if(length != phenotypes.size()){
		cerr << "There are " << length << " genotypes but " << phenotypes.size() << " phenotypes." << endl;
	}
	cout << "Kept " << length << " genotype rows." << endl;

}

/**
 *
 * Insert a new map structure into the map array at end.
 *
 * Also updates the max length for all maps.
 *
 */
void SnpData::push_map(string c, string n, long p, char ref){
	MapData m;
	m.chr = c;
	m.name = n;
	m.pos = p;
	m.refAllele = ref;
	m.flipped = false;
	map.push_back(m);

	if(n.length() > maxMapLength) maxMapLength = n.length();
}

/**
 * Some of the vector elements are going to be empty.  Fill them in.
 */
void SnpData::cleanIndividualNames(){
	vector<string> newIDs;
	vector<string>::iterator it;
	for(it = individualName.begin(); it != individualName.end(); it++){
		if((*it).length() >= 1){
			newIDs.push_back(*it);
		}
	}
	individualName = newIDs;
}

/**
 * Insert a person ID into given location and update maxLength.
 *
 */
void SnpData::setIndividualName(unsigned int i, string name){
	while(individualName.size() < i + 1)
		individualName.push_back("");

	if(name.length() > static_cast<unsigned int>(maxPersonIDLength)){
		maxPersonIDLength = name.length();
	}
	individualName.at(i) = name;
}

/**
 * Check for catagorical phenotype then
 * place it in -1,1,0 order.
 *
 * -1,1 keep sorted order.
 *
 */
void SnpData::categorical_phenotype(int low, int hi){
	double p1 = -100;
	double p2 = -100;

	for(unsigned int i=0;i < phenotypes.size();i++){
		if(!equal(phenotypes.at(i), 0)){
			if(equal(-100, p1)){ // first pass through.
				p1 = phenotypes.at(i);
			}else if( equal(p1 ,phenotypes.at(i)) ){
				// do nothing
			}else if( equal(p2 ,-100)){ // found first spot of second one.
				p2 = phenotypes.at(i);
			}else if( equal(p2, phenotypes.at(i)) ){
				// do nothing
			}else{
				// To get here, must be a third type - error.
				DataException d;
				d.message = "Non categorical phenotype when categorical phenotype expected.";
				throw d;
			}
		}
	}

	// Check for same.
	if( equal(p1 ,-100) || equal(p2 ,-100) || equal(p1 , p2)){
		cerr << "Require two distinct phenotypes." << endl;
		exit(1);
	}


	// We are okay.  Check for reorder
	if(p2 < p1){
		double t = p1;
		p1 = p2;
		p2 = t;
	}


	for(unsigned int i=0;i < phenotypes.size() ; i++){
		if(equal(phenotypes.at(i) , p1)){
			phenotypes.at(i) = low;
		}else if(equal(phenotypes.at(i) , p2)){
			phenotypes.at(i) = hi;
		}
	}
}

/**
 * Recode all threes in the data as twos for algorithms that do not require haplo information. *
 */
void SnpData::remove_haplotype(){

	long outer = snp_data.size();
	long inner = snp_data.at(0).size();

	#pragma omp parallel
	{
		#pragma omp for
		for(long i=0;i < outer;i++){
			for(long j=0;j < inner;j++){
				if(snp_data.at(i).at(j) == 3){
					snp_data.at(i).at(j) = 2;
				}
			}
		}
	} // end pragma parallel

}

/**
 * Return the name of a given SNP.
 *
 * @param l SNP index in data set.
 * @return string Given name.
 */
string SnpData::snp_name(long l){
	if(l == -1){
		return "Base ";
	}
	else if(static_cast<long>(map.size()) > l){
		return map.at(l).name;
	}else{
		stringstream ss;
		ss << "out of range: " << l << " from " << map.size() << endl;
		DataException d;
		d.message = ss.str();
		throw d;
	}

}

/**
 * Retrieve an array pointer from the snp_data set.  The result
 * corresponds to an individual's genome.
 *
 * @param i corresponds to an individual.
 * @return vector<short>* pointer to SNPs for single individual.
 */
vector<short>* SnpData::get_column(unsigned int i){
	return &(snp_data.at(i));
}

vector<double>* SnpData::covariate_column(int person){
	if(static_cast<int>(covariance.size()) > person)
		return &(covariance.at(person));
	return NULL;
}

/**
 * Retrieve the major and minor alleles as passed into the program.
 *
 * @param i SNP of interest
 * @param maj Character to hold major allele
 * @param min Character to hold minor allele.
 * @param ref Character to hold the reference allele.
 */
void SnpData::fill_allele_codes(int i, char &maj, char &min, char &ref){

	if(map.at(i).refAllele == (char) 0 || !map.at(i).flipped){
		min = character_list.at(2*i+1);
		maj = character_list.at(2*i);
		ref = min;
	}else{
		maj = character_list.at(2*i+1);
		min = character_list.at(2*i);
		ref = maj;
	}
}

////////////////////////////////////////////////////////////////////////
// Removal...

/**
 * Remove a SNP from the sample data completely.
 *
 * Remove from the vector for each individual.
 * Remove from the map vector.
 *
 * Note: this algorithm is linear in number of elements after l.
 *
 * @param l vector position
 * @return bool true if completed.
 */
void SnpData::perform_delete_of_snp(long l){

	if(l > static_cast<long>(map.size())){
		return;
	}

	// Erase from map
	this->map.erase(map.begin() + l-1);

	// Erase from each individual.
	// TODO replace this with a copy of all nonerased individuals then a clear and replace.
	for(unsigned int i=0; i < snp_data.size(); ++i){
		this->snp_data.at(i).erase(this->snp_data.at(i).begin() + l);
	}
}

/**
 * Push l onto SNP delete vector if it is not already there.
 * For removal of SNPs
 */
bool SnpData::delete_snp(long l){
	for(unsigned int i=0;i < snp_delete_records.size();++i){
		if(snp_delete_records.at(i) == l){
			return false;
		}
	}
	snp_delete_records.push_back(l);
	return true;
}

/**
 * Perform removal and clear the delete vector.
 * For removal of SNPs
 */
void SnpData::snp_flush(){
	if(snp_delete_records.size() == 0){return;}

	vector<long>::reverse_iterator it;
	sort(snp_delete_records.begin(), snp_delete_records.end());
	for(it = snp_delete_records.rbegin(); it != snp_delete_records.rend(); ++it){
		//perform_delete_of_snp(*it);
	}

	// Remove elements from the map.
	vector<MapData> map_temp;
	vector<long>::iterator lit = snp_delete_records.begin();
	for(unsigned int j=0; j < map.size(); j++){
		if(*lit == j){
			lit++;
		}else{
			map_temp.push_back(map.at(j));
		}
	}
	map = map_temp;
	map_temp.clear();

	// Remove elements from each individual.
	for(unsigned int i=0; i < snp_data.size() ; i++){
		vector<short> temp;
		vector<long>::iterator lit_inner = snp_delete_records.begin();
		for(unsigned int j=0; j < snp_data.at(i).size(); j++){
			if(*lit_inner == j){
				lit_inner++;
			}else{
				temp.push_back(snp_data.at(i).at(j));
			}
		}
		snp_data.at(i) = temp;
	}

	// Clear up the characterlist.
	vector<char> new_list;
	vector<long>::iterator cit;
	cit = snp_delete_records.begin();
	for(long i=0;i < static_cast<long>(character_list.size());i ++){
		if(*cit == i / 2){
			i++; // skip ahead two.
			cit++;
		}else{
			new_list.push_back(character_list.at(i));
		}
	}
	character_list = new_list;

	snp_delete_records.clear();
}

// Delete an indiv by inserting into buffer.
bool SnpData::delete_indiv(int l){
	for(unsigned int i=0;i < indiv_delete_records.size();++i){
		if(indiv_delete_records.at(i) == l){
			return false;
		}
	}
	indiv_delete_records.push_back(l);
	return true;
}

// Actually perform the delete.
void SnpData::indiv_flush(){
	if(indiv_delete_records.size() == 0){return;}

	vector<int>::reverse_iterator it;
	sort(indiv_delete_records.begin(), indiv_delete_records.end());
	for(it = indiv_delete_records.rbegin(); it != indiv_delete_records.rend(); ++it){
		perform_delete_of_indiv(*it);
	}
	indiv_delete_records.clear();
}

// perform deletion.
void SnpData::perform_delete_of_indiv(int l){

	// Erase from phen.
	this->phenotypes.erase(this->phenotypes.begin() + l);

	// Erase from data
	this->snp_data.erase(snp_data.begin() + l);

	// Erase from cov.
	this->covariance.erase(covariance.begin() + l);
}

/*
 * As a test method, dump all the data.  Use format
 * phen cov data
 * Space delimited.
 */
void SnpData::dump_data(){

	for(unsigned int i=0;i < snp_data.size(); i++){

		cout << phenotypes.at(i) << " " ;
		for(unsigned int j=0;j < covariance.at(i).size();j++){
			cout << covariance.at(i).at(j) << " ";
		}
		for(unsigned int j=0;j < snp_data.at(i).size();j++){
			cout << snp_data.at(i).at(j) << " ";
		}
		cout << endl;
	}

}


