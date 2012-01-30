//      dandelion.cpp
//
//      Copyright 2010 Richard T. Guy <guyrt7@wfu.edu>
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

#include "dandelion.hh"

Dandelion::Dandelion(){
	this->param_reader = ParamReader::Instance();
	this->data = new DataAccess;
	this->data->init(NULL);
	this->ld_param = new EngineParamReader;
	haveOwner = false;
}

Dandelion::~Dandelion(){
	if(!haveOwner) delete_my_innards();
}

void Dandelion::delete_my_innards(){
	
	this->param_reader = NULL;
	delete this->data;
	data = NULL;
	delete this->ld_param;
	ld_param = NULL;
}

/**
 * Performs several functions, not all documented here.
 * 
 * Main tasks: get engine specific parameters and read and prep data.
 * 
 * Param checks:
 *   map file present (does not currently exit)
 * 
 */
void Dandelion::init(){
	
	ld_param->read_parameters(param_reader->get_engine_specific_params());

	


	if(param_reader->get_input_type() == ParamReader::LINKAGE){
		reader = new LinkageReader;
	}else{
		cerr << "Other types don't work yet." << endl;
	}

	if(param_reader->get_linkage_map_file().compare("none") == 0){
		cerr << "Dandelion requires that you enter a map file.  Aborting." << endl;
		//exit(0);
	}

	reader->process(data->getDataObject(), param_reader);

	numInitSNPs = data->geno_size();
	numInitPhen = data->pheno_size();

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

	// This removes any individuals that are missing a value at any of the SNPs.
	// We don't want to use this if we are windowing.
	//data->getDataObject()->remove_indiv_with_snp_value(0);

	if(data->pheno_size() == 0){
		cerr << "Error: there are no individuals.  Program halting." << endl;
		exit(0);
	}

	numFinalPhen = data->pheno_size();
}

/*
 * Set up output.
 * Get all of the maj, min alleles.
 * Send map info to output.
 */
void Dandelion::preProcess(){

	numCase = 0;
	for(int i=0;i < data->pheno_size(); ++i)
		if(equal(data->get_phenotype(i), 2)) numCase++;

	stringstream ss;
	ss << "INDIVIDUALS READ FROM THE INPUT FILE:            " << numInitPhen << endl;
	ss << "INDIVIDUALS DELETED                              " << numInitPhen - numFinalPhen << endl;
	ss << "INDIVIDUALS LEFT FROM THE INPUT FILE:            " << numFinalPhen << "  (" << numCase << " cases and " << numFinalPhen - numCase << " controls)" << endl;

	cout << ss.str();

	if(!output.init(param_reader, ld_param, data->geno_size(), ss.str())){
		cerr << "Bad setup.  Aborting." << endl;
		exit(0);
	}
	output.setMaxPersonId(data->max_person_ID_size());
}

/**
 * Compute Dandelion throughout.
 * 
 * Two paths:
 * 	1) If window was set then multiple runs are required.
 *  2) If not, then run the whole set.
 */
void Dandelion::process(){

	int window = ld_param->get_dandelion_window();
	
	if(window < 0){

		int run_size = data->geno_size();
		// Set up window with window size.
	
		vector<DandelionSnpInfo> snps;
		char ma, mi, ref;
		string name, chr;
		int pos;
		for(int i=0;i<run_size;i++){
			DandelionSnpInfo singleSNP;
			data->get_allele_codes(i, ma, mi, ref);
			data->get_map_info(i, chr, name, pos);
			singleSNP.name = name;
			singleSNP.position = pos;
			singleSNP.majAllele = ma;
			singleSNP.minAllele = mi;
			singleSNP.refAllele = ref;
			singleSNP.index = i; // not actually used.
			snps.push_back(singleSNP);
		}
		output.writeHeadSet(snps);
		runSet(0, run_size-1);

	}else{
		
		int beg, end = data->geno_size();
		int runSize = 0;
		for(beg=0; beg < end - 1; beg++){
			runSize = min(beg + window - 1 , end -1 );
			
			vector<DandelionSnpInfo> snps;
			char ma, mi, ref;
			string name, chr;
			int pos;
			for(int i=beg;i<=runSize;i++){
				DandelionSnpInfo singleSNP;
				data->get_allele_codes(i, ma, mi, ref);
				data->get_map_info(i, chr, name, pos);
				singleSNP.name = name;
				singleSNP.position = pos;
				singleSNP.chr = chr;
				singleSNP.majAllele = ma;
				singleSNP.minAllele = mi;
				singleSNP.refAllele = ref;
				singleSNP.index = i; // not actually used.
				snps.push_back(singleSNP);
			}
			output.writeHeadSet(snps);
			runSet(beg, runSize);
		}
	}
	
	// Now create an LD program to print out the matrix of r^2 and D.
	LinkageDisequilibrium ld(data);
	ld.enslave(ld_param);

	ld.preProcess();
	ld.process();
	
	
}

/**
 * Perform dandelion computation on a single set of SNPs.  These should
 * be in order in the data set, and should include both begin and end.
 * 
 * @param begin The first SNP to use.
 * @param end The final SNP to use.  All in between begin and end are used.
 * 
 */
void Dandelion::runSet(int begin, int end){
	vector<vector<short> > emInputCs, emInputCb, emInputCn;
	vector<double> phen_nonmissing;
	vector<vector<double> > cov;

	/* Prep covariate matrix */
	vector<double> *tmp = data->get_covariates(0);
	if(tmp != NULL){
		for(unsigned int i=0;i<tmp->size();i++){
			vector<double> a;
			cov.push_back(a);
		}
	}

	vector<short> vCn, vCb, vCs; // temp vectors for each type.

	for(int s1 = begin; s1 <= end; s1++){

		if(data->getDataObject()->isUsable(s1)){
			for(int i=0; i<data->pheno_size(); i++ ){
				// push onto stacks depending on the case/cntrl status.

				if ( equal(data->get_phenotype(i) , 1)){
					vCn.push_back( data->get_data(i)->at(s1) );
					vCb.push_back( data->get_data(i)->at(s1) );
				}else if( equal(data->get_phenotype(i) , 2)){
					vCs.push_back( data->get_data(i)->at(s1) );
					vCb.push_back( data->get_data(i)->at(s1) );
				}
			}

			emInputCs.push_back(vCs);
			vCs.clear();
			emInputCn.push_back(vCn);
			vCn.clear();
			emInputCb.push_back(vCb);
			vCb.clear();
		}
	}

	for(int i=0;i < data->pheno_size();++i){
		phen_nonmissing.push_back(data->get_phenotype(i)-1);

		vector<double> *t = data->get_covariates(i);
		if(t != NULL){
			for(unsigned int j=0;j<t->size();j++){
				cov.at(j).push_back(t->at(j));
			}
		}
	}

	// Data was set up.  Run and get each hap freq.
	EM *emCs, *emCn, *emCb;
	emCs = new EM();
	emCn = new EM();
	emCb = new EM();
	vector<double> caseHapFreqs, cntrlHapFreqs, cmbdHapFreqs;
	try{
		emCs->setup(emInputCs);
		emCs->run();
		caseHapFreqs = emCs->getEMFreqs();

		emCn->setup(emInputCn);
		emCn->run();
		cntrlHapFreqs = emCn->getEMFreqs();

		emCb->setup(emInputCb);
		emCb->run();
		cmbdHapFreqs = emCb->getEMFreqs();
		delete emCs;
		delete emCn;
	}
	catch(EMAlgorithmNoSetup){
		delete emCs;
		delete emCn;
		delete emCb;
		
		cerr << "There was a problem setting up the EM Algorithm.  Please verify that data was complete." << endl;
		return;
	}

	// 
	// Analysis from here down.
	//
	vector<EMPersonalProbsResults> personalProbs = emCb->getPersonalProbabilities();
	delete emCb;
	
	Zaykin zay(ld_param);
	if(ld_param->get_haplo_thresh() >= 0)
		zay.setKeepThresh(ld_param->get_haplo_thresh());
	zay.setPhenotype(phen_nonmissing, cov);
	zay.setup(personalProbs, pow(2, end - begin+1),true);
	
	ZaykinGlobalStatsResults results = zay.runGlobal();
	double pVal = results.pvalue;
	double testStat = results.testStat;
	int degFree = results.degFreedom;
	
	//
	// Write personal probability lines.
	//
	if(ld_param->get_dandelion_pprob()){
		for(unsigned int i=0; i < personalProbs.size(); i++){
			DandelionPProbInfo dpp;
			dpp.personNum = personalProbs.at(i).personId;
			dpp.personID = data->get_person_ID(dpp.personNum);
			dpp.leftHap = prepAlleles(personalProbs.at(i).leftHap, 2, begin, end);
			dpp.rightHap = prepAlleles(personalProbs.at(i).rightHap, 2, begin, end);
			dpp.prob = personalProbs.at(i).prob;
			dpp.affectionStatus = data->get_phenotype(dpp.personNum);
			output.writePProbLine(dpp);
		}
	}


	Zaykin zayHap(ld_param);
	vector<ZaykinStatsInfo> stats;
	if(ld_param->get_haplo_thresh() >= 0)
		zayHap.setKeepThresh(ld_param->get_haplo_thresh());
	zayHap.setPhenotype(phen_nonmissing, cov);
	zayHap.setup(personalProbs, pow(2, end - begin+1),false); // RTG
	stats = zayHap.runAllHaplotypes();

	for(unsigned int i=0; i < caseHapFreqs.size(); i++){
		
		DandelionHaploInfo d;
		d.alleles = prepAlleles(i,2,begin,end);

		d.freqCs = caseHapFreqs.at(i);
		d.freqCn = cntrlHapFreqs.at(i);
		d.freqCb = cmbdHapFreqs.at(i);
		
		try{
			d.p = Statistics::chi2prob(stats.at(i).chiSqStat, stats.at(i).degFree);
			d.z = stats.at(i).chiSqStat;
			d.OR = stats.at(i).OR;
			d.UCI = stats.at(i).UCI;
			d.LCI = stats.at(i).LCI;
		}catch(...){
			
			d.p = -1;
			d.z = -1;
			d.OR = 0;
			d.UCI = 0;
			d.LCI = 0;
		}

		output.writeLine(d);
	}

	output.writeStatisticsLine(pVal, testStat, degFree);

}

/**
 * Compute the statistics for this haplotype.
 *
 * @deprecated This is currently a simple chi2 test.  Needs to be done with LR.  This code
 * is NOT tested and should never be included in production app.
 * 
 */
void Dandelion::computeHaplotypeTest(DandelionHaploInfo &d){

	double pHat, s;

	double numCntrl = numFinalPhen - numCase;

	pHat = (numCase * d.freqCs + (numCntrl) * d.freqCn) / numFinalPhen;
	s = sqrt(pHat * (1-pHat) * numFinalPhen / (numCase * numCntrl));
	d.z = (d.freqCs - d.freqCn)/s;
	try {

		d.p = Statistics::chi2prob(d.z * d.z, 1.0);

	}catch(...){

		d.p = 2.0;

	}

	d.OR = sqrt( ( 2*d.freqCs * numCase + .5)*( 2*(1-d.freqCn) * numCntrl )/(( 2*d.freqCn * numCntrl + 0.5 )*( 2*(1-d.freqCs) * numCase + 0.5 )) );
	double lor = sqrt( 1/(2*d.freqCs * numCase + .5) + 1/( 2*(1-d.freqCn) * numCntrl ) + 1/( 2*d.freqCn * numCntrl + 0.5) + 1/(2*(1-d.freqCs) * numCase + 0.5 )   );
	d.LCI =  exp(log(d.OR) - (1.96*lor));
	d.UCI =  exp(log(d.OR) + (1.96*lor));
	if(d.LCI > d.UCI){
		double t = d.LCI;
		d.LCI = d.UCI;
		d.UCI = t;
	}

}

/**
 * Prep allele list by breaking number apart and pulling either major or minor allele.
 *
 * @param row An integer that signifies a haplotype.
 * @param divisor The coding for divisor.  Default is 2, but some codings use something else.
 * @return vector<char> List of haplotypes.
 */
vector<char> Dandelion::prepAlleles(unsigned int row, int divisor, int beg, int end){

	vector<char> ret;

	char ma, mi, t;
	int cnt = 0;

	for(int i = beg;i<=end;i++){
		data->get_allele_codes(i, ma, mi, t);
		ret.push_back(ma);
	}
	
	cnt = end - beg;
	while(row > 0){
	
		data->get_allele_codes(beg + cnt, ma, mi, t);

		if(row % 2 != 0){
			ret.at(cnt) = mi;
		}
		row /= divisor;
		
		cnt--;
	}
	
	return ret;
}

void Dandelion::enslave(EngineParamReader *e){
	haveOwner = true;
	delete ld_param;
	ld_param = e;
}

void Dandelion::test(){
	// Keep this blank for now.
}



