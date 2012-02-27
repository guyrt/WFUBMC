//      snpgwa.cpp
//      
//      Copyright 2010 Richard T. Guy <guy@cs.toronto.edu>
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

#include "snpgwa.h"

Snpgwa::Snpgwa(){
	this->param_reader = ParamReader::Instance();
	this->data = new DataAccess;
	this->data->init(NULL);
	this->snp_param = new EngineParamReader;
}

Snpgwa::~Snpgwa(){
	delete_my_innards();
}

void Snpgwa::delete_my_innards(){
	
	delete this->param_reader;
	delete this->data;
	data = NULL;
	delete this->snp_param;
	snp_param = NULL;
}

void Snpgwa::enslave(EngineParamReader *e){
	cerr << "Enslavement not supported for SNPGWA." << endl;
}

void Snpgwa::test(){
	// Could contain all sorts of tests.
}

/**
 * Primary initialization function.  Assumes that the parameters have been
 * read and that the data object has not been initialized.  Reads data
 * and performs most preprocessing on that data.
 * 
 * After init() runs, the data has been loaded and is ready.
 * 
 * Specific tasks:
 * Read data
 * Make categorical
 * Remove individuals missing phenotype.
 * 		 
 */
void Snpgwa::init(){

	snp_param->read_parameters(param_reader->get_engine_specific_params());

	if(param_reader->get_linkage_map_file().compare("none") == 0){
		cerr << "SNPGWA requires that you enter a map file.  Aborting." << endl;
//		exit(0);
	}


	initializeReader(); // defined in engine.h

	Logger::Instance()->init(param_reader->get_out_file() + ".log");

	#if DB_V_SNPGWA
		cout << "Starting process" <<endl;
	#endif

	reader->process(data->getDataObject(), param_reader);
	// Read the number of individuals and genotypes.
	numInitSNPs = data->geno_size();
	numInitPhen = data->pheno_size();
	#if DB_V_SNPGWA
		cout << "Read complete.  Starting remove." <<endl;
	#endif
	data->getDataObject()->remove_phenotype(numeric_limits<double>::max());
	data->getDataObject()->remove_covariate(numeric_limits<double>::max());
	data->getDataObject()->remove_phenotype(0);
	#if DB_V_SNPGWA
		cout << "Remove complete.  Starting final prep." <<endl;
	#endif

	data->getDataObject()->prep_data(param_reader);
	/* At this point, all individuals that have phenotype of 0 or of '.' are removed from data.*/
	/* Also, the minor allele has been set (assuming that there wasn't a reference allele) */
	#if DB_V_SNPGWA
		cout << "Starting categorical." <<endl;
	#endif
	try{
		data->make_categorical(1,2);
	}catch(DataException d){
		cerr << "Data Exception with message: " << endl << d.message << endl;
		exit(0);
	}catch(...){
		cerr << "Unknown exception making categorical." << endl;
		exit(0);
	}
	delete reader; // No longer needed.

	// Read the number of phenotypes remaining and get the num of cases.
	numFinalPhen = data->pheno_size();
	#if DB_V_SNPGWA
		cout << "Finished init." <<endl;
	#endif
}

/**
 * Set up output engine, alert user about final data size, and prepare to
 * run tests.
 */
void Snpgwa::preProcess(){

	int numCase = 0;
	for(int i=0;i < data->pheno_size(); ++i)
		if(equal(data->get_phenotype(i), 2)) numCase++;

	stringstream ss;
	ss << "INDIVIDUALS READ FROM THE INPUT FILE:            " << numInitPhen << endl;
	ss << "INDIVIDUALS DELETED                              " << numInitPhen - numFinalPhen << endl;
	ss << "INDIVIDUALS LEFT FROM THE INPUT FILE:            " << numFinalPhen << "  (" << numCase << " cases and " << numFinalPhen - numCase << " controls)" << endl;

	cout << ss.str();

	if(!out.init(param_reader, snp_param, data->max_map_size(), data->geno_size(), ss.str())){
		cerr << "Error opening output files.  Aborting." <<
		endl << "If you did not expect this error, and you are on a Linux/Unix machine, please verify that you have permission to open the file requested." << endl;
		exit(0);
	}
}

/**
 * Primary entry point for SNPgwa processing.  Loop through SNPs, potentially 
 * in parallel, and write output for each SNP.  Utilizes several classes.
 * 
 * @see PopStats
 * @see GenoStats
 * @see HaploStats
 * @see SnpgwaOutput
 */
void Snpgwa::process(){

	#if RUN_IN_PARALLEL
	#pragma omp parallel
	{
	#endif
	int sz = data->geno_size();
	#if RUN_IN_PARALLEL
	#pragma omp for schedule(static, 10)
	#endif
	for(int i=0;i < sz;i++){

		SnpInfo s;
		s.index = i+param_reader->get_begin();

		if(param_reader->get_linkage_map_file().compare("none") != 0){
			data->get_map_info(i, s.chr, s.name, s.position);

			// calc diff.
			if(i < sz-1){

				string t, c;
				int p = 0;
	
				data->get_map_info(i+1, c, t, p);
				if(c.compare(s.chr) == 0){
					s.diff = p-s.position;
				}else{
					s.diff = -999;
				}

			}else{
				s.diff = -999;
			}
		}

		char ma, mi, mr;
		data->get_allele_codes(i,ma,mi, mr);
		s.majAllele = ma;
		s.minAllele = mi;
		s.refAllele = mr;
		if(s.majAllele == ' '){
			s.majAllele = '.';
		}
		if(s.minAllele == ' '){
			s.minAllele = '.';
		}
		if(s.refAllele == ' '){
			s.refAllele = '.';
		}

		GenoStatsResults ge;
		PopStatsResults p;
		HaploStatsResults hr;

		if(data->getDataObject()->isUsable(i)){

			PopStats pop_calc(data);
			GenoStats g(data, snp_param);
			LinkageDisequilibrium ld(data);
			ld.enslave(snp_param);
			LinkageMeasures lr;
			g.prepGenoStatsForOutput(i,ge);
			
			pop_calc.prepPopStatsForOutput(i,p);

			if(snp_param->get_snpgwa_dohap()){
				HaploStats h(data, snp_param);
				if(snp_param->get_haplo_thresh() >= 0)
					h.setHaploThresh(snp_param->get_haplo_thresh());
				h.prepHaploStatsForOutput(i, hr);
				if(i + 1 < data->geno_size()){
	
					ld.dprimeOnPair(i, i+1, lr);
					hr.rsquare = lr.rsquare;
					hr.dprime = lr.dPrime;
				}else{
					hr.rsquare = hr.dprime = -1.0;
				}
			}

			else{
				hr.rsquare = hr.dprime = -1.0;
				for(int ii=0;ii < 4;ii++){
					hr.twoMarkerCaseFreq.push_back(0);
					hr.twoMarkerCntrlFreq.push_back(0);
				}
				for(int ii=0;ii < 8;ii++){
					hr.threeMarkerCaseFreq.push_back(0);
					hr.threeMarkerCntrlFreq.push_back(0);
				}
			}

		}else{
			// Not usable SNP.  Throw blanks.
			initToZero(p,hr,ge);
		}

		out.writeLine(i,s,p,hr, ge);

	}

	#if RUN_IN_PARALLEL
	}
	#endif
	out.close();

}

/**
 * Initialize all members of the group to their default values.
 * 
 * @param p  A population results object to be cleared.
 * @param r  A HaplotypeStatsResults object to be cleared.
 * @param ge  A GenotypeStatsResults object to be cleared.
 */
void Snpgwa::initToZero(PopStatsResults &p, HaploStatsResults &r, GenoStatsResults &ge){

	p.caseCount = 0;
	p.cntrlCount = 0;
	p.caseRefFreq = 0;
	p.cntrlRefFreq = 0;
	p.pMissingCombined = 1.0;
	p.pMissingCase = 1.0;
	p.pMissingCntrl = 1.0;
	p.missingPVal = 2.0;
	p.missingOR = -1.0;
	p.cntrlPP = 0;
	p.casePP = 0;
	p.expcntrlPP = 0;
	p.expcasePP = 0;
	p.expcmbdPP = 0;
	p.cntrlPQ = 0;
	p.casePQ = 0;
	p.expcntrlPQ = 0;
	p.expcasePQ = 0;
	p.expcmbdPQ = 0;
	p.cntrlQQ = 0;
	p.caseQQ = 0;
	p.expcntrlQQ = 0;
	p.expcaseQQ = 0;
	p.expcmbdQQ = 0;

	p.cmbdTestStat  = -1;
   	p.caseTestStat = -1;
   	p.cntrlTestStat = -1;

	p.cmbdPVal = 2;
	p.casePVal = 2;
	p.cntrlPVal = 2;
	p.cmbdExactPVal  = 2;
   	p.caseExactPVal = 2;
	p.cntrlExactPVal = 2;

	ge.twodegTestStat = -1;
	ge.twodegPVal = 2;

	ge.domTestStat = -1;
	ge.domPVal = 2;
	ge.domOR = -1;
	ge.domLCI = -1;
	ge.domUCI = -1;
	ge.domSens = -1;
	ge.domSpec = -1;
	ge.domCStat = -1;

	ge.addTestStat = -1;
	ge.addPVal = 2;
	ge.addOR = -1;
	ge.addLCI = -1;
	ge.addUCI = -1;
	ge.addSensNNRN = -1;
	ge.addSpecNNRN = -1;
	ge.addSensNNRR = -1;
	ge.addSpecNNRR = -1;
	ge.addSensNRRR = -1;
	ge.addSpecNRRR = -1;
	ge.addCStat = -1;

	ge.recTestStat = -1;
	ge.recPVal = 2;
	ge.recOR = -1;
	ge.recLCI = -1;
	ge.recUCI = -1;
	ge.recSens = -1;
	ge.recSpec = -1;
	ge.recCStat = -1;

	ge.lofTestStat = -1;
	ge.lofPVal = 2;

	r.dprime = -1;
	r.rsquare = -1;

	r.allelicChiS = -1;
	r.allelicDF = -1;
	r.allelicPval = 2.0;

	r.twoMarkerChiS = -1;
	r.twoMarkerDF = -1;
	r.twoMarkerPval = 2;

	r.threeMarkerChiS = -1;
	r.threeMarkerDF = -1;
	r.threeMarkerPval = 2;

	for(int i=0;i < 4;i++){
		r.twoMarkerCaseFreq.push_back(0);
		r.twoMarkerCntrlFreq.push_back(0);
	}
	for(int i=0;i < 8;i++){
		r.threeMarkerCaseFreq.push_back(0);
		r.threeMarkerCntrlFreq.push_back(0);
	}

}
