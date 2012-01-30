#include "haplostats.hh"

HaploStats::HaploStats(DataAccess *d, EngineParamReader *p){
	
	params = p;
	data = d;
	twoMarkerPval = threeMarkerPval = 2.0;
	twoMarkerChiS = threeMarkerChiS = -1;
	twoMarkerDF = threeMarkerDF = -1;
	allelicPval = 2.0;
	allelicChiSq = -1;
	allelicDF = -1;
	
	haploThresh = -1;
}

void HaploStats::prepHaploStatsForOutput(int snp, HaploStatsResults &res){
	
	calculateAllelic(snp);
	
	#if DEBUG_HAPL_PROGRESS
	cout << "Two marker started on " << snp << endl;
	#endif
	if(snp+1 < data->geno_size() && data->getDataObject()->isUsable(snp+1)){
		calculateTwoMarker(snp, snp+1);
	}else{
		twoMarkerPval = 2.0;
		for(int i=0;i<4;i++) twoMarkerCaseHapFreq.push_back(0.0);
		twoMarkerCmbdHapFreq = twoMarkerCaseHapFreq;
		twoMarkerCntrlHapFreq = twoMarkerCaseHapFreq;	
	}
	#if DEBUG_HAPL_PROGRESS
	cout << "Two marker finished on " << snp << endl;
	cout << "Three marker started on " << snp << endl;
	#endif
	if(snp+2 < data->geno_size() && data->getDataObject()->isUsable(snp+1) && data->getDataObject()->isUsable(snp+2)){
		calculateThreeMarker(snp, snp+1, snp+2);
	}else{
		threeMarkerPval = 2.0;
		for(int i=0;i<8;i++) threeMarkerCaseHapFreq.push_back(0.0);
		threeMarkerCmbdHapFreq = threeMarkerCaseHapFreq;
		threeMarkerCntrlHapFreq = threeMarkerCaseHapFreq;		
	}
	#if DEBUG_HAPL_PROGRESS
	cout << "Three marker ended on " << snp << endl;
	#endif
	res.allelicPval = allelicPval;
	res.allelicChiS = allelicChiSq;
	res.allelicDF = allelicDF;

	res.twoMarkerChiS = twoMarkerChiS;
	res.twoMarkerDF = twoMarkerDF;
	res.twoMarkerPval = twoMarkerPval;
	res.twoMarkerCaseFreq = twoMarkerCaseHapFreq;
	res.twoMarkerCntrlFreq = twoMarkerCntrlHapFreq;
	
	res.threeMarkerChiS = threeMarkerChiS;
	res.threeMarkerDF = threeMarkerDF;
	res.threeMarkerPval = threeMarkerPval;
	res.threeMarkerCaseFreq = threeMarkerCaseHapFreq;
	res.threeMarkerCntrlFreq = threeMarkerCntrlHapFreq;
	
}

/**
 * Calculate allelic test.  This can be conceived as a single SNP haplotype test
 * where all of hte haplotpyes are known with probability one.  It leverages the 
 * code built for larger hapltoypes, but EM is never called.
 * 
 */
void HaploStats::calculateAllelic(int snp){
	
	// Build a series of haplotype results.
	vector<EMPersonalProbsResults> haps;
	
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
	
	int cnt = 0;
	for(int i=0; i<data->pheno_size(); i++ ){
		if(data->get_data(i)->at(snp) != 0){
			EMPersonalProbsResults e;
			e.personId = cnt++;
			e.prob = 1.0;
			
			switch (data->get_data(i)->at(snp)){
				case 1:
					e.leftHap = 0;
					e.rightHap = 0;
				break;
				case 2:
					e.leftHap = 0;
					e.rightHap = 1;
				break;
				case 3:
					e.leftHap = 1;
					e.rightHap = 0;
				break;
				case 4:
					e.leftHap = 1;
					e.rightHap = 1;
				break;
				default:
				break;	
			}
			haps.push_back(e);
			
			phen_nonmissing.push_back(data->get_phenotype(i)-1);
			vector<double> *t = data->get_covariates(i);
			if(t != NULL){
				for(unsigned int j=0;j<t->size();j++){
					cov.at(j).push_back(t->at(j));
				}
			}
		}
	}

	Zaykin zay(params);
	if(haploThresh >= 0)
		zay.setKeepThresh(haploThresh);
	zay.setPhenotype(phen_nonmissing, cov);
	zay.setErrorInformation(snp, cov.size(), data);
	zay.setup(haps, 2, true);
	
	try{
		ZaykinGlobalStatsResults results = zay.runGlobal();
		allelicPval = results.pvalue;
		allelicChiSq = results.testStat;
		allelicDF = results.degFreedom;
	}catch(...){
		stringstream ss;
		ss << "Failed to compute global allelic test for SNP " << snp << endl;
		Logger::Instance()->writeLine(ss.str());
	}
	if(allelicPval != allelicPval){
		allelicPval = 2.0; // check for isNAN
		allelicChiSq = -1;
		allelicDF = -1;
	}
	
}

/**
 * Calculate the two marker haplotype test.
 * 
 * Steps:
 * 	1) Prep and run the EM algorithm for cases and for controls.
 *  2) Calculate haplotypic associatoin
 * 
 */
void HaploStats::calculateTwoMarker(int s1, int s2){
	
	/* Set up three EM algorithms and run them. */
	EM emCase, emCntrl, emCmbd;
	int degreesOfFreedom = 0;
	
	vector<short> v1Cs, v2Cs, v1Cn, v2Cn, v1Cb, v2Cb;
	vector<vector<double> > cov;
	vector<double> phen_nonmissing;
	
	phen_nonmissing.clear();
	/* Prep covariate matrix */
	vector<double> *tmp = data->get_covariates(0);
	if(tmp != NULL){
		for(unsigned int i=0;i<tmp->size();i++){
			vector<double> a;
			cov.push_back(a);
		}
	}
	for(int i=0; i<data->pheno_size(); i++ ){
		// push onto stacks depending on the case/cntrl status.
		if(data->get_data(i)->at(s1) != 0 && data->get_data(i)->at(s2) != 0){
			// Data to be used only if not missing in both.
			if (data->get_phenotype(i) == 1){
				v1Cn.push_back( data->get_data(i)->at(s1) );
				v2Cn.push_back( data->get_data(i)->at(s2) );
				v1Cb.push_back( data->get_data(i)->at(s1) );
				v2Cb.push_back( data->get_data(i)->at(s2) );
			}else if(data->get_phenotype(i) == 2){
				v1Cs.push_back( data->get_data(i)->at(s1) );
				v2Cs.push_back( data->get_data(i)->at(s2) );
				v1Cb.push_back( data->get_data(i)->at(s1) );
				v2Cb.push_back( data->get_data(i)->at(s2) );
			}
			phen_nonmissing.push_back(data->get_phenotype(i)-1);
			
			vector<double> *t = data->get_covariates(i);
			if(t != NULL){
				for(unsigned int j=0;j<t->size();j++){
					cov.at(j).push_back(t->at(j));
				}
			}
		}
	}
	
	vector<vector<short> > vCs;
	vCs.push_back(v1Cs);
	vCs.push_back(v2Cs);
	emCase.setup(vCs);
	//emCase.setup(v1Cs, v2Cs);
	vCs.clear();
	vector<vector<short> > vCn;
	vCn.push_back(v1Cn);
	vCn.push_back(v2Cn);
	emCntrl.setup(vCn);
	vCn.clear();
	vector<vector<short> > vCb;
	vCb.push_back(v1Cb);
	vCb.push_back(v2Cb);
	emCmbd.setup(vCb);
	vCb.clear();
	
	emCase.run();
	emCntrl.run();
	emCmbd.run();

	try{
		twoMarkerCaseHapFreq = emCase.getEMFreqs();
		twoMarkerCntrlHapFreq = emCntrl.getEMFreqs();
		twoMarkerCmbdHapFreq = emCmbd.getEMFreqs();
	}catch (EMAlgorithmFailureException){
		
		// Log message
		stringstream ss;
		ss << "Error SNP " << s1 << ".  Unable to calculate two marker haplotype test." << endl;
		Logger::Instance()->writeLine(ss.str());
		twoMarkerPval = 2.0;
		twoMarkerCaseHapFreq.clear();
		twoMarkerCntrlHapFreq.clear();
		twoMarkerCmbdHapFreq.clear();
		for(int i=0;i < 4;i++){
			twoMarkerCaseHapFreq.push_back(0);
			twoMarkerCntrlHapFreq.push_back(0);
			twoMarkerCmbdHapFreq.push_back(0);
		}
		return;
	}
	
	// Calculate degrees of freedom
	for(unsigned int i=0; i<twoMarkerCaseHapFreq.size(); i++){
		if(twoMarkerCaseHapFreq.at(i) > EPS() || twoMarkerCntrlHapFreq.at(i) > EPS()) degreesOfFreedom++;
	}
	
	vector<EM *> ems;
	ems.push_back(&emCmbd);ems.push_back(&emCase); ems.push_back(&emCntrl);
	twoMarkerPval = computeGlobalZaykinStatisic(emCmbd, 2, twoMarkerChiS, twoMarkerDF, phen_nonmissing, cov, s1);
	
}

/**
 * Calculate the three marker haplotype test.
 * 
 * Steps:
 * 	1) Prep and run the EM algorithm for cases and for controls.
 *  2) Calculate haplotypic associatoin
 * 
 */
void HaploStats::calculateThreeMarker(int s1, int s2, int s3){
	
	/* Set up three EM algorithms and run them. */
	EM emCase, emCntrl, emCmbd;
		
	vector<short> v1Cs, v2Cs, v3Cs, v1Cn, v2Cn, v3Cn, v1Cb, v2Cb, v3Cb;
	vector<vector<double> > cov;
	vector<double> phen_nonmissing;
	
	phen_nonmissing.clear();
	/* Prep covariate matrix */
	vector<double> *tmp = data->get_covariates(0);
	if(tmp != NULL){
		for(unsigned int i=0;i<tmp->size();i++){
			vector<double> a;
			cov.push_back(a);
		}
	}
	for(int i=0; i<data->pheno_size(); i++ ){
		// push onto stacks depending on the case/cntrl status.
		if(data->get_data(i)->at(s1) != 0 && data->get_data(i)->at(s2) != 0 && data->get_data(i)->at(s3) != 0){
			// Data to be used only if not missing in both.
			if (data->get_phenotype(i) == 1){
				v1Cn.push_back( data->get_data(i)->at(s1) );
				v2Cn.push_back( data->get_data(i)->at(s2) );
				v3Cn.push_back( data->get_data(i)->at(s3) );
				v1Cb.push_back( data->get_data(i)->at(s1) );
				v2Cb.push_back( data->get_data(i)->at(s2) );
				v3Cb.push_back( data->get_data(i)->at(s3) );
			}else if(data->get_phenotype(i) == 2){
				v1Cs.push_back( data->get_data(i)->at(s1) );
				v2Cs.push_back( data->get_data(i)->at(s2) );
				v3Cs.push_back( data->get_data(i)->at(s3) );
				v1Cb.push_back( data->get_data(i)->at(s1) );
				v2Cb.push_back( data->get_data(i)->at(s2) );
				v3Cb.push_back( data->get_data(i)->at(s3) );
			}
			phen_nonmissing.push_back(data->get_phenotype(i)-1);
			vector<double> *t = data->get_covariates(i);
			if(t != NULL){
				for(unsigned int j=0;j<t->size();j++){
					cov.at(j).push_back(t->at(j));
				}
			}
		}
	}

	#if DEBUG_HAPL_PROGRESS
	cout << "Three marker run all start " << s1 << endl;
	#endif
	
	vector<vector<short> > temp;
	temp.push_back(v1Cs);
	temp.push_back(v2Cs);
	temp.push_back(v3Cs);
	emCase.setup(temp);
	emCase.run();
		
	temp.clear();
	temp.push_back(v1Cn);
	temp.push_back(v2Cn);
	temp.push_back(v3Cn);
	emCntrl.setup(temp);
	emCntrl.run();
		
	temp.clear();
	temp.push_back(v1Cb);
	temp.push_back(v2Cb);
	temp.push_back(v3Cb);
	emCmbd.setup(temp);
	emCmbd.run();
	
	#if DEBUG_HAPL_PROGRESS
	cout << "Three marker run all end " << s1 << endl;
	#endif
	
	try{
		threeMarkerCaseHapFreq = emCase.getEMFreqs();
		threeMarkerCntrlHapFreq = emCntrl.getEMFreqs();
		threeMarkerCmbdHapFreq = emCmbd.getEMFreqs();
	}catch(EMAlgorithmFailureException){
		// Log message
		stringstream ss;
		ss << "Error SNP " << s1 << ".  EM algorithm failure during three marker haplotype test." << endl;
		
		Logger::Instance()->writeLine(ss.str());
		
		threeMarkerPval = 2.0;
		threeMarkerCaseHapFreq.clear();
		threeMarkerCntrlHapFreq.clear();
		threeMarkerCmbdHapFreq.clear();
		for(int i=0;i < 9;i++){
			threeMarkerCaseHapFreq.push_back(0);
			threeMarkerCntrlHapFreq.push_back(0);
			threeMarkerCmbdHapFreq.push_back(0);
		}
		return;
	}catch(...){
		stringstream ss;
		ss << "Error SNP " << s1 << ".  Unable to calculate three marker haplotype test." << endl;
		
		Logger::Instance()->writeLine(ss.str());
		
		threeMarkerPval = 2.0;
		threeMarkerCaseHapFreq.clear();
		threeMarkerCntrlHapFreq.clear();
		threeMarkerCmbdHapFreq.clear();
		for(int i=0;i < 9;i++){
			threeMarkerCaseHapFreq.push_back(0);
			threeMarkerCntrlHapFreq.push_back(0);
			threeMarkerCmbdHapFreq.push_back(0);
		}
		return;
	}
	
	#if DEBUG_HAPL_PROGRESS
		cout << "Three marker dof s " << s1 << endl;
	#endif

	threeMarkerPval = computeGlobalZaykinStatisic(emCmbd, 3, threeMarkerChiS, threeMarkerDF, phen_nonmissing, cov, s1);
	#if DEBUG_HAPL_PROGRESS
		cout << "Three marker zay e " << s1 << endl;
	#endif
}

/*
 * Get Zaykin's p-value.
 */
double HaploStats::computeGlobalZaykinStatisic(EM &em, int size, double &testStat, int &degFreedom, 
		const vector<double> &phenotype, const vector<vector<double> > &cov, int snp){
	
	vector<EMPersonalProbsResults> r;
	try{
		r = em.getPersonalProbabilities();
	}catch(...){
		testStat = 2.0;
		degFreedom = -1.0;
		
		return 2.0;
	}
	
	
	Zaykin zay(params);
	if(haploThresh >= 0)
		zay.setKeepThresh(haploThresh);
	zay.setErrorInformation(snp, cov.size(), data);
	zay.setPhenotype(phenotype, cov);
	try{
		zay.setup(r, static_cast<int>(pow(2, size)), true);
		ZaykinGlobalStatsResults results = zay.runGlobal();
		
		testStat = results.testStat;
		degFreedom = results.degFreedom;
		return results.pvalue;
	}catch(...){
		cout << "Failed to compute global stat." << endl;
		return 2.0;
	}
	
}


