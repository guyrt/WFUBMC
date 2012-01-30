#include "popstats.hh"

PopStats::PopStats(DataAccess *d){
	data = d;
	lastBuild = -1;
	reinit();
}

void PopStats::reinit(){
	
	numMissingCase = numMissingCntrl = numMissingTotal = 0;
	numTotalCases = numTotalCntrls = 0;
	casesNonMissing = cntrlsNonMissing = 0;
	minorAlleleFreqCases = minorAlleleFreqCntrls = 0.0;
	ppCases = pqCases = qqCases = ppCntrls = pqCntrls = qqCntrls = 0;
	caseExpTotal[0] = caseExpTotal[1] = caseExpTotal[2] = 0;
	caseChiSquare = -1;
	caseChiSquarePValue = caseExactTestPVal = 2.0;
	cntrlExpTotal[0] = cntrlExpTotal[1] = cntrlExpTotal[2] = 0;
	cntrlChiSquare = -1;
	cntrlChiSquarePValue = cntrlExactTestPVal = 2.0;
	cmbdExpTotal[0] = cmbdExpTotal[1] = cmbdExpTotal[2] = 0;
	cmbdChiSquare = -1;
	cmbdChiSquarePValue = cmbdExactTestPVal = 2.0;
	snp_data.clear();
	phen_data.clear();
	percentMissingTotal = 0.0;
}

/*
 * Perform all computations in the populationStatistics framework and send
 * their results to the output mechanism o
 */
bool PopStats::prepPopStatsForOutput(int snp, PopStatsResults &results){
	bool ret = true;
	pullData(snp);
	ret &= minorAlleleFreq(snp);
	ret &= numCategories(snp);
	ret &= missingStats(snp);
	ret &= hwEquilibruim(snp);

//	results.caseCount = numTotalCases;
//	results.cntrlCount = numTotalCntrls;
	results.caseCount = numTotalCases - numMissingCase;  //or = casesNonMissing;
	results.cntrlCount = numTotalCntrls - numMissingCntrl;// or = cntrlsNonMissing;
	results.caseRefFreq = minorAlleleFreqCases;
	results.cntrlRefFreq = minorAlleleFreqCntrls;
	results.pMissingCombined = percentMissingTotal;
	results.pMissingCase = percentMissingCases;
	results.pMissingCntrl = percentMissingControls;
	results.missingPVal = differentialMissingPval;
	results.missingOR = differentialMissingOR;
	results.cntrlPP = ppCntrls;
	results.casePP = ppCases;
	results.expcmbdPP = cmbdExpTotal[0];
	results.expcasePP = caseExpTotal[0];
	results.expcntrlPP = cntrlExpTotal[0];
	results.cntrlPQ = pqCntrls;
	results.casePQ = pqCases;
	results.expcmbdPQ = cmbdExpTotal[1];
	results.expcasePQ = caseExpTotal[1];
	results.expcntrlPQ = cntrlExpTotal[1];
	results.cntrlQQ = qqCntrls;
	results.caseQQ = qqCases;
	results.expcmbdQQ = cmbdExpTotal[2];
	results.expcaseQQ = caseExpTotal[2];
	results.expcntrlQQ = cntrlExpTotal[2];
	results.cmbdPVal = cmbdChiSquarePValue;
	results.casePVal = caseChiSquarePValue;
	results.cntrlPVal = cntrlChiSquarePValue;
	results.cmbdExactPVal = cmbdExactTestPVal;
	results.caseExactPVal = caseExactTestPVal;
	results.cntrlExactPVal = cntrlExactTestPVal;
	
	results.caseTestStat = caseChiSquare;
	results.cntrlTestStat = cntrlChiSquare;
	results.cmbdTestStat = cmbdChiSquare;
	
	return ret;
}

/**
 * Repopulate the vector with the new SNP if necessary.
 * If we run all our statistics on a single SNP at one time,
 * we save serious time.
 * 
 * Return true if it was recomputed.
 */
bool PopStats::pullData(int snp){
	bool retVal;
	#pragma omp critical
	{

		if(lastBuild != snp){
			
			#if DEBUG_POPSTATS_PULL
			cout << "Rebuilding with SNP " << snp << endl;
			#endif
			
			retVal = true;
			lastBuild = snp;
			reinit();
			for(int i=0; i<data->pheno_size(); i++ ){
				if(data->get_data(i)->at(snp) > 0){
					snp_data.push_back(data->get_data(i)->at(snp));
					phen_data.push_back(data->get_phenotype(i));
				}else{
					if(data->get_phenotype(i) == 1){
						numMissingCntrl++;
					}else if(data->get_phenotype(i) == 2){
						numMissingCase++;
					}
					numMissingTotal++;
				}
			}
		}else{
			retVal = false;
		}
	}
	return retVal;
}

/**
 * Calculate the minor allele frequency for a given SNP.
 */
bool PopStats::minorAlleleFreq(int snp){
	double ma1 = 0.0, ma2 = 0.0;
	double sum1 = 0.0, sum2 = 0.0;

	/* Compile data. */
	pullData(snp);
	for(unsigned int i = 0;i<snp_data.size();++i){
		if(snp_data.at(i) > 0){

			if( phen_data.at(i) == 1 ){
				sum1+=2;
				switch(snp_data.at(i)){
					case 1:
						ppCntrls++;
					break;
					case 2:
						pqCntrls++;
						ma1++;
					break;
					case 3:
						pqCntrls++;
						ma1++;
					break;
					case 4:
						qqCntrls++;
						ma1+=2;
					break;
					default:
					//empty
					break;
				}
			}else if(phen_data.at(i) == 2){
				sum2+=2;
				switch(snp_data.at(i)){
					case 1:
						ppCases++;
					break;
					case 2:
						pqCases++;
						ma2++;
					break;
					case 3:
						pqCases++;
						ma2++;
					break;
					case 4:
						qqCases++;
						ma2+=2;
					break;
					default:
					//empty
					break;
				}
			}
		}
	}
	
	minorAlleleFreqCases = ma2/sum2;
	minorAlleleFreqCntrls = ma1/sum1;
	return true;
}

/**
 * Compute number of cases and controls with information at this locus.
 */
bool PopStats::numCategories(int snp){

	pullData(snp);
	numTotalCases = numTotalCntrls = 0;  //Do these need to be set to zero again?
	/*for(unsigned int i=0;i<snp_data.size();++i){
		if(snp_data.at(i) != 0){
			if(phen_data.at(i) == 1){
				numTotalCntrls++;
			}else if(phen_data.at(i) == 2){
				numTotalCases++;
			}
		}
	}*/
	
	casesNonMissing = qqCases+pqCases+ppCases;
	cntrlsNonMissing = qqCntrls+pqCntrls+ppCntrls;
	numTotalCases = casesNonMissing+numMissingCase;
	numTotalCntrls = cntrlsNonMissing+numMissingCntrl;
	
	return true;
}

/*
 * Calculate the percent missing for combined, case, and control.
 * Also calculate an odds ratio and p-value for differential missing.
 * 
 * This function makes use of the numbers computed in the pullData function.
 * 
 */
bool PopStats::missingStats(int snp){
	
	if(pullData(snp)){
		numCategories(snp);
	}
	int numTotalIndiv = data->pheno_size();
	percentMissingTotal = static_cast<double>(numMissingTotal)/numTotalIndiv;
	percentMissingControls = static_cast<double>(numMissingCntrl)/numTotalCntrls; // Note: this assumes missing phenotype people included.
	percentMissingCases = static_cast<double>(numMissingCase)/numTotalCases;
	
	// Compute p-val and odds ratio.
	if(percentMissingCases < 1.0 && percentMissingControls < 1.0){
		
		double fNotMissCase = numTotalCases-numMissingCase + 0.5;  // or = casesNonMissing + 0.5;
		double fNotMissCntrl = numTotalCntrls-numMissingCntrl + 0.5; // or = cntrlsNonMissing + 0.5;
		double fTotCase = numTotalCases + 1.0;
		double fTotCntrl = numTotalCntrls + 1.0;
		double fMissCase = numMissingCase + 0.5 ;
		double fMissCntrl = numMissingCntrl + 0.5 ;

		
		double chi_orc = fNotMissCase * log(fNotMissCase) + fMissCase * log(fMissCase);
		chi_orc += fNotMissCntrl * log(fNotMissCntrl) + fMissCntrl * log(fMissCntrl);
		chi_orc -= (fTotCase * log( fTotCase) + fTotCntrl * log(fTotCntrl) + (fNotMissCase + fNotMissCntrl)*log(fNotMissCase + fNotMissCntrl) + (fMissCase+fMissCntrl)*log(fMissCase+fMissCntrl));
		chi_orc += (fTotCase+fTotCntrl)*log(fTotCase+fTotCntrl);
		chi_orc *=2;
	
		if(chi_orc < -0.001){
			
			stringstream ss;
			ss << "Missing statistics for SNP " << snp << " Chi value very negative in differential missing statistic: " << chi_orc << "." << endl;
			Logger::Instance()->writeLine(ss.str());
			
			differentialMissingPval = 2.0;
			differentialMissingOR = -1.0;
			return false;
		}
		if(chi_orc < 0) chi_orc = 0; // Fixes problem of very small but negative numbers.
		
		try{
			differentialMissingPval = Statistics::chi2prob(chi_orc,1);
		}catch(StatsException){
			differentialMissingPval = 2.0;
			stringstream ss;
			ss << "Error, SNP " << snp << " missing p-value calculation not computable." << endl;
			Logger::Instance()->writeLine(ss.str());
		}
		//caseCntlUnknownOdds = (1-nmiss_casec) / nmiss_casec / ((1-nmiss_cntc) / nmiss_cntc);
		
		fMissCase /= fTotCase;
		fMissCntrl /= fTotCntrl;
		differentialMissingOR = (1-fMissCase) / (fMissCase) / ((1-fMissCntrl) / (fMissCntrl));
//		differentialMissingOR = (1-percentMissingCases+.5/numTotalIndiv) / (percentMissingCases+.5/numTotalIndiv) / ((1-percentMissingControls+.5/numTotalIndiv) / (percentMissingControls+.5/numTotalIndiv));
		if(differentialMissingOR < 1) differentialMissingOR = 1 / differentialMissingOR;
		if(differentialMissingOR > 30) differentialMissingOR = 30;
		
	}else{
		differentialMissingPval = 2.0;
		differentialMissingOR = -1.0;
	}
	
	
	return true;
}

/* 
 * Perform test of HW equilibrium.  
 * 
 * Uses exact test.  Calculates separately for Cases and Controls.
 * Note that this method is split into three completely separate tests.  This
 * could be pulled out into one test with parameters, but it would be a lot
 * of parameters.  
 */
bool PopStats::hwEquilibruim(int snp){
	
	if(pullData(snp)){
		numCategories(snp);
		minorAlleleFreq(snp); 
	}
	
	if(casesNonMissing > 0){ //if(numTotalCases > 0){
		
		double majAlleleFreqCase = 1.0 - minorAlleleFreqCases;
		
		#if DEBUG_HWEQUIL
			cout << "DEBUG: snp " << snp << " majAlleleFreqCase " << majAlleleFreqCase << endl;
		#endif
		
		
		caseExpTotal[0] = majAlleleFreqCase*majAlleleFreqCase*static_cast<double>(casesNonMissing);
		caseExpTotal[1] = 2*majAlleleFreqCase*minorAlleleFreqCases*static_cast<double>(casesNonMissing);
		caseExpTotal[2] = minorAlleleFreqCases*minorAlleleFreqCases*static_cast<double>(casesNonMissing);
		
		if((caseExpTotal[0] > 0.0) && (caseExpTotal[1] > 0.0) && (caseExpTotal[2] > 0.0)){
			caseChiSquare = (pow((static_cast<double>(ppCases) - caseExpTotal[0]), 2)/caseExpTotal[0]) +
							(pow((static_cast<double>(pqCases) - caseExpTotal[1]), 2)/caseExpTotal[1]) +
							(pow((static_cast<double>(qqCases) - caseExpTotal[2]), 2)/caseExpTotal[2]);
			
			#if DEBUG_HWEQUIL
				cout << "DEBUG: snp " << snp << " caseChiSquare " << caseChiSquare << endl;
			#endif
			
			try{
				caseChiSquarePValue = Statistics::chi2prob(caseChiSquare, 1);
			}catch(InvalidChiSquareException){
				caseChiSquarePValue = 2.0;
			}
		}else{
			caseChiSquarePValue = 2.0;
		}
		caseExactTestPVal = exactTest(ppCases, pqCases, qqCases);
		
	}else{
		stringstream ss;
		ss << "HW equil snp " << snp << " has no total cases." << endl;
		Logger::Instance()->writeLine(ss.str());
	}
	
	if(cntrlsNonMissing > 0){ //if(numTotalCntrls > 0){
		
		double majAlleleFreqCntrl = 1.0 - minorAlleleFreqCntrls;
		cntrlExpTotal[0] = majAlleleFreqCntrl*majAlleleFreqCntrl*static_cast<double>(cntrlsNonMissing);
		cntrlExpTotal[1] = 2*majAlleleFreqCntrl*minorAlleleFreqCntrls*static_cast<double>(cntrlsNonMissing);
		cntrlExpTotal[2] = minorAlleleFreqCntrls*minorAlleleFreqCntrls*static_cast<double>(cntrlsNonMissing);
		
		if((cntrlExpTotal[0] > 0.0) && (cntrlExpTotal[1] > 0.0) && (cntrlExpTotal[2] > 0.0)){
			cntrlChiSquare = (pow((static_cast<double>(ppCntrls) - cntrlExpTotal[0]), 2)/cntrlExpTotal[0]) +
							(pow((static_cast<double>(pqCntrls) - cntrlExpTotal[1]), 2)/cntrlExpTotal[1]) +
							(pow((static_cast<double>(qqCntrls) - cntrlExpTotal[2]), 2)/cntrlExpTotal[2]);
			try{
				cntrlChiSquarePValue = Statistics::chi2prob(cntrlChiSquare, 1);
			}catch(InvalidChiSquareException){
				cntrlChiSquarePValue = 2.0;
			}
		}else{
			cntrlChiSquarePValue = 2.0;
		}
		cntrlExactTestPVal = exactTest(ppCntrls, pqCntrls, qqCntrls);
	}else{
		stringstream ss;
		ss << "HW equil snp " << snp << " has no total controls." << endl;
		//Logger::Instance()->writeLine(ss.str());
	}
	
	// Calculate combined HW.
	if(casesNonMissing + cntrlsNonMissing > 0){
		int numTotCmbd = casesNonMissing+cntrlsNonMissing;
		
		int numMajAlleles = 2*ppCases+pqCases+2*ppCntrls+pqCntrls;
		//int numMinAlleles = 2*qqCases+pqCases+2*qqCntrls+pqCntrls;
		double majAlleleFreq = static_cast<double>(numMajAlleles)/(2*numTotCmbd);
		
		#if DEBUG_HWEQUIL
			cout << "numMajAlleles " << numMajAlleles << endl;
			cout << "majAlleleFreq " << majAlleleFreq << endl;
		#endif
		
		cmbdExpTotal[0] = majAlleleFreq * majAlleleFreq*static_cast<double>(numTotCmbd);
		cmbdExpTotal[1] = 2*majAlleleFreq*(1-majAlleleFreq)*static_cast<double>(numTotCmbd);
		cmbdExpTotal[2] = (1-majAlleleFreq)*(1-majAlleleFreq)*static_cast<double>(numTotCmbd);
		
		#if DEBUG_HWEQUIL
			cout << "Expected: " << cmbdExpTotal[0] << " - " << cmbdExpTotal[1] << " - " << cmbdExpTotal[2] << endl;
		#endif
		
		if((cmbdExpTotal[0] > 0.0) && (cmbdExpTotal[1] > 0.0) && (cmbdExpTotal[2] > 0.0)){
			cmbdChiSquare = (pow((static_cast<double>(ppCntrls + ppCases) - cmbdExpTotal[0]), 2)/cmbdExpTotal[0]) +
							(pow((static_cast<double>(pqCntrls + pqCases) - cmbdExpTotal[1]), 2)/cmbdExpTotal[1]) +
							(pow((static_cast<double>(qqCntrls + qqCases) - cmbdExpTotal[2]), 2)/cmbdExpTotal[2]);
			
			#if DEBUG_HWEQUIL
				cout << "Combined chi square: " << cmbdChiSquare << endl;
			#endif
			
			try{
				cmbdChiSquarePValue = Statistics::chi2prob(cmbdChiSquare, 1);
			}catch(StatsException){
				cmbdChiSquarePValue = 2.0;
			}
		}else{
			stringstream ss;
			ss << "HW for SNP " << snp << ": unable to compute combined HW analysis." << endl;
			Logger::Instance()->writeLine(ss.str());
			
			cmbdChiSquarePValue = 2.0;
		}
		cmbdExactTestPVal = exactTest(ppCntrls+ppCases,pqCntrls+pqCases,qqCntrls+qqCases);
		
	}else{
		stringstream ss;
		ss << "Hardy-Weinberg Analysis for SNP " << snp << ": There is no Case/Control Hardy-Weinberg Analysis because there are no cases or controls" << endl;
		Logger::Instance()->writeLine(ss.str());
		
	}
	
	return true;
}

/*
 * Exact test of HW equilibrium.  
 * 
 * Written by Matt Steigart (?) and taken verbatim.  
 * 
 * Static function!
 */
double PopStats::exactTest(double numPP, double numPQ, double numQQ){
	
	double doubleNumA, doubleNumB;
	int numA, numB;
	double exactValue;
	
	if(numPP + numPQ + numQQ > 0.0){
		
		doubleNumA = (2*numPP) + numPQ;
		doubleNumB = (2*numQQ) + numPQ;
		
		numA = static_cast<int>(floor(doubleNumA + 0.50));
		numB = static_cast<int>(floor(doubleNumB + 0.50));
		
		bool odd = true;
		if((numA % 2) == 0)
			odd = false;
			
			
		double freqA, freqB;
		freqA = static_cast<double>(numA)/(2*(numPP + numPQ + numQQ));
		freqB = static_cast<double>(numB)/(2*(numPP + numPQ + numQQ));
		
		int initNumAA, initNumAB, initNumBB;
		initNumAB = static_cast<int>((2*freqA*freqB)*(numPP + numPQ + numQQ));
		if((odd) && ((initNumAB % 2) == 0)){
			initNumAB += 1;
		}
		else if(!odd && ((initNumAB % 2) != 0)){
			initNumAB += 1;
		}
		initNumAA = static_cast<int>((numA - initNumAB)/2);
		initNumBB = static_cast<int>((numB - initNumAB)/2);
		
		int numAA, numAB, numBB;
		numAA = initNumAA;
		numAB = initNumAB;
		numBB = initNumBB;
		
		deque<double> probNAB;
		probNAB.resize(0);
		probNAB.push_back(100.0);
		deque<double> numHeterozygotes;
		numHeterozygotes.resize(0);
		numHeterozygotes.push_back(numAB);
		
		do{
			double ProbNMinusTwoAB;
			ProbNMinusTwoAB = probNAB.back()*
				((static_cast<double>(numAB)*(static_cast<double>(numAB) - 1))/
				(4*(static_cast<double>(numAA) + 1)*(static_cast<double>(numBB) + 1)));
			probNAB.push_back(ProbNMinusTwoAB);
			numAB -= 2;
			numAA++;
			numBB++;
			numHeterozygotes.push_back(numAB);
		}while((numAB - 2) >= 0);
		
		numAA = initNumAA;
		numAB = initNumAB;
		numBB = initNumBB;
		
		do{
			double ProbNPlusTwoAB;
			ProbNPlusTwoAB = probNAB.front()*
				((4*static_cast<double>(numAA)*static_cast<double>(numBB))/
				((static_cast<double>(numAB) + 2)*(static_cast<double>(numAB) + 1)));
			probNAB.push_front(ProbNPlusTwoAB);
			numAB += 2;
			numAA--;
			numBB--;
			numHeterozygotes.push_front(numAB);
		}while(numAA > 0);
		
		double sumWeights = 0.0;
		int iweight;
		for(iweight = 0; iweight < static_cast<int>(probNAB.size()); iweight++){
			sumWeights += probNAB[iweight];
		}
		
		for(iweight = 0; iweight < static_cast<int>(probNAB.size()); iweight++){
			probNAB[iweight] = probNAB[iweight]/sumWeights;
		}
		
		int jweight;
		for(jweight = 0; jweight < static_cast<int>(probNAB.size()); jweight++){
			if(equal(numHeterozygotes[jweight] , numPQ)){
				break;
			}
		}
		
		exactValue = 0.0;
		for(iweight = 0; iweight < static_cast<int>(probNAB.size()); iweight++){
			if(probNAB[jweight] >= probNAB[iweight]){
				exactValue += probNAB[iweight];
			}
		}
		return exactValue;
	}else{
		return 2.0;
	}
	
}
