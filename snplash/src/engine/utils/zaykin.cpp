//      zaykin.cpp
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

#include "zaykin.hh"

Zaykin::Zaykin(EngineParamReader *p){
	params = p;
	keepThresh = 20;
	data = 0;
}

Zaykin::~Zaykin(){
}

void Zaykin::setPhenotype(const vector<double> &p, const vector<vector<double> > &c){
	phenotype = p;
	cov = c;
}

void Zaykin::setErrorInformation(int snp, int numCovariates, DataAccess *data){
	
	this->snp = snp;
	this->numCovariates = numCovariates;
	this->data = data;
	
}

/**
 * Retrieve data and build internal representation.
 *
 * Given a series of haplogenotypes, for a given individual, calculate
 * the prob of each haplotype involved.
 *
 * The class variable matrix will be correctly filled with only those haplotypes that
 * have greater than 0.05 global frequency.  It will have been rescaled.
 * 
 * RTG updated Oct 14, 2010 to change the thresholding from 0.05% global expression
 * to showing up in at least n chromosomes (2m chr for m people).  Default is 20.
 *
 * @param input vector of EMPersonalProbsResults note that the person index is 0 based.
 * @param numHaps total number of haplotypes.
 * @param reweight If true, reweight the kept haplotypes.
 */
void Zaykin::setup(const vector<EMPersonalProbsResults> &input, int numHaps, bool reweight){

	if(phenotype.size() == 0){
		RunOrderException e;
		e.message = "Zaykin: called setup without setting a phenotype.";
		throw e;
	}

	vector<vector<double> > matrixT; // temp used internally.
	vector<int> used(numHaps, 0);
	
	matrixT = vecops::getDblVec(phenotype.size() , numHaps);

	for(unsigned int i=0; i < input.size(); ++i){
		double t2 = input.at(i).prob;
		if(isnan(t2)) t2 = 0;
		// compute the actual haplotype number
		// place half the percentage in that element in this individual.
		int ti = input.at(i).leftHap;
		matrixT.at(input.at(i).personId).at(ti) += t2 / 2;
		used.at(ti)++;
		ti = input.at(i).rightHap;
		matrixT.at(input.at(i).personId).at(ti) += t2 / 2;
		used.at(ti)++;
	}

	vector<unsigned int> keepIndices;
	usedHaplotype.clear();
	usedHaplotype.reserve(used.size());
	for(unsigned int i=0; i < used.size(); i++){
		
		if(used.at(i) > keepThresh){
			
			usedHaplotype.push_back(true);
			keepIndices.push_back(i);
		}else{
			
			usedHaplotype.push_back(false);
		}
	}

	matrix = vecops::getDblVec(phenotype.size() , keepIndices.size());
	for(unsigned int i=0;i < phenotype.size(); i++){
		
		for(unsigned int j=0;j < keepIndices.size(); j++){
			matrix.at(i).at(j) = matrixT.at(i).at(keepIndices.at(j));
		}
		if(reweight){
			double cumSum = vecops::vecCumSum(matrix.at(i));
			vecops::vecDiv(matrix.at(i) , cumSum);
		}
	}
}

/**
 * Prepare haplotypes.
 */
void Zaykin::prepHaplotypes(vector<double> &ones, vector<vector<double> > &haps){

	ones.clear();
	haps.clear();
	// prep haps
	haps = vecops::getDblVec(matrix.at(0).size() , phenotype.size());

	#if DEBUG_ZAY_PROGRESS
		cout << "Zaykin start genotype and cov prep" << endl;
	#endif

	/* Prep haplotypes */
	for(unsigned int i=0; i<phenotype.size(); i++ ){

		ones.push_back(1.0);

		for(unsigned int j=0; j < matrix.at(i).size(); j++){
			if(matrix.at(i).at(j) != matrix.at(i).at(j)){
				haps.at(j).at(i) = 0;
			}else{
				haps.at(j).at(i) = matrix.at(i).at(j);
			}
		}
	}
}

/**
 * Run the global analysis.  Put all data and covariates in a logistic regression model.
 *
 * Compute likelihood ratio statistic.
 * This uses a likelihood statistic where n-1 haplotypes are tested.
 * 
 * @return pvalue.
 */
ZaykinGlobalStatsResults Zaykin::runGlobal(){

	ZaykinGlobalStatsResults stats;

	vector<vector<double> > haps;
	vector<double> ones;

	prepHaplotypes(ones, haps);

	#if DEBUG_ZAY_PROGRESS
		cout << "Zaykin start LR portion" << endl;
	#endif
	LogisticRegression with(params->getRegressionConditionNumberThreshold()), without(params->getRegressionConditionNumberThreshold());
	/*
	 * Run without the haplotypes:
	 */
	vector<vector<double> > inv_infmatrixWithOut;
	vector<double> betasWithOut;

	vector<vector<double> > inWithout;
	inWithout = cov;
	inWithout.push_back(ones);

	inv_infmatrixWithOut = vecops::getDblVec(inWithout.size() , inWithout.size());

	int retry = 0;
	double startVal = 0;  // value to start betas with.

	while(retry < 3){
		try{
			betasWithOut = without.newtonRaphson(inWithout, phenotype, inv_infmatrixWithOut, startVal);
			break;
		}catch(NewtonRaphsonFailureEx){
			handleException(stats, startVal, retry, "Unable to compute reduced model: Newton-Raphson setup failure.");
		}catch(NewtonRaphsonIterationEx){
			handleException(stats, startVal, retry, "Unable to compute reduced model: max iterations hit.");
		}catch(SingularMatrixEx){
			handleException(stats, startVal, retry, "Unable to compute reduced model: information matrix was singular.");
		}catch(ConditionNumberEx err){
			
			LogisticRegression lr;
			int separableVariable = lr.dataIsSeparable(inWithout, phenotype);
			
			string message;
			if (separableVariable < 0){
				// Error: poor conditioning.
				stringstream ss;
				ss << "Unable to compute reduced model: Poor conditioning in information matrix. ";
				ss << "Condition number (1-norm) is " << err.conditionNumber;
				message = ss.str();
			}else{
				stringstream ss;
				ss << "Unable to compute reduced model: Separable data matrix.";
				ss << "Condition number (1-norm) is " << err.conditionNumber;
				message = ss.str();
			}
			handleException(stats, startVal, retry, message);
			if (retry >= 3) return stats;	
		}catch(ADTException e){
			// This one is generic.
			string message = "Unable to compute reduced model: Newton-Raphson error.";
			handleException(stats, startVal, retry, message);
			if (retry >= 3) return stats;
		}catch(alglib::ap_error err){
			stringstream ss;
			ss << "Unable to compute reduced model due to linalg exception: " << err.msg;
			handleException(stats, startVal, retry, ss.str());
			if (retry >= 3) return stats;
		}
	}

	/*
	 * Run with the haplotypes:
	 */
	vector<vector<double> > inv_infmatrixWith, inWith;
	vector<double> betasWith;

	inWith = inWithout;
	for (unsigned int i=0; i < haps.size()-1; i++){ // NOTE: Don't push the very last haplotype.
		inWith.push_back(haps.at(i));
	}

	inv_infmatrixWith = vecops::getDblVec(inWith.size() , inWith.size());

	retry = 0;
	startVal = 0;  // value to start betas with.

	while(retry < 3){
		try{
			betasWith = with.newtonRaphson(inWith, phenotype, inv_infmatrixWith, startVal);
			break;
		}catch(NewtonRaphsonFailureEx){
			handleException(stats, startVal, retry, "Unable to compute full model: Newton-Raphson setup failure.");
		}catch(NewtonRaphsonIterationEx){
			handleException(stats, startVal, retry, "Unable to compute full model: max iterations hit.");
		}catch(SingularMatrixEx){
			handleException(stats, startVal, retry, "Unable to compute full model: information matrix was singular.");
		}catch(ConditionNumberEx err){
			
			LogisticRegression lr;
			int separableVariable = lr.dataIsSeparable(inWith, phenotype);
			
			stringstream ss;
			if (separableVariable < 0){
				// Error: poor conditioning.
				ss << "Unable to compute reduced model: Poor conditioning in information matrix. ";
				ss << "Condition number (1-norm) is " << err.conditionNumber;
			}else{
				ss << "Unable to compute reduced model: Separable data matrix.";
				ss << "Condition number (1-norm) is " << err.conditionNumber;
			}
			string message = ss.str();
			handleException(stats, startVal, retry, message);
			if (retry >= 3) return stats;	
		}catch(ADTException e){
			string message = "Unable to compute full model: Newton-Raphson error.";
			handleException(stats, startVal, retry, message);
			if (retry >= 3) return stats;
		}catch(alglib::ap_error err){
			stringstream ss;
			ss << "Unable to compute full model due to linalg exception: " << err.msg;
			handleException(stats, startVal, retry,ss.str());
			if (retry >= 3) return stats;
		}
	}

	double likeRatio =  with.likelihoodRatio(betasWithOut, inWithout, betasWith, inWith, phenotype);
	try{
		stats.pvalue = Statistics::chi2prob(likeRatio, betasWith.size() - betasWithOut.size());
		stats.testStat = likeRatio;
		stats.degFreedom = betasWith.size() - betasWithOut.size();
	}catch(...){

		stringstream ss;
		ss << "Zaykin's method: unable to compute chi square: " << likeRatio << " " << betasWith.size() - betasWithOut.size() << endl;
		Logger::Instance()->writeLine(ss.str());

		stats.fillDefault();
		return stats;
	}
	return stats;
}

/**
 * Run each haplotype using Zaykin's method.
 *
 * Assumes that usedHaplotype was set.
 *
 */
vector<ZaykinStatsInfo>  Zaykin::runAllHaplotypes(){
	
	vector<ZaykinStatsInfo> ret;
	ZaykinStatsInfo stats;
	
	vector<double> ones;
	vector<vector<double> > haps;

	prepHaplotypes(ones, haps); // set up haps vector and ones vector.

	// Run each haplotype alone.
	int actualHaplotypeCntr = 0; // Counts through the haplotypes used.
	for(unsigned int i=0; i < usedHaplotype.size(); i++){ // for each haplotype.

		if(usedHaplotype.at(i)){

			stats = runLikelihoodRatio(haps.at(actualHaplotypeCntr), ones);
			
			actualHaplotypeCntr++;
		}else{
			stats.fillDefault();
		}
		ret.push_back(stats);
	}
	return ret;
}

/**
 * Run a single likelihood ratio test.  Assumes covariates and phenotype are set as global variables.
 * The test is performed on the haps vector.
 *
 * @param haps Haplotype variable.  H_0 is that this vector is independant.
 * @param ones Set of ones.  Precomputed for speed.
 * @return ZaykingStatsInfo should contain all information for this likelihood ratio test.
 */
ZaykinStatsInfo Zaykin::runLikelihoodRatio(const vector<double> &haps, const vector<double> &ones){

	ZaykinStatsInfo stats;
	vector<vector<double> > testVecWithout = cov;
	vector<vector<double> > testVecWith; // wait to fill this one.
	testVecWithout.push_back(ones);

	LogisticRegression lrWith, lrWithout;
	vector<vector<double> > inv_infmatrixWithOut, inv_infmatrixWith;
	vector<double> betasWith, betasWithOut;

	inv_infmatrixWithOut = vecops::getDblVec(testVecWithout.size() , testVecWithout.size());

	int retry = 0;
	double startVal = 0;  // value to start betas with.

	while(retry < 3){
		try{
			betasWithOut = lrWithout.newtonRaphson(testVecWithout, phenotype, inv_infmatrixWithOut, startVal);
			break;
		}catch(NewtonRaphsonFailureEx){
			handleException(stats, startVal, retry, "Unable to compute reduced model in single haplotype test: Newton-Raphson setup failure.");
		}catch(NewtonRaphsonIterationEx){
			handleException(stats, startVal, retry, "Unable to compute reduced model in single haplotype test: max iterations hit.");
		}catch(SingularMatrixEx){
			handleException(stats, startVal, retry, "Unable to compute reduced model in single haplotype test: information matrix was singular.");
		}catch(ConditionNumberEx){
			
			LogisticRegression lr;
			int separableVariable = lr.dataIsSeparable(testVecWithout, phenotype);
			
			string message;
			if (separableVariable < 0){
				// Error: poor conditioning.
				message = "Unable to compute reduced model in single haplotype test: Poor conditioning in information matrix.";
			}else{
				message = "Unable to compute reduced model in single haplotype test: Separable data matrix.";
			}
			handleException(stats, startVal, retry, message);
			if (retry >= 3) return stats;	
		}catch(alglib::ap_error err){
			stringstream ss;
			ss << "Unable to compute reduced model  single haplotype test due to linalg exception: " << err.msg;
			handleException(stats, startVal, retry, ss.str());
			if (retry >= 3) return stats;
		}
	}

	retry = 0;
	startVal = 0;  // value to start betas with.

	testVecWith = testVecWithout;
	testVecWith.push_back(haps);

	inv_infmatrixWith = vecops::getDblVec(testVecWith.size() , testVecWith.size());

	while(retry < 3){
		try{
			betasWith = lrWith.newtonRaphson(testVecWith, phenotype, inv_infmatrixWith, startVal);
			break;
		}catch(NewtonRaphsonFailureEx){
			handleException(stats, startVal, retry, "Unable to compute full model in single haplotype test: Newton-Raphson setup failure.");
		}catch(NewtonRaphsonIterationEx){
			handleException(stats, startVal, retry, "Unable to compute full model in single haplotype test: max iterations hit.");
		}catch(SingularMatrixEx){
			handleException(stats, startVal, retry, "Unable to compute full model in single haplotype test: information matrix was singular.");
		}catch(ConditionNumberEx){
			
			LogisticRegression lr;
			int separableVariable = lr.dataIsSeparable(testVecWith, phenotype);
			
			string message;
			if (separableVariable < 0){
				// Error: poor conditioning.
				message = "Unable to compute full model in single haplotype test: Poor conditioning in information matrix.";
			}else{
				message = "Unable to compute full model in single haplotype test: Separable data matrix.";
			}
			handleException(stats, startVal, retry, message);
			if (retry >= 3) return stats;	
		}catch(alglib::ap_error err){
			stringstream ss;
			ss << "Unable to compute full model single haplotype test due to linalg exception: " << err.msg;
			handleException(stats, startVal, retry, ss.str());
			if (retry >= 3) return stats;
		}
	}

	stats.chiSqStat = lrWith.likelihoodRatio(betasWithOut, testVecWithout, betasWith, testVecWith, phenotype);
	stats.degFree = betasWith.size() - betasWithOut.size();
	double beta = betasWith.at(betasWith.size() - 1);
	double stderr = sqrt(inv_infmatrixWith.at(inv_infmatrixWith.size() - 1).at(inv_infmatrixWith.size() - 1));
	stats.OR = exp(beta);
	stats.LCI = exp(beta - 1.96*stderr);
	stats.UCI = exp(beta + 1.96*stderr);
	return stats;

}

/**
 * Run all individual tests to test dom, add, rec relationship with each phenotype that had enough elements.
 *
 * @return vector of pValues, one per person.
 *
 * @bug This method not done yet.
 */
vector<double> Zaykin::runAllIndiv(vector<double> &testStats, vector<int> &degsFreedom){

	testStats.clear(); degsFreedom.clear();
	vector<double> pValues;



	return pValues;
}


/**
 * Dump a vector to the screen for testing purposes only.
 *
 * @param m A double vector to be printed.
 */
void Zaykin::dump(vector<vector<double> > m){

	cout << "Matrix dump: " << endl;
	for(unsigned int i=0;i < m.size();i++){
		for(unsigned int j=0;j < m.at(i).size();j++){
			cout << m.at(i).at(j) << " ";
		}
		cout << endl;
	}

}

// recode a number using 1,2 to 0,1.
int Zaykin::recode12_01(int i){
	vector<int> t;
	int ret;
	while(i > 0){
		t.push_back(i % 2);
		i /= 2;
	}

	ret = 0;
	for(unsigned int ii=0;ii < t.size(); ii++){
		ret = ret*2 + t.at(ii)-1;
	}
	return ret;
}

/*
 * Test the Zaykin package.
 */
int Zaykin::test(){

	// Test 1: Input.
	vector<EMPersonalProbsResults> testV;

	EMPersonalProbsResults t;
	t.personId = 0;
	t.leftHap = 0;
	t.rightHap = 1;
	t.prob = .5;
	testV.push_back(t);

	t.personId = 0;
	t.leftHap = 1;
	t.rightHap = 0;
	t.prob = .5;
	testV.push_back(t);


	return 1;
}

void Zaykin::handleException(StatsFillable &stats, double &startVal, int &retry, const string &message){
	
	retry++;
	if(retry == 1) {startVal = 0.5;}
	else if(retry == 2){ startVal = -0.5;}
	else{
		retry = 3;
		if(retry == 1) {startVal = 0.5;}
		else if(retry == 2){ startVal = -0.5;}
		else{

			stringstream ss;
			string tempS;
			int tempP;
			string name1;
			if (data != 0){
				(*data).get_map_info(snp, tempS, name1, tempP);
				ss << "Zaykin for SNP " << name1 << ": " << message << endl;
			}else{
				ss << "Zaykin: " << message << endl;
			}

			Logger::Instance()->writeLine(ss.str());
	
			stats.fillDefault();
		}
	}
}
