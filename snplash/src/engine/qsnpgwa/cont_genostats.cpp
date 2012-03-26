//      cont_genostats.cpp
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

#include "cont_genostats.hh"
#include "../linalg/specialfunctions.h"

ContGenoStats::ContGenoStats(DataAccess *d){
	data = d;
	ranMeans = false;
}

/**
 * 
 * Perform all genotypic association operations.
 * 
 * Makes heavy use of the linear regression and ANCOVA engines.
 * 
 * Calculates
 * 	* 2 df test
 * 	* dom, add, rec test
 * 	* means and variance for dom, rec groupings
 * 	* lof test
 * 
 * @author Richard T. Guy
 */
void ContGenoStats::prepGenoStatsForOutput(int snp, ContGenoStatsResults &results){
	covariateAdjust(snp);
	computeMeanAndSD(snp, results);
	computeLinRegStats(snp, results);
	computeLackOfFit(snp, results);
}

/**
 * Compute the mean and standard deviation statistics for the following groups of individuals:
 * 
 * AA,Aa,aa,AA/Aa,Aa/aa
 * 
 */
void ContGenoStats::computeMeanAndSD(int snp, ContGenoStatsResults &results){
	
	vector<double> phenAA,phenAa,phenaa;
	double cumAA, cumAa, cumaa;
	double cumAAsq, cumAasq, cumaasq;
	double numAA, numAa, numaa;
	cumAA = cumAa = cumaa = 0;
	cumAAsq = cumAasq = cumaasq = 0;
	numAA = numAa = numaa = 0;
	double t;
	
	for(int i = 0; i < data->pheno_size(); i++){
		if(data->get_phenotype(i) != 0){
			t = data->get_phenotype(i);
			switch(data->get_data(i)->at(snp)){
				case 1:
					
					phenAA.push_back(t);
					cumAA+= t;
					cumAAsq+= t*t;
					numAA++;
				break;
				case 2:
					
					phenAa.push_back(t);
					cumAa+= t;
					cumAasq+= t*t;
					numAa++;
				break;
				case 3:
					phenAa.push_back(t);
					cumAa+= t;
					cumAasq += t*t;
					numAa++;
				break;
				case 4:
					phenaa.push_back(t);
					cumaa += t;
					cumaasq += t*t;
					numaa++;
				break;
				default:
					
				break;		
			}			
		}
	}
	
	/* AA */
	if(numAA >= 2){
		results.meanAA = cumAA / numAA;
		results.sdAA = pow( numAA/(numAA-1) * (cumAAsq/(numAA) - pow(cumAA/ (numAA),2)) , 0.5); 
	}else{
		results.meanAA = -999;
		results.sdAA = -1;	
	}

	/* Aa */
	if(numAa >= 2){
		results.meanAa = cumAa / numAa;
		results.sdAa = pow( numAa/(numAa-1) * (cumAasq/numAa - pow(cumAa / numAa,2)), 0.5); 
	}else{
		results.meanAa = -999;
		results.sdAa = -1;	
	}
	
	/* aa */
	if(numaa >= 2){
		results.meanaa = cumaa / numaa;
		results.sdaa = pow( numaa/(numaa-1) * (cumaasq/numaa - pow(results.meanaa,2)), 0.5);
	}else{
		results.meanaa = -999;
		results.sdaa = -1;	
	}
	
	/* AA Aa */
	if(numAa + numAA >= 2){
		results.meanAA_Aa = (cumAa + cumAA) / (numAa + numAA);
		results.sdAA_Aa = pow( (numAA+numAa)/(numAA+numAa-1) * ((cumAAsq + cumAasq)/(numAa+numAA) - pow(results.meanAA_Aa,2)), 0.5); 
	}else{
		results.meanAA_Aa = -999;
		results.sdAA_Aa = -1;	
	}

	/* Aa aa */
	if(numaa + numAa >= 2){
		results.meanAa_aa = (cumAa + cumaa) / (numAa + numaa);
		results.sdAa_aa = pow( (numaa+numAa)/(numaa+numAa-1) * ((cumAasq + cumaasq)/(numAa+numaa) - pow(results.meanAa_aa,2)), 0.5); 
	}else{
		results.meanAa_aa = -999;
		results.sdAa_aa = -1;	
	}

	/* class variable set*/
	if(numAA+numAa+numaa > 0){
		responseMean = ( cumAA+cumAa+cumaa ) / (numAA+numAa+numaa);
	}else{
		responseMean = -1;
	}

	ranMeans = true;
}

/**
 * Use the linear regression engine to compute statistics
 * 	  beta
 *    SE
 * 	  p-val 
 */
void ContGenoStats::computeLinRegStats(int snp, ContGenoStatsResults &results){

	if(!ranMeans) throw QSnpgwaException();
	
	if(residuals.size() == 0){
		blankLinRegStats(results);
		return;
	}
	vector<double> add_betas, dom_betas, rec_betas,  lof_betas; // return vals from lr engine.
	vector<double> dom, add, rec, ones;

	vector<vector<double> > add_in, dom_in, rec_in, lof_in;


	/* Prep genotypes */
	for(int i=0; i<data->pheno_size(); i++ ){
		// push onto stacks depending on the case/cntrl status.
		if(data->get_data(i)->at(snp) != 0){
			ones.push_back(1);
			switch(data->get_data(i)->at(snp)){
				case 1:
					add.push_back(-1);
					dom.push_back(0);
					rec.push_back(0);
				break;
				case 2:
					add.push_back(0);
					dom.push_back(1);
					rec.push_back(0);
				break;
				case 3:
					add.push_back(0);
					dom.push_back(1);
					rec.push_back(0);
				break;
				case 4:
					add.push_back(1);
					dom.push_back(1);
					rec.push_back(1);
				break;
				default:
				// won't happen.
				break;
			}
		}
	}
	
	add_in.push_back(ones);
	add_in.push_back(add);
	
	dom_in.push_back(ones);	
	dom_in.push_back(dom);
	
	rec_in.push_back(ones);
	rec_in.push_back(rec);
	
	/*
	 * F = SSR / p / ( SSE / (n - p - 1) ) for n obs and p vars not including intercept. 
	 */
	#if CONT_GENOSTATS_STATS
		cout << "Debug for qsnpgwa: " << snp << endl;
	#endif
	statisticsOutput stats = computeSingleStats(add_in, residuals, meanResidual, "Additive test", snp+1);
	results.add_pval = stats.pVal;
	results.add_beta = stats.beta;
	results.add_se 	 = stats.se;
	
	#if CONT_GENOSTATS_STATS
		cout << "Beta: " << setprecision(9) << stats.beta << endl;
		cout << "Pval: " << stats.pVal << endl;
		cout << "SE:   " << stats.se << endl;
	#endif
	
	stats = computeSingleStats(dom_in, residuals, meanResidual, "Dominant test", snp+1);
	results.dom_pval = stats.pVal;
	results.dom_beta = stats.beta;
	results.dom_se 	 = stats.se;
	
	
	stats = computeSingleStats(rec_in, residuals, meanResidual, "Recessive test", snp+1);
	results.rec_pval = stats.pVal;
	results.rec_beta = stats.beta;
	results.rec_se 	 = stats.se;
	
	// Compute the two degree of freedom test.
	double means[4];
	means[0] = results.meanAA; means[1] = results.meanAa; means[2] = results.meanaa;
	means[3] = responseMean;
	if(means[0] == -999 || means[1] == -999 || means[2] == -999 || means[3] == -999){
		results.twodegfree_pval = 2.0;
	}else{
		LinearRegression lr;
		results.twodegfree_pval = lr.anova(residuals, add, means);
	}
	
}

/**
 * Compute the lack of fit test statistic.
 *
 * @param snp The SNP index.
 * @param results The stats results object that will hold the response.
 */
void ContGenoStats::computeLackOfFit(int snp, ContGenoStatsResults &results){
	
	vector<double> lofScore, ones;
	vector<double> betasLOF;

	/* Prep genotypes */
	for(int i=0; i<data->pheno_size(); i++ ){
		// push onto stacks depending on the case/cntrl status.
		if(data->get_data(i)->at(snp) != 0){
			ones.push_back(1.0);
			switch(data->get_data(i)->at(snp)){
				case 1:
					lofScore.push_back(1);
					
				break;
				case 2:
					lofScore.push_back(-2);
				break;
				case 3:
					lofScore.push_back(-2);
				break;
				case 4:
					lofScore.push_back(1);
				break;
				default:
				// won't happen.
				break;
			}
		}
	}

	LinearRegression lr;
	vector<vector<double> > in;
	in.push_back(ones);
	in.push_back(lofScore);

	statisticsOutput stats = computeSingleStats(in, residuals, meanResidual, "Lack of fit test", snp+1);
	results.lof_pval = stats.pVal;
	
	
}

/**
 * Compute a single set of statistics including p_val, beta, se.
 * 
 * NOTE: The pvalues ect will be reported based on the last beta.
 * Order your variables accordingly.
 * 
 *	@param in The input to LR engine
 *	@param response The response variables.
 *	@param meanResidual of the response variable
 *	@return Struct with results for last beta.
 */
ContGenoStats::statisticsOutput ContGenoStats::computeSingleStats(const vector<vector<double> > &in, 
			const vector<double> &response, double meanResidual, string failMessage, int currentSNP)
{
	
	statisticsOutput ret;
	LinearRegression lr;
	try {
		
		vector<double> betas = lr.leastSquares(in, response);
		
		#if CONT_GENOSTATS_STATS
		for (unsigned int idb=0; idb < betas.size(); idb++)
			cout << betas[idb] << " ";
		cout << endl;
		#endif
		
		double sse, ssr, f;
		lr.sumSquaredStats(in, betas, response, meanResidual, sse, ssr);
		
		f = ssr / (sse / (response.size() - 2 ));
		
		ret.beta = betas.at(betas.size()-1);
		double lxx = sumSquaresTotal(in.at(in.size()-1)); // Corrected sumsquares for x.
		ret.se = sqrt(sse / lxx) / (response.size() - 2 );
		
		#if CONT_GENOSTATS_STATS
		cout << "Sum squared stats: " << endl;
		cout << "meanResidual: " << meanResidual << endl;
		cout << "sse:          " << sse << endl;
		cout << "ssr:          " << ssr << endl;
		cout << "f:            " << f << endl;
		cout << "deg freedom:  " << 1 << " " << response.size() - 2 << endl;
		#endif
		
		ret.pVal = 1.0 - alglib::fdistribution(1, response.size() - 2, f);
		//ret.pVal = Statistics::fdist(f, 1, response.size() - 2);
	}catch(ConditionNumberEx){
		stringstream ss;
		ss << "Error on snp " << currentSNP << " " << failMessage << " condition number exception in single stats." << endl;
		Logger::Instance()->writeLine(ss.str());

		ret.pVal = 2.000;
		ret.se = -1;
		ret.beta = -1;
	}catch(LinearRegressionException){
		stringstream ss;
		ss << "Error on snp " << currentSNP << " " << failMessage << " singular matrix exception in single stats." << endl;
		Logger::Instance()->writeLine(ss.str());

		ret.pVal = 2.000;
		ret.se = -1;
		ret.beta = -1;
	}catch(...){
		stringstream ss;
		ss << "Error on snp " << currentSNP << " " << failMessage << " unknown exception in single stats." << endl;
		Logger::Instance()->writeLine(ss.str());
		ret.pVal = 2.000;
		ret.se = -1;
		ret.beta = -1;
	}

	return ret;

}

/**
 * Default all pvals, se, and betas.
 * 
 * @param results A Genostats results struct that we will turn to 
 * all default values.
 */
void ContGenoStats::blankLinRegStats(ContGenoStatsResults &results){
	
	results.twodegfree_pval = 2;
	
	results.add_pval = 2;
	results.add_beta = -1;
	results.add_se = -1;
	
	results.dom_pval = 2;
	results.dom_beta = -1;
	results.dom_se = -1;
	
	results.rec_pval = 2;
	results.rec_beta = -1;
	results.rec_se = -1;
	
	results.lof_pval = 2.0;
	
}

/**
 * Return the total sum of squares.
 * Equal to sum of norm of difference to mean.
 * 
 */
double ContGenoStats::sumSquaresTotal(vector<double> response){
	
	
	double sum, sumsq;
	sum = sumsq = 0.0;
	double sst = 0;
	for(unsigned int i=0;i<response.size();i++){
		sum += response.at(i);
		sumsq += pow(response.at(i),2);
	}
	sst = (sumsq - pow(sum , 2)/response.size()) / response.size();
	return sst;
}

/**
 * Fill the residuals column.
 * 
 * All other tests are designed to ignore covariates because they are preadjusted away here.
 * 
 * NOTE: The residuals are stored for all individuals in the data set for
 * whom the current SNP's value is not 0.  Therefore, it can be used directly
 * rather than recomputed.
 * 
 * @param snp The SNP we are operating on.
 */
void ContGenoStats::covariateAdjust(int snp){
	
	meanResidual = 0;
	vector<double> ones, phen_vec;
	vector<vector<double> > cov;
	/* Prep covariate matrix */
	vector<double> *tmp = data->get_covariates(0);
	for(unsigned int i=0;i<tmp->size();i++){
		vector<double> a;
		cov.push_back(a);
	}
	
	/* Prep genotypes and fill cov matrix */
	for(int i=0; i<data->pheno_size(); i++ ){
		if(data->get_data(i)->at(snp) != 0){
			phen_vec.push_back(data->get_phenotype(i));
			ones.push_back(1.0);
			vector<double> *t = data->get_covariates(i);
			for(unsigned int j=0;j<t->size();j++){
				cov.at(j).push_back(t->at(j));
			}
		}
	}
	
	cov.push_back(ones);
	
	#if CONT_GENOSTATS_STATS
	
		cout << "Covariate size for computing residuals: " << cov.size() << endl;
	
	#endif
	
	LinearRegression lr;
	
	try {
		vector<double> betas = lr.leastSquares(cov, phen_vec);
		
		residuals = lr.residuals(cov, betas, phen_vec);
	}catch(...){
		cout << "Error in residuals." << endl;
		residuals.clear();
		return;
	}
	
	
	for(unsigned int i=0;i < residuals.size(); ++i)
		meanResidual += residuals.at(i);
	meanResidual /= residuals.size();

}
