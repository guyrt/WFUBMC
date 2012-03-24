//      linear_regression.cpp
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

#include "linear_regression.hh"
#include "../utils/exceptions.h"
#include "../utils/float_ops.hh"
#include "../linalg/linalg.h" // used for condition number.
#include "../linalg/alglibinternal.h"
#include "../linalg/blas.h"
#include "../linalg/specialfunctions.h"

LinearRegression::LinearRegression(){
}

/**
 * Compute the ordinary least squares estimation of the solution to a linear regression
 * 
 * @param data The regression input.  Each inner vector should be a single variable (so column maj order)
 * @param response The response variable.
 * @return betas The estimate of regression coefficients
 */
vector<double> LinearRegression::leastSquares(const vector<vector<double> > &data, const vector<double> &response){
	
	int sz = data.size();
	if(sz < 0) throw LinearRegressionException();
	
	int indivSz = response.size();
	
	vector<double> betas;
	
	alglib::real_1d_array vecResp;
	alglib::real_1d_array vecBetas;
	alglib::real_2d_array vecData;
	alglib::real_2d_array tempSquare;
	alglib::real_1d_array tempWork;
	
	vecResp.setlength(indivSz);
	vecBetas.setlength(sz);
	vecData.setlength(indivSz,sz);
	tempSquare.setlength(sz,sz);
	tempWork.setlength(indivSz+1);
	
	// Transfer into the matrices.
	for(unsigned int i=0;i<data.size();i++){
		for(unsigned int j=0;j < data.at(0).size();j++){
			vecData(j,i) = data.at(i).at(j);
		}
	}
	for(unsigned int i=0;i<response.size();i++)
		vecResp(i) = response.at(i);
	
	
	/// X'*X
	matrixmatrixmultiply(vecData,0,indivSz-1,0,sz-1,true,
						vecData,0,indivSz-1,0,sz-1,false,1.0,
						tempSquare,0,sz-1,0,sz-1,0.0,tempWork);
	
	alglib::matinvreport report;

	alglib::ae_int_t reportInfo ;
	rmatrixinverse(tempSquare, reportInfo, report);
	
	if(reportInfo != 1){
		throw LinearRegressionException();
	}
	if (report.r1 > LINEAR_REGRESSION_CONDITION_NUMBER_LIMIT){
		throw ConditionNumberEx(1.0/report.r1);
	}
		
	/// X'*beta
	matrixvectormultiply(vecData,0,indivSz-1,0,sz-1,true,
						vecResp,0,indivSz-1,1.0,
						tempWork,0,sz-1,0.0);
						
	/// tempSquare * tempWork
	matrixvectormultiply(tempSquare,0,sz-1,0,sz-1,false,
						tempWork,0,sz-1,1.0,
						vecBetas,0,sz-1,0.0);
	
	for(int i=0;i<sz;i++)
		betas.push_back(vecBetas(i));
	
	return betas;
}

/**
 * Compute the sum of squared errors for a given set of beta coefficients.
 * 
 * @param data The regression input.  Each inner vector is a single variable.
 * @param betas The betas to be compared.
 * @param response The response to compare against.
 * @param mean The mean of the response variable
 * @return sse The sum squared error
 * @return ssr The sum squared residual.
 */
void LinearRegression::sumSquaredStats(const vector<vector<double> > &data, const vector<double> &betas, 
				const vector<double> &response, double mean, double &sse, double &ssr)
{
	
	sse = 0;
	ssr = 0;
	double temp = 0;
	for(unsigned int i=0;i<response.size();i++){
		
		temp = 0;
		for(unsigned int j=0;j<betas.size();j++){
			temp += betas.at(j) * data.at(j).at(i);
		}
		sse += pow( response.at(i) - temp ,2);
		ssr += pow( mean - temp,2);
	}
}

/**
 * Compute the residuals of a given parameter set betas.
 * 
 * @param data vector of vectors.
 * @param betas vector of parameters
 * @param response vector of dependent variable values.
 */
vector<double> LinearRegression::residuals(const vector<vector<double> > &data, const vector<double> &betas, const vector<double> &response){
	
	vector<double> resid;

	double temp = 0;
	for(unsigned int i=0;i < response.size(); ++i){
		temp = 0;
		for(unsigned int j=0;j<betas.size();j++){
			temp += betas.at(j) * data.at(j).at(i);
		}
		resid.push_back(response.at(i) - temp);
	}
	return resid;
}

/**
 * Anova of the response variable.  In practice, this is called from
 * the genostats::compLinRegStats method using the calculated genotype
 * from the additive model.  That is why we expect -1,0,1 coding.
 * 
 * @param response Vector of response variables.
 * @param division Vector of SNP variables.  Expects -1,0,1 coding.
 * @param means[3] holds mean(-1), mean(0), mean(1), mean(tot) where mean is of response variable.
 */
double LinearRegression::anova(const vector<double> &response, const vector<double> &division, double means[4]){
	
	double partialAA = 0, partialAa = 0, partialaa = 0;
	double wss = 0, bss = 0;
	double stat;
	double retval = 0;
	int AA = 0, Aa = 0, aa = 0;
	
	for(unsigned int i=0;i<response.size();i++){
		
		if(equal(division.at(i) ,-1)){
			partialAA += pow( means[0] - response.at(i) ,2);
			AA++;
		}else if(equal(division.at(i) ,0)){
			partialAa += pow( means[1] - response.at(i) ,2);
			Aa++;
		}else if(equal(division.at(i), 1)){
			partialaa += pow( means[2] - response.at(i) ,2);
			aa++;
		}
	}
	
	bss = static_cast<double>(AA) * pow( means[0] - means[3] , 2) + Aa * pow(means[1]-means[3],2) + aa * pow(means[2]-means[3],2);
	wss = partialAA + partialAa + partialaa;
	
	stat = ( (division.size() - 3) / 2 ) * bss / wss;
	try{
		retval = alglib::fdistribution(2, response.size() - 3, stat);
	}catch(...){
		retval = 2.0;
	}
	return retval;
}

void LinearRegression::test(){
	/* Build a set of test data and get the betas.  Do they make sense? */
	
	vector<vector<double> > inMat;
	vector<double> m1;
	m1.push_back(0);
	m1.push_back(2);
	m1.push_back(1);
	m1.push_back(0);
	m1.push_back(-1);
	m1.push_back(0);
	
	vector<double> m2;
	m2.push_back(1);
	m2.push_back(1);
	m2.push_back(1);
	m2.push_back(-1);
	m2.push_back(0);
	m2.push_back(1);
	
	
	vector<double> m3;
	m3.push_back(0);
	m3.push_back(1);
	m3.push_back(0);
	m3.push_back(1);
	m3.push_back(0);
	m3.push_back(1);
		
	inMat.push_back(m1);
	inMat.push_back(m2);
	
//	vector<double> ans = leastSquares(inMat, m3); 
	
}

