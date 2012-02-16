//      lr.cpp
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

#include "lr.hh"

#include "../linalg/linalg.h" // used for condition number.
#include "../linalg/alglibinternal.h"
#include "../linalg/blas.h"

#include "../../logger/log.hh"

LogisticRegression::LogisticRegression(){
     condition_number_limit = 1e-12;
}
LogisticRegression::LogisticRegression(double conditionNumber){
    condition_number_limit = 1 / conditionNumber;
}
/**
 * Check whether any column in the data completely separates the response variable. 
 * If so, return the index.
 */
int LogisticRegression::dataIsSeparable(const vector<vector<double> > &data, const vector<double> &response){
    
    vector<double> covariance(data.size(), 0);
    vector<double> means(data.size(), 0);
    vector<double> variances(data.size(), 0);
    
    double varResponse = vecops::vecVariance(response);
    double meanResponse = vecops::vecCumSum(response);
    meanResponse /= response.size();
    
    for (unsigned int i=0; i < data.size(); i++){
        means[i] = vecops::vecCumSum(data[i]);
        variances[i] = vecops::vecVariance(data[i]);
    }
    vecops::vecDiv<double>(means, data[0].size());
    
    // Now compute E[XY]
    for (unsigned int i=0; i < data.size(); i++){
        for (unsigned int j=0; j < data[i].size(); j++){
            covariance[i] += data[i][j] * response[j];
        }
    }
    
    vecops::vecDiv<double>(covariance, data[0].size());
    
    for (unsigned int i=0; i < data.size(); i++){
        covariance[i] -= (means[i]*meanResponse);
        covariance[i] = abs(covariance[i]);
    }
    
    for (unsigned int i=0; i < data.size(); i++){
        covariance[i] /= ( sqrt(varResponse * variances[i] ));
    }
    
    double mx = -1.0;
    int mxLocation = -1;
    for (unsigned int i=0; i < covariance.size(); i++){
        if (covariance[i] > mx){
            mx = covariance[i];
            mxLocation = i;
        }
    }
    if (mx >LogisticRegression::SEPARABLE_THRESHOLD)
        return mxLocation;
    
    return -1;
}

/**
 * Compute p-value and SE.
 * 
 * P-value uses Z score.
 */
LRStats LogisticRegression::getSingleStats(const vector<double> &betas, const vector<vector<double> > &invInfMatrix, int index){
    
    LRStats ret;
    double beta = betas.at(index);
    ret.OR = exp(beta);
    
    double invInf = invInfMatrix[index][index];
    
    ret.invInf = sqrt(invInf);
    ret.LCI = exp(beta - 1.96*sqrt(invInf));
    ret.UCI = exp(beta + 1.96*sqrt(invInf));
    
    if(ret.OR != ret.OR){
        ret.OR = -1;
        ret.UCI = -1;
        ret.LCI = -1;
    }else if(ret.OR >= 10000.0){
        stringstream ss;
        ss << "Large odds ratio detected and rewritten: OR was " << ret.OR << "." << endl;
        Logger::Instance()->writeLine(ss.str());
        ret.OR = 9999.0;
        ret.UCI = 9999.0;
        ret.LCI = 9999.0;
    }
    
    double z;
    z = beta/sqrt(invInf);
    if(!equal(z , z )){
        ret.pVal = 2.0;
        ret.testStat = -1;
    }
    else{
        ret.pVal = Statistics::normalPValue(z);
        ret.testStat = z;
        //Computing normalPvalue sets pValue to 0 if |z|>12.7
        //That is a problem, so if it happens, recompute pValue based on 
        //relationship between chi-sq and normal distribution
        if(ret.pVal <= 0){
            try{
                ret.pVal = Statistics::chi2prob(z*z,1.0);
                
            }catch(StatsException){
                ret.pVal = 2;
                ret.testStat = -1;
            }
        }
    }
    
    return ret;
    
}

/*
 * Compute p-val using Wald test.
 * 
 * @param betas A vector of beta values
 * @param invInf Inverse of corresponding information matrix.
 * @return chiS The test statistic
 * @return p-value.
 */
double LogisticRegression::getStats(const vector<double> &betas, const vector<vector<double> > invInf, double &chiS){
    
    if(betas.size() != invInf.size())
        return 2.0;
    
    
    int sz = betas.size()-1;
    
    // Make vector and matrix.
    alglib::real_1d_array vBetas;
    alglib::real_1d_array vTemp;
    alglib::real_1d_array singleNum;
    alglib::real_2d_array vInvInf;
    
    vBetas.setlength(sz+1);
    vTemp.setlength(sz+1);
    singleNum.setlength(1);
    vInvInf.setlength(sz+1, sz+1);
    
    for(unsigned int i=0;i < betas.size();i++){
        vBetas(i) = betas.at(i);
        for(unsigned int j=0;j < betas.size() ; j++){
            vInvInf(i,j) = invInf.at(i).at(j);
        }
    }
    
    #if DEBUG_STATS
    cout << "Before:" << endl;
    cout << vInvInf(0,0) << " " << vInvInf(0,1) << endl;
    cout << vInvInf(1,0) << " " << vInvInf(1,1) << endl;
    #endif
    
    alglib::matinvreport report;
    alglib::ae_int_t reportInfo;
    rmatrixinverse(vInvInf, reportInfo, report);
    if(reportInfo != 1) return 2.0;
    #if DEBUG_STATS
    cout << "After:" << endl;
    cout << vInvInf(0,0) << " " << vInvInf(0,1) << endl;
    cout << vInvInf(1,0) << " " << vInvInf(1,1) << endl;
    #endif
    
    // Check condition number.
    if (report.r1 < condition_number_limit){
        throw ConditionNumberEx(1.0/report.r1);
    }
    
    // betas' * vInvInf * betas  (vInvInf has been inverted)
    matrixvectormultiply(vInvInf,0,sz,0,sz,false,
                        vBetas,0,sz,1.0,vTemp,0,sz,0.0);
    
    alglib::real_2d_array vBetasTemp;
    vBetasTemp.setlength(sz+1,1);
    for(int i=0;i<=sz;i++)
        vBetasTemp(i,0) = vBetas(i);
    
    matrixvectormultiply(vBetasTemp,0,sz,0,0,true,
                        vTemp,0,sz,1.0,singleNum,0,0,0.0);
    
    chiS = singleNum(0);
    
    #if DEBUG_STATS
    cout << "x2val: " << chiS << endl;
    #endif
    
    try{
        return Statistics::chi2prob(chiS, betas.size());
    }catch(StatsException){
        chiS = -1;
        return 2.0;
    }
    
}

vector<vector<double> > LogisticRegression::invertMatrix(const vector<vector<double> > &mat){
    
    alglib::real_2d_array vInvInf;
    vInvInf.setlength(mat.size(), mat.size());
    
    for(unsigned int i=0;i < mat.size();i++){
        for(unsigned int j=0;j < mat.at(i).size() ; j++){
            vInvInf(i,j) = mat.at(i).at(j);
        }
    }
    
    alglib::matinvreport report;
    alglib::ae_int_t reportInfo;
    rmatrixinverse(vInvInf, reportInfo, report);
    if(reportInfo != 1){
        throw SingularMatrixEx();
    }
    // Check condition number.
    if (report.r1 < condition_number_limit){
        throw ConditionNumberEx(1.0/report.r1);
    }
    vector<vector<double> > returnMatrix = vector<vector<double> >(mat.size(), vector<double>(mat.size(), 0));
    
    for(unsigned int i=0;i<mat.size();i++){
        for(unsigned int j=0;j<mat.at(i).size();j++){
            returnMatrix[i][j] = vInvInf(i,j);
        }
    }
    
    return returnMatrix;
}

/**
 * Compute the likelhood ratio between null model (beta1) and non-null model (beta2)
 * 
 * D = -2 * sum ( y_i ln (pi_i / y_i) + (1 - y_i) ln ((1-pi_i) / (1 - y_i) ) ) for given model pi.
 * 
 * Compute G = D model 1 - D model 2
 * 
 * @return The likelihood statistic.
 * 
 */
double LogisticRegression::likelihoodRatio(const vector<double> &betas1, const vector<vector<double> > &data1, 
        const vector<double> &betas2, const vector<vector< double> > &data2, const vector<double> phen)
{
    
    // Compute each likelhood.
    double d1 = 0, d2 = 0;
    
    vector<double> expectedOne = expectedScores(betas1, data1);
    vector<double> expectedTwo = expectedScores(betas2, data2);
    
    for(unsigned int i=0; i < expectedOne.size(); i++){
        
        if(equal(phen.at(i) ,0)){
            d1 += log(1-expectedOne.at(i));
            d2 += log(1-expectedTwo.at(i));
        }else{
            d1 += log(expectedOne.at(i));
            d2 += log(expectedTwo.at(i));
        }
    }
            
    return -2 * ( d1 - d2 );
}

/*
 * Compute the expected scores for a set of people.
 * 
 * f(z) = exp(z) / (exp(z) + 1) where z = beta_1 x_1 + ... + beta_n x_n
 * 
 * @param betas The vector of betas.
 * @param data The vector of data.
 * @param vector<double> the expected scores.
 */
vector<double> LogisticRegression::expectedScores(const vector<double> &betas, const vector<vector<double> > &data){
    
    vector<double> ret;
    ret.reserve(data.at(0).size());
    for(unsigned int i=0; i < data.at(0).size(); ++i){
        
        double d = 0;
        for(unsigned int j=0; j < betas.size(); ++j){
            d += betas.at(j) * data.at(j).at(i);
        }
        double ed = exp(d);
        ret.push_back(ed/ (ed + 1));
    }
    return ret;
}

/*
 * Compute the inverse fisher information matrix, which is an approximation to the
 * covariance of the logistic model.
 * 
 * @param vector<vector<double>> data
 * @param vector<double> betas 
 * @return vector<vector<double>> the computed inverse matrix.
 * @return bool True if inverse worked.
 */
bool LogisticRegression::invFisherInformation(const vector<vector<double> > &data, const vector<double> &betas, vector<vector<double> > &returnMatrix){
    
    vector<double> expVals = expectedScores(betas, data);
    
    int sz = betas.size()-1;
    
    alglib::real_2d_array vInvInf;
    vInvInf.setlength(sz+1,sz+1);
    
    for(unsigned int indiv = 0; indiv < data.at(0).size(); indiv++){
        for(int i=0;i <= sz; i++){
            for(int j=0; j <= sz; j++){
                vInvInf(i,j) += expVals.at(indiv)*(1.0 - expVals.at(indiv))*data.at(i).at(indiv)*data.at(j).at(indiv);
            }
        }
    }
    
    alglib::matinvreport report;
    alglib::ae_int_t reportInfo;
    rmatrixinverse(vInvInf, reportInfo, report);
    if(reportInfo != 1) return false;
    // Check condition number.
    if (report.r1 < condition_number_limit){
        throw ConditionNumberEx(1.0/report.r1);
    }
    returnMatrix.clear();
    vector<double> temp;
    for(int i=0;i <= sz; i++){
        temp.push_back(0);
    }
    for(int i=0;i <= sz; i++){
        returnMatrix.push_back(temp);
    }
    
    for(int i=0;i <= sz; i++){
        for(int j=0; j <= sz; j++){
            returnMatrix.at(i).at(j) = vInvInf(i,j);
        }
    }
    
    return true;
}

/*
 * @depricated
 * 
 * NR implementation by RTG.
 *
 * Note on matrix vector multiplication:
 *  void matrixmatrixmultiply(const alglib::real_2d_array& a,
 *    int ai1,
 *    int ai2,
 *    int aj1,
 *    int aj2,
 *    bool transa,
 *    const alglib::real_2d_array& b,
 *    int bi1,
 *    int bi2,
 *    int bj1,
 *    int bj2,
 *    bool transb,
 *    double alpha,
 *    alglib::real_2d_array& c,
 *    int ci1,
 *    int ci2,
 *    int cj1,
 *    int cj2,
 *    double beta,
 *    alglib::real_1d_array& work);
 * 
 * transa is true if transposed.
 * Operation is: c = A * alpha * b + beta * C
 * 
 *  
 * @input const vector<vector<double>>  data holds explantory variables.  each inner vector is a variable.
 * @input const vector<double>          response holds response variable.
 * @input       vector<double>          Place to return the information matrix inverse.  Column major order.
 *                                      Expects correct size.
 */
vector<double> LogisticRegression::newtonRaphsonFast(const vector<vector<double> > &data, const vector<double> &response
                                                    , vector<vector<double> > &invInfMatrix, double startVal)
{
    
    
    vector<double> betas;
    
    
    // Variables used in the computation:
    alglib::real_1d_array tempBetas;
    alglib::real_2d_array tempData;
    alglib::real_1d_array oldExpY;
    alglib::real_2d_array tempDeriv; // holds data' * w * data.  Returned in last input param.
    alglib::real_1d_array expY;
    
    alglib::real_1d_array adjy;
    
    
    double stop_var = 1e-10;

    int iter = 0;
    int maxIter = 200;
    
    int numVars = data.size();
    if(numVars < 1){
        throw NewtonRaphsonFailureEx();
    }
    
    int numSamples = data.at(0).size();
    if(numSamples < 1){
        throw NewtonRaphsonFailureEx();
    }
    
    tempBetas.setlength(numVars);
    tempData.setlength(numSamples,numVars);
    oldExpY.setlength(numSamples);
    expY.setlength(numSamples);
    tempDeriv.setlength(numVars,numVars);
    adjy.setlength(numSamples);
    
    for(int i=0;i < numVars; i++){
        for(int j=0;j < numSamples; j++){
            tempData(j,i) = data.at(i).at(j);
        }
        tempBetas(i) = startVal;
    }
    for(int i=0;i < numSamples;i++){
        oldExpY(i) = -1;
        adjy(i) = 0; // makes valgrind happier.
    }

    while(iter < maxIter){
    
        //data * tempBetas
        matrixvectormultiply(tempData, 0, numSamples-1,0,numVars-1,false, 
                            tempBetas, 0, numVars-1, 1.0, 
                            adjy, 0, numSamples-1, 0.0);
        
        for(int i=0;i < numSamples; i++){
            expY(i) = 1 / (1 + exp(-adjy(i)));
        }
        
        // build deriv.
        double deriv = -100000;
        for(int i=0;i < numSamples; i++){
            if(expY(i) * (1 - expY(i)) > deriv)
                deriv = expY(i) * (1 - expY(i));
        }
        if(stop_var * 0.001 > deriv) deriv = stop_var * 0.001;
        
        // adjy = adjy + (y-expy) ./ deriv
        for(int i=0;i<numSamples;i++){
            adjy(i) = adjy(i) + (response.at(i) - expY(i)) / deriv;
        }
        
        // build data' * w * data
        alglib::real_1d_array work;
        work.setlength(numSamples); // This temporary workspace must be one larger than the longest
                                    // row or column that will be seen.  Otherwise, the program
                                    // crashes or - worse! - corrupts data.
        matrixmatrixmultiply(tempData,0,numSamples-1,0,numVars-1,true,
                            tempData,0,numSamples-1,0,numVars-1,false,deriv,
                            tempDeriv,0,numVars-1,0,numVars-1,0.0,work);

        #if DEBUG_NR
            cout << "A' * w * A " << endl;
            cout << tempDeriv(0,0) << " " << tempDeriv(0,1) << endl;
            cout << tempDeriv(1,0) << " " << tempDeriv(1,1) << endl;
            cout << endl;
        #endif
        
        alglib::matinvreport report;
        alglib::ae_int_t reportInfo;
        rmatrixinverse(tempDeriv, reportInfo, report);
        
        if( reportInfo != 1 ){
            throw SingularMatrixEx();
        }

        #if DEBUG_NR
            cout << "inv(A' * w * A) " << endl;
            cout << tempDeriv(0,0) << " " << tempDeriv(0,1) << endl;
            cout << tempDeriv(1,0) << " " << tempDeriv(1,1) << endl;
            cout << endl;
        #endif

        matrixvectormultiply(tempData,0,numSamples-1,0,numVars-1,true,
                            adjy,0,numSamples-1,deriv,
                            work,0,numVars-1,0.0);
        matrixvectormultiply(tempDeriv,0,numVars-1,0,numVars-1,false,
                            work,0,numVars-1,1.0,
                            tempBetas,0,numVars-1,0.0);
        
        #if DEBUG_NR
            cout << "Betas ";
            for(int i=0;i < numVars;i++) cout << tempBetas(i) << "  " ;
            cout << endl;
        #endif
        
        double stop = 0.0;
        for(int i=0;i < numSamples;i++){
            stop += abs(expY(i) - oldExpY(i));
        }
    
        if (stop < numSamples*stop_var){
            break;
        }
        
        oldExpY = expY;
        
        iter++;
    }
    
    if(iter == maxIter){
        throw NewtonRaphsonIterationEx();
    }
    
    betas.clear();
    for(int i=0;i<numVars;i++){
        betas.push_back(tempBetas(i));
    }

    for(int i=0;i<numVars;i++){
        for(int j=0;j<numVars;j++){
            invInfMatrix.at(i).at(j) = tempDeriv(i,j);
        }
    }
    
    return betas;
}
// Perform exact (and slower) NR test using the exact fisher information matrix.
vector<double> LogisticRegression::newtonRaphson(const vector<vector<double> > &data, const vector<double> &response, vector<vector<double> > &invInfMatrix, double startVal)
{

    vector<double> betas;
    
    // Variables used in the computation:
    alglib::real_1d_array tempBetas;
    alglib::real_2d_array tempData;
    alglib::real_2d_array tempDataTrans;
    alglib::real_1d_array oldExpY;
    alglib::real_2d_array hessian; // holds data' * w * data.  Returned in last input param.
    alglib::real_1d_array expY;
    alglib::real_1d_array W; // holds diagonal of the w matrix above.
    alglib::real_1d_array adjy;
    
    alglib::real_1d_array work;
    
    double stop_var = 1e-10;

    int iter = 0;
    int maxIter = 200;
    
    int numVars = data.size();
    if(numVars < 1){
        throw NewtonRaphsonFailureEx();
    }
    
    int numSamples = data.at(0).size();
    if(numSamples < 1){
        throw NewtonRaphsonFailureEx();
    }
    
    tempBetas.setlength(numVars);
    tempData.setlength(numSamples,numVars);
    tempDataTrans.setlength(numVars, numSamples);
    oldExpY.setlength(numSamples);
    expY.setlength(numSamples);
    hessian.setlength(numVars, numVars);
    adjy.setlength(numSamples);
    W.setlength(numSamples);
    
    work.setlength(numVars);
    
    for(int i=0;i < numVars; i++){
        for(int j=0;j < numSamples; j++){
            tempData(j,i) = data.at(i).at(j);
        }
        tempBetas(i) = startVal;
    }
    for(int i=0;i < numSamples;i++){
        oldExpY(i) = -1;
        adjy(i) = 0; // makes valgrind happier.
    }

    // End initial setup.
    // In each iteration, create a hessian and a first derivative.
    while(iter < maxIter){
        
        //adjy <- data * tempBetas (get new y guess)
        matrixvectormultiply(tempData, 0, numSamples-1,0,numVars-1,false, 
                            tempBetas, 0, numVars-1, 1.0, 
                            adjy, 0, numSamples-1, 0.0);
        
        // adjy = 1 / (1 + exp(-adjy))  
        for(int i=0;i < numSamples; i++){
            expY(i) = 1 / (1 + exp(-adjy(i)));
        }
        
        // build deriv.
        for(int i=0;i < numSamples; i++){
            W(i) = expY(i) * (1 - expY(i));
        }
        
        // adjy = adjy + (y-expy) ./ deriv
        for(int i=0;i<numSamples;i++){
            adjy(i) = adjy(i) + (response.at(i) - expY(i)) / W(i);
        }
        
        // build data' * w * data
        // set to hessian.  
        // Also doing secondary computation (see inside)
        for(int i=0; i < numVars; i++){
            for(int j=0; j < numVars; j++)
                hessian(i,j) = 0.0;
        }
        
        for(int indiv=0; indiv < numSamples; ++indiv){
            for(int i = 0; i < numVars; i++){
                for(int j = 0; j < numVars; j++){
                    hessian(i,j) += W(indiv) * tempData(indiv, j) * tempData(indiv, i);
                }
                
                // NOTE: as a speedup, I'm also computing X' * W
                tempDataTrans(i, indiv) = W(indiv) * tempData(indiv, i);
                
            }
        }
        
        alglib::matinvreport report;
        alglib::ae_int_t reportInfo;
        rmatrixinverse(hessian, reportInfo, report);
        if(reportInfo != 1 ){
            throw SingularMatrixEx();
        }
        
        // Check condition number.
        if (report.r1 < condition_number_limit){
            throw ConditionNumberEx(1.0/report.r1);
        }
        
        // work <- X'W * adjy
        matrixvectormultiply(tempDataTrans,0,numVars-1,0,numSamples-1,false,
                            adjy,0,numSamples-1,1,
                            work,0,numVars-1,0.0);
        // tempBetas <= invHessian * work
        matrixvectormultiply(hessian,0,numVars-1,0,numVars-1,false,
                            work,0,numVars-1,1.0,
                            tempBetas,0,numVars-1,0.0);
        
        
        #if DEBUG_NR
            cout << "Betas ";
            for(int i=0;i < numVars;i++) cout << tempBetas(i) << "  " ;
            cout << endl;
        #endif
        
        double stop = 0.0;
        // Could be computed as a 1-norm.
        // This should be done as a sum of abs diff.
        for(int i=0;i < numSamples;i++){
            stop += abs(expY(i) - oldExpY(i));
        }

        if (stop < numSamples*stop_var){
            break;
        }
        
        oldExpY = expY;
        
        iter++;
    }
    
    if(iter == maxIter){
        
        throw NewtonRaphsonIterationEx();
    }
    
    betas.clear();
    for(int i=0;i<numVars;i++){
        betas.push_back(tempBetas(i));
    }

    for(int i=0;i<numVars;i++){
        for(int j=0;j<numVars;j++){
            invInfMatrix.at(i).at(j) = hessian(i,j);
        }
    }
    
    //dumpMatrix(invInfMatrix);
    
    return betas;
    
}

void LogisticRegression::dumpMatrix(const vector<vector<double> > &data){
    
    for (unsigned int i=0; i < data.size(); ++i){
        for (unsigned int j =0; j < data.at(i).size(); ++j){
            cout << data[i][j] << " ";
        }
        cout << endl;
    }
    
}

// Utilities
double LogisticRegression::abs(double x){
    if(x > 0) return x;
    return -1 * x;
}

