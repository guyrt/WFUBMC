#include "genostats.hh"
#include "../linalg/linalg.h" // used for condition number.
#include "../linalg/alglibinternal.h"
#include "../linalg/blas.h"
/**
 * Assign the data access.
 *
 * Also, prep covariate matrix by zero meaning all covariates and preloading.
 */
GenoStats::GenoStats(DataAccess *d, EngineParamReader *p){
    data = d;
    params = p;
}

/*
 * This is the primary entry point for this class.
 *
 * Compute all genotypic statistics.
 *
 * I've made some use of outside methods, but this one is still pretty long.
 * Many methods are so entertwined for efficiency that it would take 5 or 6 parameters to
 * break out any other pieces.  Still, feel free to try.
 */
void GenoStats::prepGenoStatsForOutput(int snp, GenoStatsResults &results){

    bool hasCov = false;
    currentSNP = snp;

    vector<double> phen_vec;
    vector<double> dom, add, rec, lof, twodegfree1, twodegfree2;

    vector<vector<double> > add_in, dom_in, rec_in, twodegfree_in, lof_in;
    vector<vector<double> > cov;
    vector<double> ones;
    vector<double> cov_cum_sum;

    double geno_bins[6]; // Used to bin the data by genotype.
    geno_bins[0] = geno_bins[1] = geno_bins[2] = geno_bins[3] = geno_bins[4] = geno_bins[5] = 0;

    /* Prep covariate matrix */
    vector<double> *tmp = data->get_covariates(0);
    if(tmp != NULL){
        hasCov = true;
        for(unsigned int i=0;i<tmp->size();i++){
            vector<double> a;
            cov.push_back(a);
            cov_cum_sum.push_back(0);
        }
    }

    double temp;
    /* Prep genotypes and fill cov matrix
     *
     * This pattern is used often enough to potentially extract.
     * Pattern:
     *      Loop through nonmissing data
     *      Apply some function depending on the case data.
     *      Those functions include:
     *          pushing something onto a vector
     *          increasing an array value.
     *      Push covariate matrix on.
     *      Optionally, push nonzero phenotypes onto vector.
     *
     * Used by: right here, cont_genostats.cpp, haplotypes testing
     * ld.cpp
     *
     * Must: include way to consider multiple SNPs.
     */
    for(int i=0; i<data->pheno_size(); i++ ){
        // push onto stacks depending on the case/cntrl status.
        if(data->get_data(i)->at(snp) != 0){
            temp = data->get_phenotype(i)-1;
            phen_vec.push_back(temp);
            ones.push_back(1.0);
            switch(data->get_data(i)->at(snp)){
                case 1:
                    add.push_back(-1);
                    dom.push_back(0);
                    rec.push_back(0);

                    twodegfree1.push_back(0);
                    twodegfree2.push_back(0);

                    lof.push_back(1.0);

                    geno_bins[static_cast<int>(temp) * 3]++;

                break;
                case 2:
                    add.push_back(0);
                    dom.push_back(1);
                    rec.push_back(0);

                    twodegfree1.push_back(0);
                    twodegfree2.push_back(1);

                    lof.push_back(-2.0);

                    geno_bins[static_cast<int>(temp) * 3+1]++;
                break;
                case 3:
                    add.push_back(0);
                    dom.push_back(1);
                    rec.push_back(0);

                    twodegfree1.push_back(0);
                    twodegfree2.push_back(1);

                    lof.push_back(-2.0);

                    geno_bins[static_cast<int>(temp) * 3+1]++;
                break;
                case 4:
                    add.push_back(1);
                    dom.push_back(1);
                    rec.push_back(1);

                    twodegfree1.push_back(1);
                    twodegfree2.push_back(0);

                    lof.push_back(1.0);

                    geno_bins[static_cast<int>(temp) * 3 + 2]++;
                break;
                default:
                // won't happen.
                break;
            }
            if(hasCov){
                vector<double> *t = data->get_covariates(i);
                for(unsigned int j=0;j<t->size();j++){
                    cov.at(j).push_back(t->at(j));
                    cov_cum_sum.at(j)+=t->at(j);
                }
            }
        }
    }

    /* Make all matrices
     */
    if(hasCov){
        add_in = cov;
        dom_in = cov;
        rec_in = cov;
        twodegfree_in = cov;
        lof_in = cov;
    }

    add_in.push_back(ones);
    add_in.push_back(add);

    dom_in.push_back(ones);
    dom_in.push_back(dom);

    rec_in.push_back(ones);
    rec_in.push_back(rec);

    twodegfree_in.push_back(twodegfree2);
    twodegfree_in.push_back(twodegfree1);
    twodegfree_in.push_back(ones);

    lof_in.push_back(ones);
    lof_in.push_back(lof);

    vector<double> dom_beta, add_beta, rec_beta, lof_beta;

    errorInformation errorData = {snp, cov.size(), "Additive test "};

    LRStats add_l = runSingleLRTest(add_in, phen_vec, add_beta, errorData);
    results.addTestStat = add_l.testStat;
    results.addOR = add_l.OR;
    results.addUCI = add_l.UCI;
    results.addLCI = add_l.LCI;
    results.addPVal = add_l.pVal;
    #if DEBUG_GENO_SINGLE
    if(add_beta.size() > 0){
    cout << "Add size: " << add_in.size() << " " << add_in.at(0).size() << endl;
    cout << "Additive beta: " ;
    for(unsigned int i=0;i < add_beta.size(); i++){
        cout << add_beta.at(i) << " ";
    }cout << endl;
    }
    #endif

    errorData.message = "Dominant test ";
    LRStats dom_l;

    dom_l = runSingleLRTest(dom_in, phen_vec, dom_beta, errorData);
    results.domTestStat = dom_l.testStat;
    results.domOR = dom_l.OR;
    results.domUCI = dom_l.UCI;
    results.domLCI = dom_l.LCI;
    results.domPVal = dom_l.pVal;


    errorData.message = "Recessive test ";
    LRStats rec_l;
    rec_l = runSingleLRTest(rec_in, phen_vec, rec_beta, errorData);
    results.recTestStat = rec_l.testStat;
    results.recOR = rec_l.OR;
    results.recUCI = rec_l.UCI;
    results.recLCI = rec_l.LCI;
    results.recPVal = rec_l.pVal;


    errorData.message = "Lack of fit test ";
    LRStats lof_l;

    lof_l = runSingleLRTest(lof_in, phen_vec, lof_beta, errorData);
    results.lofTestStat = lof_l.testStat;
    results.lofPVal = lof_l.pVal;

    errorData.message = "Two deg freedom test ";

    results.twodegPVal = runLRTest(twodegfree_in, phen_vec, results.twodegTestStat, errorData);

    calculateSensSpec(snp, results, geno_bins);
}

/**
 * Calculate sensitivity and specificity and CStatistic (all unadjusted)
 */
void GenoStats::calculateSensSpec(int snp, GenoStatsResults &results, double *geno_bins){

    if(geno_bins[0] + geno_bins[1] > 0 && geno_bins[0] + geno_bins[2] > 0
        && geno_bins[3] + geno_bins[4] > 0 && geno_bins[3] + geno_bins[5] > 0
        && geno_bins[4] + geno_bins[5] > 0 && geno_bins[1] + geno_bins[2] > 0)
    {

        // Compute sens/spec/CStat.
        // WARNING: This is an unadjusted statistic
        results.domSens = (geno_bins[1] + geno_bins[2]) / (geno_bins[0] + geno_bins[1] + geno_bins[2]);
        results.domSpec = (geno_bins[3]) / (geno_bins[3] + geno_bins[4] + geno_bins[5]);
        results.domCStat = ( results.domSens + results.domSpec ) / 2;

        results.recSens = (geno_bins[2]) / (geno_bins[0] + geno_bins[1] + geno_bins[2]);
        results.recSpec = (geno_bins[3] + geno_bins[4]) / (geno_bins[3] + geno_bins[4] + geno_bins[5]);
        results.recCStat = ( results.recSens + results.recSpec ) / 2;

        results.addSensNNRN = (geno_bins[1]) / (geno_bins[0] + geno_bins[1]);
        results.addSpecNNRN = (geno_bins[3]) / (geno_bins[3] + geno_bins[4]);

        results.addSensNNRR = (geno_bins[2]) / (geno_bins[0] + geno_bins[2]);
        results.addSpecNNRR = (geno_bins[3]) / (geno_bins[3] + geno_bins[5]);

        results.addSensNRRR = (geno_bins[2]) / (geno_bins[1] + geno_bins[2]);
        results.addSpecNRRR = (geno_bins[4]) / (geno_bins[4] + geno_bins[5]);

        double p11, p12, p13, p21, p22, p23;
        p11 = (geno_bins[0]) / (geno_bins[0] + geno_bins[1] + geno_bins[2]);
        p12 = (geno_bins[1]) / (geno_bins[0] + geno_bins[1] + geno_bins[2]);
        p13 = (geno_bins[2]) / (geno_bins[0] + geno_bins[1] + geno_bins[2]);
        p21 = (geno_bins[3]) / (geno_bins[3] + geno_bins[4] + geno_bins[5]);
        p22 = (geno_bins[4]) / (geno_bins[3] + geno_bins[4] + geno_bins[5]);
        p23 = (geno_bins[5]) / (geno_bins[3] + geno_bins[4] + geno_bins[5]);

        results.addCStat = 0.5*p23*p13 + 0.5*p22*(2*p13 + p12) + 0.5*p21*(2*p13 + 2*p12 + p11);

    }else{
        results.addSensNNRN = results.addSpecNNRN = results.addSensNNRR = results.addSpecNNRR = -1;
        results.addSensNRRR = results.addSpecNRRR = -1;
        results.addCStat = -1;
        results.domSens = results.domSpec = results.domCStat = -1;
        results.recSens = results.recSpec = results.recCStat = -1;

        // Write log messages.
        string tempS;
        int tempP;
        string name1;
        data->get_map_info(snp, tempS, name1, tempP);
        stringstream ss;
        ss << "Snpgwa for SNP " << name1 << " test for sensitivity and specificity: zero denominator detected. Tests not run.";
        Logger::Instance()->writeLine(ss.str());
    }
}

// Compute a single LR Wald test and return results.
LRStats GenoStats::runSingleLRTest(const vector<vector<double> > &in, const vector<double> &phen, vector<double> &betas, errorInformation errorData){

    vector<vector<double> > inv_infmatrix;

    LRStats l;
    l.fillDefault();
    LogisticRegression lr(params->getRegressionConditionNumberThreshold());

    inv_infmatrix = vecops::getDblVec(in.size(), in.size());

    int retry = 0;
    double startVal = 0;  // value to start betas with.

    while(retry < 3){
        try{
            betas = lr.newtonRaphson(in, phen, inv_infmatrix, startVal);
            l = lr.getSingleStats(betas, inv_infmatrix, betas.size()-1);

            if(l.OR != l.OR){
                stringstream ss;
                ss << "Odds ratio was NaN";
                retry = 4;
                handleException(l, startVal, retry, ss.str(), errorData);
                l.OR = -1;
                l.UCI = -1;
                l.LCI = -1;
            }else if(l.OR >= 10000.0){
                stringstream ss;
                ss << "Large odds ratio detected and rewritten: OR was " << l.OR << ".";
                retry = 4;
                handleException(l, startVal, retry, ss.str(), errorData);

                l.OR = 9999.0;
                l.UCI = 9999.0;
                l.LCI = 9999.0;
            }


            break;
        }catch(NewtonRaphsonFailureEx){
            handleException(l, startVal, retry, "newton-raphson setup failure.", errorData);
        }catch(NewtonRaphsonIterationEx){
            handleException(l, startVal, retry, "max iterations hit.", errorData);
        }catch(SingularMatrixEx){
            handleException(l, startVal, retry, "information matrix was singular.", errorData);
        }catch(alglib::ap_error err){
            stringstream ss;
            ss << "linalg exception: " << err.msg;
            handleException(l, startVal, retry, ss.str(), errorData);
        }catch(ConditionNumberEx err){
            //Handle special.
            retry = 10;
            // The condition number of the information matrix is large. Check for separation.
            int separableVariable = lr.dataIsSeparable(in, phen);

            string tempS;
            int tempP;
            string name1;
            data->get_map_info(errorData.snp, tempS, name1, tempP);

            if (separableVariable < 0){
                // Error: poor conditioning.
                stringstream ss;
                ss << "Poor conditioning in information matrix. " ;
                ss << "Condition number (1-norm) is " << err.conditionNumber;
                handleException(l, startVal, retry, ss.str(), errorData);

            }else{
                // Error: the poor conditioning is due to separable variable.
                string message;
                if (separableVariable == errorData.numCovariates - 1){
                    message = "SNP separates variable.";
                }else if (separableVariable < errorData.numCovariates - 2){
                    message = "Separable by covariate.";
                }
                handleException(l, startVal, retry, message, errorData);
            }
            return l;
        }

    }
    return l;
}

/*
 * Compute Wald statistic for two coefficients.
 * Return p-value.
 *
 * TODO: make this work for > 2 coefficients.
 */
double GenoStats::runLRTest(const vector<vector<double> > &in, const vector<double> &phen, double &chiS, errorInformation errorData){

    vector<vector<double> > inv_infmatrix;
    vector<double> betas;

    int retry = 0;
    double startVal = 0.0;

    inv_infmatrix = vecops::getDblVec(in.size(), in.size());

    while(retry < 3){
        try{
            LogisticRegression lr(params->getRegressionConditionNumberThreshold());
            betas = lr.newtonRaphson(in, phen, inv_infmatrix, startVal);

            vector<double> tdfb;
            tdfb.push_back(betas.at(betas.size()-3));
            tdfb.push_back(betas.at(betas.size()-2));
            vector<vector<double> > td;
            td.push_back(tdfb);
            td.push_back(tdfb);

            td.at(0).at(0) = inv_infmatrix.at(betas.size()-3).at(betas.size()-3);
            td.at(1).at(0) = inv_infmatrix.at(betas.size()-3).at(betas.size()-2);
            td.at(0).at(1) = inv_infmatrix.at(betas.size()-2).at(betas.size()-3);
            td.at(1).at(1) = inv_infmatrix.at(betas.size()-2).at(betas.size()-2);

            #if DEBUG_GENO_TDF
            cout << "Results: " << currentSNP << endl;
            cout << td.at(0).at(0) << " " << td[1][0] << endl;
            cout << td.at(0).at(1) << " " << td[1][1] << endl;
            cout << endl;
            cout << "Using " << tdfb.at(0) << " " << tdfb.at(1) << endl;
            cout << "Covar:" << endl;
            for(unsigned int ij=0;ij < inv_infmatrix.size();ij++){
                for(unsigned int i=0;i < inv_infmatrix.at(0).size();i++)
                    cout << inv_infmatrix.at(ij).at(i) << " ";
                cout << endl;
            }
            #endif

            return lr.getStats(tdfb, td, chiS);
        }catch(NewtonRaphsonFailureEx){


            string tempS;
            int tempP;
            string name1;
            data->get_map_info(errorData.snp, tempS, name1, tempP);

            stringstream ss;
            ss << "SNPGWA for SNP " << name1 << " " << errorData.message << "Setup error in Logistic Regression" << endl;
            Logger::Instance()->writeLine(ss.str());

            return 2.0;
        }
        catch(NewtonRaphsonIterationEx){
            retry++;
            if(retry == 1) {startVal = 0.5;}
            else if(retry == 2){ startVal = -0.5;}
            else{

                string tempS;
                int tempP;
                string name1;
                data->get_map_info(errorData.snp, tempS, name1, tempP);
                stringstream ss;
                ss << "SNPGWA for SNP " << name1 << " " << errorData.message << "Failed to compute: Max iterations hit." << endl;
                Logger::Instance()->writeLine(ss.str());
                return 2.0;
            }
        }
        catch(SingularMatrixEx){

            retry++;
            if(retry == 1) {startVal = 0.5;}
            else if(retry == 2){ startVal = -0.5;}
            else{

                string tempS;
                int tempP;
                string name1;
                data->get_map_info(errorData.snp, tempS, name1, tempP);
                stringstream ss;
                ss << "SNPGWA for SNP " << name1 << " " << errorData.message << "Failed to compute: the matrix is singular.  Three singular matrix exceptions." << endl;
                Logger::Instance()->writeLine(ss.str());
                return 2.0;
            }
        }catch(ConditionNumberEx err){

            // The condition number of the information matrix is large. Check for separation.
            LogisticRegression lr(params->getRegressionConditionNumberThreshold());
            int separableVariable = lr.dataIsSeparable(in, phen);

            string tempS;
            int tempP;
            string name1;
            data->get_map_info(errorData.snp, tempS, name1, tempP);


            if (separableVariable < 0){
                // Error: poor conditioning.
                stringstream ss;
                ss << "Snpgwa for SNP " << name1 << " " << errorData.message << ": Poor conditioning in information matrix.";
                ss << "Condition number (1-norm) is " << err.conditionNumber << endl;
                Logger::Instance()->writeLine(ss.str());
            }else{
                // Error: the poor conditioning is due to separable variable.
                stringstream ss;
                ss << "Snpgwa for SNP " << name1 << " " << errorData.message << ": ";
                if (separableVariable == errorData.numCovariates - 1){
                    ss << "SNP separates variable.";
                }else if (separableVariable < errorData.numCovariates - 2){
                    ss << "Separable by covariate.";
                }
                ss << endl;
                Logger::Instance()->writeLine(ss.str());
            }
            return 2.0;
        }catch(alglib::ap_error err){
            retry++;
            if(retry == 1) {startVal = 0.5;}
            else if(retry == 2){ startVal = -0.5;}
            else{
                string tempS;
                int tempP;
                string name1;
                data->get_map_info(errorData.snp, tempS, name1, tempP);

                stringstream ss;

                ss << "SNPGWA for SNP " << name1 << " " << errorData.message << "Failed to compute: linear algebra exception: " << err.msg << endl;
                Logger::Instance()->writeLine(ss.str());
                return 2.0;
            }
        }
    }

    return 2.0;
}

void GenoStats::handleException(LRStats &l, double &startVal, int &retry, const string &message, errorInformation errorData){

    l.fillDefault();
    retry++;
    if(retry == 1) {startVal = 0.5;}
    else if(retry == 2){ startVal = -0.5;}
    else{
        string tempS;
        int tempP;
        string name1;
        data->get_map_info(errorData.snp, tempS, name1, tempP);
        stringstream ss;
        ss << "SNPGWA for SNP " << name1 << " " << errorData.message << "Failed to compute: " << message << endl;
        Logger::Instance()->writeLine(ss.str());

    }

}

