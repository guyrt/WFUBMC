#include <gtest/gtest.h>

#include "../engine/utils/lr.hh"
#include "../engine/utils/vecops.hh"

class LR_Engine_Test : public ::testing::Test {

	protected:

	vector<vector<double> > inMat;
	vector<double> phenotype;
	vector<double> ones;
	LogisticRegression lr;	

	virtual void SetUp(){
	ones = vector<double>(6,1.0);
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
	
	phenotype.push_back(0);
	phenotype.push_back(1);
	phenotype.push_back(0);
	phenotype.push_back(1);
	phenotype.push_back(0);
	phenotype.push_back(1);
	
	inMat.push_back(ones);
	inMat.push_back(m1); 
	inMat.push_back(m2);
	
	}

};



TEST_F(LR_Engine_Test, test_newton_raphson) {
  
	vector<vector<double> > invInfMatrix = vecops::getDblVec(inMat.size(), inMat.size());
	vector<double> betas = lr.newtonRaphson(inMat, phenotype, invInfMatrix);
	ASSERT_DOUBLE_EQ(betas[0], 0.36750847348763349);
}


