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

class LR_Separable_Test : public ::testing::Test{

	protected: 

	vector<vector<double> > cov;
        vector<double> response;

	void SetUp(){

        for (int i=0;i < 4; i++)
                response.push_back(1.0);
        for (int i=0;i < 4; i++)
                response.push_back(0.0);
                
        for (int i=0; i < 2; i++)
                cov.push_back( vector<double>() );
        for (int i=0; i < 4; i++){
                cov[0].push_back(1.0);
        }
        for (int i=0; i < 4; i++){
                cov[0].push_back(2.0);
        }
        
        for (int i=0; i < 8; i++){
                cov[1].push_back(i);
        }
	}

}; 

//--------------------------------------------------------------------------        

TEST_F(LR_Engine_Test, test_newton_raphson) {
  
	vector<vector<double> > invInfMatrix = vecops::getDblVec(inMat.size(), inMat.size());
	vector<double> betas = lr.newtonRaphson(inMat, phenotype, invInfMatrix);
	ASSERT_DOUBLE_EQ(betas[0], 0.36750847348763349);
	ASSERT_DOUBLE_EQ(betas[1], 1.4467808518431764);
	ASSERT_DOUBLE_EQ(betas[2], -1.5618476151485159);
}

TEST_F(LR_Separable_Test, test_separable) {
	
	LogisticRegression lr;
        int separableColumn = lr.dataIsSeparable(cov, response);
        ASSERT_EQ( separableColumn, 0);  

}
