#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE LogisticRegression

#include <boost/test/unit_test.hpp>

#include "../../../engine/utils/lr.hh"
#include "../../../engine/utils/vecops.hh"

BOOST_AUTO_TEST_SUITE(LR_ENGINE)


BOOST_AUTO_TEST_CASE( test_get_single_stats )
{
	
	
}

BOOST_AUTO_TEST_CASE( test_newton_raphson)
{
		
	vector<vector<double> > inMat;
	
	vector<double> ones(6,1.0);
	
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
	
	inMat.push_back(ones);
	inMat.push_back(m1); 
	inMat.push_back(m2);
	
	LogisticRegression lr;
	vector<vector<double> > invInfMatrix = vecops::getDblVec(inMat.size(), inMat.size());
	vector<double> betas = lr.newtonRaphson(inMat, m3, invInfMatrix);
	
	BOOST_CHECK_CLOSE( betas[0], 0.367508472, 0.0001);
	BOOST_CHECK_CLOSE( betas[1], 1.44678, 0.0001);
	BOOST_CHECK_CLOSE( betas[2], -1.5618476, 0.0001);

}

BOOST_AUTO_TEST_CASE( test_separable )
{
	vector<vector<double> > cov;
	vector<double> response;
	
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
	
	LogisticRegression lr;
	int separableColumn = lr.dataIsSeparable(cov, response);
	
	BOOST_CHECK_EQUAL( separableColumn, 0);
	
}

BOOST_AUTO_TEST_SUITE_END()
