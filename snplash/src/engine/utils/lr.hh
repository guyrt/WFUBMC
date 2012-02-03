//      lr.hh
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

/**
 * @class LogisticRegression
 * 
 * @author Richard T. Guy
 * @date July, 2010
 * 
 * Logistic regression computation engine.  
 */

#include "../engine.h"
#include "../utils/exceptions.h"
#include "../utils/statistics.h"
#include "../utils/float_ops.hh"
#include "vecops.hh"

#include <math.h>



#define DEBUG_NR 0
#define DEBUG_STATS 0

#ifndef LOGISTICREG_H
#define LOGISTICREG_H

using namespace std;

class StatsFillable{
	public:
		virtual void fillDefault() = 0;
		virtual ~StatsFillable(){};
};

class LRStats : public StatsFillable {
	
	public:
	double OR;
	double UCI;
	double LCI;
	double testStat;
	double pVal;
	double invInf;
	
	inline void fillDefault(){
		OR = -1;
		UCI = -1;
		LCI = -1;
		pVal = 2.0;
		invInf = -1;
		testStat = -1;
	}
};

class LogisticRegression {
	
	public :
		LogisticRegression();
		LogisticRegression(double conditionNumber);
	
		LRStats getSingleStats(const vector<double> &betas, const vector<vector<double> > &invInf, int index);
		double getStats(const vector<double> &betas, const vector<vector<double > > invInf, double &chiS);
	
		double likelihoodRatio(const vector<double> &betas1, const vector<vector<double> > &data1, 
			const vector<double> &betas2, const vector<vector< double> > &data2, const vector<double> phen);
	
		vector<double> expectedScores(const vector<double> &betas, const vector<vector<double> > &data);
	
		/* Calculate the Newton-Raphson algorithm for matrix of explanatory variables
		 * and vector of reponse variables.  The last matrix is the infomation matrix inverse
		 * used to compute stats on the model. */
		vector<double> newtonRaphsonFast(const vector<vector<double> > &data, const vector<double> &response, vector<vector<double> > &invInfMatrix, double startVal = 0.0);
		vector<double> newtonRaphson(const vector<vector<double> > &data, const vector<double> &response, vector<vector<double> > &invInfMatrix, double startVal = 0.0);
		
		bool invFisherInformation(const vector<vector<double> > &data, const vector<double> &betas, vector<vector<double> > &returnMatrix);
		
		int dataIsSeparable(const vector<vector<double> > &data, const vector<double> &response);
		
		inline double getConditionNumberLimit(){return condition_number_limit;}
		inline void setConditionNumberLimit(double l){condition_number_limit = l;}

	private :
		
		vector<vector<double> > invertMatrix(const vector<vector<double> > &mat);
		double abs(double);
		void dumpMatrix(const vector<vector<double> > &data);
		double variance(const vector<double> &data);
		
		double condition_number_limit; // Stored as inverse!
		static const double SEPARABLE_THRESHOLD = 0.98;
	
};

#endif
