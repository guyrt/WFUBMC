//      zaykin.hh
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
 * @class Zaykin
 * 
 * Computation engine for Zaykin's method for haplo-genotype analysis.
 * 
 * Performs dom, add, rec tests per haplotype.
 * Performs global test.
 * 
 * Requires that a phenotype and covariates be entered.
 * Will cull haplotypes with global freq < 0.05 and will either reweight
 * the rest or keep the same depending on user flags.
 * 
 * @author Richard T. Guy
 */

#ifndef ZAYKIN_H
#define ZAYKIN_H

#include "statistics.h"
#include "../em/em.h"
#include "../../param/param_reader.h"
#include "../../param/engine_param_reader.h"
#include "../engine.h"
#include "vecops.hh"
#include "../../logger/log.hh"

#include "lr.hh"

using namespace std;

#define DEBUG_ZAY_PROGRESS 0



/**
 * Holds information that is returned from a likelihood ratio test.
 */
class ZaykinStatsInfo : public StatsFillable {

	public:
		double chiSqStat;
		int degFree;
		double OR;
		double LCI;
		double UCI;
		
	inline void fillDefault(){
		chiSqStat = -1;
		degFree = -1;
		OR = -1;
		LCI = -1;
		UCI = -1;
	};
};

class ZaykinGlobalStatsResults : public StatsFillable {
	
	public:
		int degFreedom;
		double testStat;
		double pvalue;
		
		inline void fillDefault(){
			degFreedom = -1;
			testStat = -1.0;
			pvalue = 2.0;
		}; 
	
};


class Zaykin {
	
	friend class HaploStats;
	
	public : 
	
		Zaykin(EngineParamReader *p);
		~Zaykin();
		
		void setup(const vector<EMPersonalProbsResults>  &, int numHaps, bool reweight);
		
		ZaykinGlobalStatsResults runGlobal(); // Return pvalue.
		vector<double> runAllIndiv(vector<double> &testStats, vector<int> &degsFreedom); // Run all tests.
		vector<ZaykinStatsInfo> runAllHaplotypes();
		
		static int test();
		
		void setPhenotype(const vector<double> &p, const vector<vector<double> > &c);
		/**
		 * Sets the keepThresh limit, which is min number of chromosomes
		 * that must contain a haplotype for it to be used in testing.
		 * @param k New threshold.
		 */
		void setKeepThresh(int k){keepThresh = k;}
		
		void setErrorInformation(int snp, int numCovariates, DataAccess *data);
		
	protected : 
	
		const EngineParamReader *params;
	
		vector<double> phenotype;
		vector<vector<double> > cov;	

		vector<vector<double> > matrix; // each inner vector is a single individual's haplogenotype.
										// outer vector loops through individuals/
		vector<bool> usedHaplotype; // Holds true if the haplotype was used.

		int recode12_01(int i); // recode a number using 1,2 to 0,1.
		
		void dump(vector<vector<double> > m);
		
		void prepHaplotypes(vector<double> &ones, vector<vector<double> > &haps);
		ZaykinStatsInfo runLikelihoodRatio(const vector<double> &haps, const vector<double> &ones);

		int keepThresh;

		// Error data:
		int snp;
		int numCovariates;
		DataAccess *data;
		void handleException(StatsFillable &stats, double &startVal, int &retry, const string &message);
};

#endif

