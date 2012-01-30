//      genostats.hh
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

/*
 * Compute genotype statistics for single SNP with or without covariates.
 *
 * Computes:
 * 	Lack of fit.
 * 	Dom, Add, Rec
 *  2 Deg Freedom test.
 */

#include "../../param/engine_param_reader.h"

#include "../engine.h"
#include "../data_plugin.h" // only for error printing.
#include "../utils/lr.hh"
#include "../utils/statistics.h"
#include "../utils/vecops.hh"
#include "../output/snpgwa_out.hh" // this gives us access to the writeout format

using namespace std;

#define DEBUG_ADD_STAT 0
#define DEBUG_GENO_TDF 0
#define DEBUG_GENO_SINGLE 0

class GenoStats {

	public :
		GenoStats(DataAccess *d, EngineParamReader *p);
		void prepGenoStatsForOutput(int snp, GenoStatsResults &g);

	private :
	
		struct errorInformation {
			int snp;
			int numCovariates; // used to test for separation.
			string message;
		};
	
		DataAccess *data; // pointer to data.
		const EngineParamReader *params;
		
		void calculateSensSpec(int snp, GenoStatsResults &result, double *geno_bins);
		LRStats runSingleLRTest(const vector<vector<double> > &in, const vector<double> &phen, vector<double> &betas, errorInformation failMessage);
		double runLRTest(const vector<vector<double> > &in, const vector<double> &phen, double &chiS, errorInformation failMessage);

		/**
		 * Handle an exception. This is a convenience method.
		 * 
		 */
		void handleException(LRStats &l, double &startVal, int &retry, const string &message, errorInformation errorData);

		int currentSNP;

};
