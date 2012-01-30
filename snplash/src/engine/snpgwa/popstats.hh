//      harwein.hh
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
 * Contains all the programs for computing population statistics:
 * 		* Missing
 * 		* Num case/control
 * 		* Missing statistics
 * 
 */

#include <deque>
#include "../engine.h"
#include "../utils/statistics.h"
#include "../utils/vecops.hh"
#include "../output/snpgwa_out.hh" // this gives us access to the writeout format

#ifndef POPSTATS_H
#define POPSTATS_H

#define DEBUG_HWEQUIL 0
#define DEBUG_POPSTATS_PULL 0

using namespace std;

class PopStats {
	
	public:
		PopStats(DataAccess *);
	
		/* Main handle that will print everything to a print object */
		bool prepPopStatsForOutput(int snp, PopStatsResults &results);
	
		/* Calculate MAF for a given SNP */
		bool minorAlleleFreq(int snp);
		/* Calculate number in each category for given SNP */
		bool numCategories(int snp);
		/* Calculate all missing statistics */
		bool missingStats(int snp);
		/* Calculate Hardy-Weinburg Equilibrium */
		bool hwEquilibruim(int snp);
		/* Compute exact HWE test for counts listed. */
		static double exactTest(double numPP, double numPQ, double numQQ);
	protected:
		
		DataAccess *data;
		void reinit(); // Initialize all variables.
	
		/* OPTIMIZING STEP */
		int lastBuild;
		vector<short> snp_data;
		vector<double> phen_data;
		bool pullData(int build); // builds data.
		
		
		/* All of the following values are candidates for computing */
		/* Following values are computed because it's cheap in the pullData step. */
		int numMissingCase, numMissingCntrl;
		int numMissingTotal; // Note: this includes individuals with missing phenotype.
		int casesNonMissing, cntrlsNonMissing;
		/* Following values are computed in numCategories */
		int numTotalCases, numTotalCntrls;
		/* Following values are computed in minorAlleleFreq */
		double minorAlleleFreqCases, minorAlleleFreqCntrls;
		int ppCases, pqCases, qqCases, ppCntrls, pqCntrls, qqCntrls;
		/* Following values are computed in missingStatistics */
		double percentMissingTotal, percentMissingCases, percentMissingControls;
		double differentialMissingPval, differentialMissingOR;
		/* Following values are computed in hwEquilibrium */
		double caseExpTotal[3]; 
		double caseChiSquare, caseChiSquarePValue, caseExactTestPVal;
		double cntrlExpTotal[3];
		double cntrlChiSquare, cntrlChiSquarePValue, cntrlExactTestPVal;
		double cmbdExpTotal[3];
		double cmbdChiSquare, cmbdChiSquarePValue, cmbdExactTestPVal;
};
#endif
