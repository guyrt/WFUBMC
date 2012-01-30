//      cont_genostats.hh
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

#ifndef CONT_GENOSTATS_H
#define CONT_GENOSTATS_H

using namespace std;

// turn on to get more information about stats
#define CONT_GENOSTATS_STATS 0

#include "../engine.h"
#include "../utils/statistics.h"
#include "../output/qsnpgwa_out.hh" // this gives us access to the writeout format
#include "../utils/linear_regression.hh"

#include "../../logger/log.hh"

#if CONT_GENOSTATS_STATS
#include <iomanip> // reformat a few numbers.
#endif


class ContGenoStats{
	
	public :
		ContGenoStats(DataAccess *d);
		void prepGenoStatsForOutput(int snp, ContGenoStatsResults &results);

	protected :

		// used in output from the statistical hypothesis testing operation.
		struct statisticsOutput {
			double pVal;
			double beta;
			double se;
		};

		DataAccess *data;
		/* Mean of the response variable (phenotype) over all individuals.  Set in computeMeanAndSD.  */
		double responseMean; 
		bool ranMeans;
		double meanResidual; // holds mean of last set of residuals.
		vector<double> residuals; // Holds last set of residuals.
		
		
		double sumSquaresTotal(vector<double> response);

		// These are protected rather than public because their order is important.
		// See prepGenoStatsForOutput() implementation.
		/* Compute mean and stand. dev. statistics */
		void computeMeanAndSD(int snp, ContGenoStatsResults &results);
		/* Compute linear regression statistics */
		void computeLinRegStats(int snp, ContGenoStatsResults &results);
		/* Compute lack of fit stats */
		void computeLackOfFit(int snp, ContGenoStatsResults &results);
		/* Compute a single set of stats */
		statisticsOutput computeSingleStats(const vector<vector<double> > &input, const vector<double> &response, double meanResponse,
					string message, int snp);
		

		/* Adjust for covariates */
		void covariateAdjust(int snp);
		/* Zero out a ContGenoStatsResults struct's linear reg portion. */
		void blankLinRegStats(ContGenoStatsResults &results);
};

#endif
