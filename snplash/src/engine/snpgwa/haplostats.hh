//      haplostats.hh
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
 * Haplotype test statistics.
 */

#include "../em/em.h"
#include "../engine.h"
#include "../utils/statistics.h"
#include "../output/snpgwa_out.hh" // this gives us access to the writeout format
#include "../utils/zaykin.hh"

#define DEBUG_GLOBSTAT 0
#define DEBUG_HAPL_PROGRESS 0

using namespace std;

class HaploStats {
	
	public :
		HaploStats(DataAccess *, EngineParamReader *p);
		void prepHaploStatsForOutput(int, HaploStatsResults &);
		void setHaploThresh(int k){haploThresh = k;}
		
	protected :
	
		DataAccess *data;
		EngineParamReader *params;
	
		void calculateAllelic(int);
		void calculateTwoMarker(int, int);
		void calculateThreeMarker(int, int, int);
		
		double computeGlobalZaykinStatisic(EM &em, int size, double &testStat, int &degFreedom);
		double computeGlobalZaykinStatisic(EM &em, int size, double &testStat, int &degFreedom, 
						const vector<double> &phenotype, const vector<vector<double> > &cov, int snp);
			
		/* The following are computed in calculateAllelic */
		double allelicPval, allelicChiSq;
		int allelicDF;
		/* The following is computed in calculateTwoMarker */
		double twoMarkerPval, twoMarkerChiS;
		int twoMarkerDF;
		vector<double> twoMarkerCaseHapFreq, twoMarkerCntrlHapFreq, twoMarkerCmbdHapFreq;
		/* The following are computed in calculateThreeMarker */
		double threeMarkerPval, threeMarkerChiS;
		int threeMarkerDF;
		vector<double> threeMarkerCaseHapFreq, threeMarkerCntrlHapFreq, threeMarkerCmbdHapFreq;
		
		int haploThresh; /// Used by Zaykin.  Passed in at start.
	
		static double EPS () {return 0.00001;}
};
