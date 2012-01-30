//      cont_popstats.cpp
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

#include "cont_popstats.hh"

ContPopStats::ContPopStats(DataAccess *d){
	data = d;
	numTotal = pp = pq = qq = 0;
	missingQuant = nonMissingQuant = 0.0;
}

ContPopStats::~ContPopStats(){}

/*
 * Calculate the missing stats and the HWE stats.
 */
void ContPopStats::prepPopStatsForOutput(int snp, ContPopStatsResults &results){
	
	countIndividuals(snp);
	
	results.totalIndiv = pp+pq+qq;
	results.maf = (2*qq+pq) / (2*(pp+pq+qq));
	results.perMissing = 100 * ( 1 - (pp+pq+qq)/numTotal );
	
	results.numPP = pp; results.numPQ = pq; results.numQQ = qq;
	
	performMissingtTest(snp, results);
	performHWE(results);
}

/*
 * Count individuals, missing, and genotypes.
 */
void ContPopStats::countIndividuals(int snp){
	
	if(data->geno_size() <= snp){
		throw QSnpgwaException();
	}
	
	for(int i = 0;i<data->pheno_size() ;++i){
		numTotal++;
		
		if(data->get_data(i)->at(snp) == 0){
			missingQuant+=data->get_phenotype(i);
		}else{
			nonMissingQuant+=data->get_phenotype(i);
			switch(data->get_data(i)->at(snp)){
				case 1:
					pp++;
				break;
				case 2:
					pq++;
				break;
				case 3:
					pq++;
				break;
				case 4:
					qq++;
				break;
				default:
				//empty
				break;
			}
		}
	}	

}

/**
 * Perform a 2-tailed t-test of missingness.
 */
void ContPopStats::performMissingtTest(int snp, ContPopStatsResults &r){
	
	double meanMissing, meanNonMissing;
	int numMissing, numNonMissing;
	numNonMissing = pp+pq+qq;
	numMissing = numTotal - numNonMissing;
	
	if(numNonMissing * numMissing == 0){
		r.perMissingPVal = 2.0;
		return;
	}
	
	meanMissing = missingQuant / numMissing;
	meanNonMissing = nonMissingQuant / numNonMissing;
	
	double missingSampVar, nonMissingSampVar;
	missingSampVar = nonMissingSampVar = 0;
	// Get sample variance.
	for(int i = 0;i<data->pheno_size() ;++i){
		if(data->get_phenotype(i) != 0){
			int val = data->get_data(i)->at(snp);
			if(val == 0){
				missingSampVar += pow( (data->get_phenotype(i) - meanMissing  )   ,2);
			}else if(val < 5){
				nonMissingSampVar += pow( (data->get_phenotype(i) - meanNonMissing  )   ,2);
			}
		}
	}
	missingSampVar /= (numMissing - 1.0);
	nonMissingSampVar /= (numNonMissing - 1.0);
	
	double stdDev = sqrt((((numMissing - 1.0)*missingSampVar) + ((numNonMissing - 1.0)*nonMissingSampVar))/(numMissing + numNonMissing - 2.0));

	double tStat = (meanMissing - meanNonMissing)/(stdDev*sqrt((1.0/numMissing) + (1.0/numNonMissing)));
	
	int degFree;
	degFree = numTotal - 2;
	
	try{
		r.perMissingPVal = Statistics::tdist(tStat, degFree);
	}catch(StatsException){
		r.perMissingPVal = 2.0;
	}
}


/*
 * Perform HWE computations and store in the results.
 */
void ContPopStats::performHWE(ContPopStatsResults &results){
	
	int numTot = pp+pq+qq; // num with nonmissing.
	int numMajAlleles = 2*pp+pq;
	double majAlleleFreq = static_cast<double>(numMajAlleles)/(2*(pp+pq+qq));
	
	results.expPP = majAlleleFreq * majAlleleFreq*static_cast<double>(numTot);
	results.expPQ = 2*majAlleleFreq*(1-majAlleleFreq)*static_cast<double>(numTot);
	results.expQQ = (1-majAlleleFreq)*(1-majAlleleFreq)*static_cast<double>(numTot);
	
	if(results.expQQ > 0 && results.expPQ > 0 && results.expQQ > 0 ){
		double cmbdChiSquare = (pow(static_cast<double>(pp - results.expPP), 2)/results.expPP) +
		(pow(static_cast<double>(pq - results.expPQ), 2)/results.expPQ) +
		(pow(static_cast<double>(qq - results.expQQ), 2)/results.expQQ);
	
		try{
			results.chiSqPval = Statistics::chi2prob(cmbdChiSquare, 1);
		}catch(StatsException s){
			results.chiSqPval = 2.0;
		}
	}else{
		results.chiSqPval = 2.0;
	}
	results.pHWE = PopStats::exactTest(pp,pq,qq);		
	
}
