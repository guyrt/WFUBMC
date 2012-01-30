//      snpgwa_out.hh
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

#ifndef SNPGWA_OUT_H
#define SNPGWA_OUT_H

/**
 * Output class for SNPGWA program.
 * 
 * There are several output files for SNPGWA including (some optional)
 * 		- main file
 * 		- geno1,2,3
 * 		- 3 hwe files
 * 		- log file
 * 		- ref allele file
 * 
 * All but the log file are handled in this class.
 */
#include "output.h"
#include "snpinfo_out.hh"
#include "../../param/param_reader.h"
#include "../../param/engine_param_reader.h"
#include <sstream> // Used to create and manage the string.
#include <map>
#include "../utils/stringutils.h"
#include "../../logger/log.hh"

#include <time.h>


using namespace std;

/*
 * These structures hold the data that will be eventually outputted.
 */
struct PopStatsResults{
	int caseCount;
	int cntrlCount;
	double caseRefFreq;
	double cntrlRefFreq;
	double pMissingCombined;
	double pMissingCase;
	double pMissingCntrl;
	double missingPVal;
	double missingOR;
	int cntrlPP;
	int casePP;
	double expcntrlPP;
	double expcasePP;
	double expcmbdPP;
	int cntrlPQ;
	int casePQ;
	double expcntrlPQ;
	double expcasePQ;
	double expcmbdPQ;
	int cntrlQQ;
	int caseQQ;
	double expcntrlQQ;
	double expcaseQQ;
	double expcmbdQQ;
	
	double cmbdTestStat ;
    double caseTestStat; 
    double cntrlTestStat;
	
	double cmbdPVal;
	double casePVal;
	double cntrlPVal;
	double cmbdExactPVal ;
    double caseExactPVal; 
    double cntrlExactPVal;
};

struct GenoStatsResults{
	
	double twodegTestStat;
	double twodegPVal;
	
	double domTestStat;
	double domPVal;
	double domOR;
	double domLCI;
	double domUCI;
	double domSens;
	double domSpec;
	double domCStat;
	
	double addTestStat;
	double addPVal;
	double addOR;
	double addLCI;
	double addUCI;
	
	double addSensNNRN;
	double addSpecNNRN;
	double addSensNNRR;
	double addSpecNNRR;
	double addSensNRRR;
	double addSpecNRRR;
	double addCStat;
	
	double recTestStat;
	double recPVal;
	double recOR;
	double recLCI;
	double recUCI;
	double recSens;
	double recSpec;
	double recCStat;
	
	double lofTestStat;
	double lofPVal;
};

struct HaploStatsResults{
	
	double dprime;
	double rsquare;
	
	double allelicChiS;
	double allelicDF;
	double allelicPval;
	
	double twoMarkerChiS;
	double twoMarkerDF;
	double twoMarkerPval;
	
	double threeMarkerChiS;
	double threeMarkerDF;
	double threeMarkerPval;
	
	vector<double> twoMarkerCaseFreq;
	vector<double> threeMarkerCaseFreq;
	vector<double> twoMarkerCntrlFreq;
	vector<double> threeMarkerCntrlFreq;
};

class SnpgwaOutput {
	
	friend class Snpgwa;
		
	public:
		SnpgwaOutput();
		~SnpgwaOutput();
		
		bool init(ParamReader *, EngineParamReader *, int maxMapSize, int numSNPs, string message);
		void close();
	
	protected:
		Output outMain, outGeno1, outGeno2, outGeno3, outHWEcase, outHWEcntrl, outHWEcomb, outRefAllele;
		Output outHap1, outHap2, outHap3;
		Output outVal;
		
		int maxMapSize, totalNumSNPs;

		void writeLine(int ids, SnpInfo &, const PopStatsResults &p, const HaploStatsResults &h, const GenoStatsResults &g);

		/* Extra writers */
		void writeHWELine(int ids, const SnpInfo &, const PopStatsResults &p);
		void writeGenoLine(int idx, const SnpInfo &s, const PopStatsResults &p, const GenoStatsResults &g, const HaploStatsResults &h);
		void writeValLine(int ids, const SnpInfo &, const PopStatsResults &p, const HaploStatsResults &h, const GenoStatsResults &g);
		void writeHaploLine(int idx, const SnpInfo &s, const HaploStatsResults &p);

		/* Output file type options */
		bool writeGenoFiles, writeHWEFiles, writeRefFile, writeValFile, writeHaploFiles;
		/* Output file options */
		bool writeMainFileMap, writeMainFileHap;
	
		/*  header writers */
		void writeMainHeader(Output &,ParamReader *);
		void writeMainLegend(Output &);
		void writeHWELegend(Output &, Output &, Output &);
		void writeValLegend(Output &);
		void writeGenoLegend(Output &o1, Output &o2, Output &o3);
		void writeHaploLegend(Output &o1, Output &o2, Output &o3);
		void writeLogHead(ParamReader *param); // write to log file.

};

#endif
