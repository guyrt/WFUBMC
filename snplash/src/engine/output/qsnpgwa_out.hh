/*
 *      untitled
 *      
 *      Copyright 2010 Richard T. Guy <guyrt7@wfu.edu>
 *      
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 2 of the License, or
 *      (at your option) any later version.
 *      
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *      
 *      You should have received a copy of the GNU General Public License
 *      along with this program; if not, write to the Free Software
 *      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *      MA 02110-1301, USA.
 */

#ifndef QSNPGWA_OUT_H
#define QSNPGWA_OUT_H

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

using namespace std;

struct ContPopStatsResults{
	
	int totalIndiv;
	double maf;
	double perMissing;
	double perMissingPVal;
	
	int numPP;
	int numPQ;
	int numQQ;
	
	double expPP;
	double expPQ;
	double expQQ;
	
	double chiSqPval;
	double pHWE;
	
};

struct ContGenoStatsResults{
	
	double twodegfree_pval;
	
	double add_pval;
	double add_beta;
	double add_se;
	
	double dom_pval;
	double dom_beta;
	double dom_se;
	
	double rec_pval;
	double rec_beta;
	double rec_se;
	
	double lof_pval;
	
	double meanAA;
	double sdAA;
	double meanAa;
	double sdAa;
	double meanaa;
	double sdaa;
	double meanAA_Aa;
	double sdAA_Aa;
	double meanAa_aa;
	double sdAa_aa;
	
	double rsquare;
	double dprime;
};

class QSnpgwaOutput {
	
	friend class QSnpgwa;
	
	public :
	
		QSnpgwaOutput();
		~QSnpgwaOutput();
		
		bool init(ParamReader *param, EngineParamReader *eparams, int maxMapSize, string message);
		void close();
	
	protected :

		Output outMain, outGeno1, outGeno2, outGeno3, outHWE, outRefAllele;

		void writeLine(int idx, SnpInfo &q, const ContPopStatsResults &cp,
			const ContGenoStatsResults &cg);

		/* HWE writers */
		void writeHWELine(int ids);

		/* Output file type options */
		bool writeGenoFiles, writeHWEFiles, writeRefFile;
		/* Output file options */
		bool writeMainFileMap, writeMainFileHap;
	
		/*  header writers */
		void writeMainHeader(Output &,ParamReader *);
		void writeMainLegend(Output &);
		void writeHWELegend(Output &);
		
		void writeHWELine(int idx, const SnpInfo &s, const ContPopStatsResults &p);
		
		int mapSize;
};

#endif
