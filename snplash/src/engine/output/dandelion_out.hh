//      dandelion_out.hh
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

#ifndef DANDELION_OUT_H
#define DANDELION_OUT_H

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
#include "../../param/param_reader.h"
#include "../../param/engine_param_reader.h"
#include <sstream> // Used to create and manage the string.
#include <map>
#include "../utils/stringutils.h"

using namespace std;

struct DandelionSnpInfo {

	int index;
	string name;
	string chr;
	int position;
	char majAllele;
	char minAllele;
	char refAllele;

};

struct DandelionHaploInfo {

	vector<char> alleles;
	double freqCs;
	double freqCn;
	double freqCb;
	double z;
	double p;
	double OR;
	double UCI;
	double LCI;

};


struct DandelionPProbInfo {

	int personNum;
	double affectionStatus;
	vector<char> leftHap;
	vector<char> rightHap;
	double prob;
	string personID;

};

class DandelionOutput {

	friend class Dandelion;

	public:
		DandelionOutput();
		~DandelionOutput();

		bool init(ParamReader *, EngineParamReader *, int numSnps, string message);
		void writeHeadSet(const vector<DandelionSnpInfo> &snps);
		void setMaxPersonId(int);
		void close();

	protected:
		Output outMain, outPProb;

		void writeLine();
		void writeLine(const DandelionHaploInfo &);
		void writeStatisticsLine(double pval, double chiSq, int DF);
		void writePProbLine(const DandelionPProbInfo &);

		void writeHeader(int);

		int personIdSpace;

};

#endif
