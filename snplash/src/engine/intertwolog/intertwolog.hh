//      intertwolog.hh
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
 * @class InterTwoLog 
 *
 * @author Richard T. Guy
 *
 *
 * Perform adjusted logistic regression search for interaction using the following model:
 * 
 * logit(y) = beta_0 + beta * covariates + beta_{n-2} * snp1 + beta_{n-1} * snp2 + beta_n * (snp1-mu1) * (snp2-mu2)
 * 
 * Test H_0: beta_n = 0 using Wald test.
 *  
 */

#ifndef INTERTWOLOG_H
#define INTERTWOLOG_H

// Very verbose: primary use is for hangups.
#define INTERTWOLOG_TESTLOOP  0
// Turn on to get a glimpse of statistics for all elements of model (instead of just intercept)
#define INTERTWOLOG_TESTSTATS 0
// Turn on to get a peak at the raw data before it enters the LR engine.
#define INTERTWOLOG_DATAPEAK 0


#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include "../../param/engine_param_reader.h"
#include "../engine.h"
#include "../output/intertwolog_out.hh"
#include "../utils/lr.hh"

#include "../../logger/log.hh"

using namespace std;

class InterTwoLog : public Engine{

	public :

		explicit InterTwoLog();
		
		~InterTwoLog();

		virtual void init();
		virtual void preProcess();
		virtual void process();
		virtual void enslave(EngineParamReader *);
		virtual void test();

	protected :

		EngineParamReader *itl_param;
		bool haveOwner; // unused.
		
		InterTwoLogOutput out;
		
		int numInitSNPs;
		int numInitPhen;
		int numFinalPhen;
		int numCase;

		void delete_my_innards();
		void processPair(int i, int j, InterTwoLogMeasures &itlm);
};

#endif

