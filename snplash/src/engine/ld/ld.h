/*
 *      ld.h
 *      
 *      Copyright 2009 Richard T. Guy <guyrt@guyrt-lappy>
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

#ifndef LINKAGE_D_H
#define LINKAGE_D_H

/**
 * @class LinkageDisequilibrium
 * 
 * Linkage Disequilibrium computation class.  Will perform diagonal or 
 * full model by calling process().  Will perform LD between specified pair 
 * of SNPs as well.  
 * 
 * @author Richard T. Guy (in present form)
 * @author Joshua Grab, Matt Steigert, Carl D. Langefeld
 */

#define RUN_IN_PARALLEL_LD 1

#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "../../param/engine_param_reader.h"
#include "../engine.h"
#include "../randwh.h"
#include "../em/em.h"
#include "../output/dprime_out.h"

using namespace std;

class LinkageDisequilibrium : public Engine{
	
	public : 
	
		explicit LinkageDisequilibrium();
		explicit LinkageDisequilibrium(DataAccess *);
		~LinkageDisequilibrium();
	
		virtual void init();
        virtual void preProcess();
        virtual void process();
        virtual void enslave(EngineParamReader *);
        virtual void test();
	
		bool dprimeOnPair(int s1, int s2, LinkageMeasures &results);
	
		// Computation engine pieces.  These are static so you could run
		// them without instantiating an LD engine if you already have
		// the EM results.
		static double computeDPrime(EM &, double aFreq[2][3]);
		static double dPrime(int iallele, int jallele, double aFreq[2][3], const vector<int> &, const vector<double> &);
		static double compDee(EM &e);
		static double compRSquare(double dee, double aFreq[2][3]);
		static double compDelta(double dee, EM &e);
	
	protected :
		
		bool haveOwner;
		int order_in_bag;
		int window;
		EngineParamReader *ld_param;
			
		bool diagonal;
		
		void delete_my_innards();
	
		LinkageOutput output;
		
		int numInitSNPs;
		int numInitPhen;
		int numFinalPhen;
	
	static double EPS(){return 0.00001;}
};

#endif

