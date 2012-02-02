//      snpgwa.h
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
 * @class Snpgwa
 *
 * Main entry point for SNPGWA, a comprehensive single-SNP association test
 * suite developed at by Carl Langefeld and group.
 * 
 * Performs tests with both categorical and quantitative phenotypes and with
 * and without covariates.
 * 
 * @author Richard T. Guy in present form
 * @author Joshua Grab
 * @author Matt Steigert
 * @author Carl D. Langefeld
 */

#ifndef SNPGWA_H
#define SNPGWA_H

#include "../../param/engine_param_reader.h"
#include "genostats.hh"
#include "popstats.hh"
#include "haplostats.hh"
#include "../engine.h"
#include "../output/snpgwa_out.hh"

#include "../../logger/log.hh"
#include "../ld/ld.h"

#define DB_V_SNPGWA 0 // CHANGE TO 1 TO GET A PLAY BY PLAY


using namespace std;

class Snpgwa : private Engine{
	
	public:
		
		explicit Snpgwa();
		
		~Snpgwa();
		
		virtual void init();		// Do setup including *** data reading ***
		virtual void preProcess();  // Do any preprocessing that isn't done in constructor.
		virtual void process();		// This is the general driver for the algorithm.
		virtual void enslave(EngineParamReader *); // Allows another engine to control this one.
		virtual void test();		// Must include a tester which can do anything that we want.
	
			
	protected:
		
		EngineParamReader *snp_param;
		
		void delete_my_innards();
		
		void initToZero(PopStatsResults &p, HaploStatsResults &r, GenoStatsResults &ge);

		SnpgwaOutput out;
		
		int numInitSNPs;
		int numInitPhen;
		int numFinalPhen;
	
};
#endif
