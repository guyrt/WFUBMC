//      dandelion.hh
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
 * Dandelion engine.
 *
 */

#ifndef DANDEL_H
#define DANDEL_H

#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include "../../param/engine_param_reader.h"
#include "../engine.h"
#include "../randwh.h"
#include "../em/em.h"
#include "../output/dandelion_out.hh"
#include "../utils/zaykin.hh"

#include "../ld/ld.h" // use this to make the r^2 output.

using namespace std;

class Dandelion : public Engine{

	public :

		explicit Dandelion();
		explicit Dandelion(DataAccess *);
		~Dandelion();

		virtual void init();
        virtual void preProcess();
        virtual void process();
        virtual void enslave(EngineParamReader *);
        virtual void test();

	protected :

		void delete_my_innards();
		vector<char> prepAlleles(unsigned int row, int divisor, int beg, int end);
		void computeHaplotypeTest(DandelionHaploInfo &d); // compute the Z stat, pval, OR, and CI for a single haplotype.

		void runSet(int begin, int end);

		EngineParamReader *ld_param;

		bool haveOwner;

		DandelionOutput output;

		int numInitSNPs;
		int numInitPhen;
		int numFinalPhen;
		int numCase;

	static double EPS() {return 0.00001;}
};

#endif

