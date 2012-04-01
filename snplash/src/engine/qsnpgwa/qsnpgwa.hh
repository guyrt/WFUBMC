//      qsnpgwa.hh
//
//      Copyright 2010 Richard T. Guy <richardtguy84@gmail.com>
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
 * @class QSnpgwa
 *
 *
 * Perform single locus tests using a quantitative phenotype and zero or
 * more covariates.
 *
 * @author Richard T. Guy
 * @author Joshua Grab
 * @author Matt Steigert
 * @author Carl D. Langefeld
 *
 */

#ifndef QSNPGWA_H
#define QSNPGWA_H

// A 1 means spit out checkpoints (for infinite loop debugging)
#define DB_V_QSNP 0
// A

#include "../../param/engine_param_reader.h"
#include "../engine.h"
#include "../output/qsnpgwa_out.hh"
#include "cont_popstats.hh"
#include "cont_genostats.hh"
#include "../ld/ld.h"

using namespace std;

class QSnpgwa : public Engine {

	public:

		explicit QSnpgwa();

		~QSnpgwa();

		virtual void init();		// Do setup including *** data reading ***
		virtual void preProcess();  // Do any preprocessing that isn't done in constructor.
		virtual void process();		// This is the general driver for the algorithm.
		virtual void enslave(EngineParamReader *); // Allows another engine to control this one.
		virtual void test();		// Must include a tester which can do anything that we want.


	protected:

		EngineParamReader *snp_param;
		void delete_my_innards();

		void initToZero(ContPopStatsResults &p, ContGenoStatsResults &ge);

		QSnpgwaOutput out;

		int numInitSNPs;
		int numInitPhen;
		int numFinalPhen;

};

#endif
