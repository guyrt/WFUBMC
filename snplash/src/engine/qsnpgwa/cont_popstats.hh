//      cont_popstats.hh
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
 * Contains all the programs for computing population statistics
 * for a continuous phenotype:
 * 		* % Missing
 * 		* HWE
 */

#include "../utils/exceptions.h"
#include <deque>
#include "../engine.h"
#include "../utils/statistics.h"
#include "../snpgwa/popstats.hh"
#include "../output/qsnpgwa_out.hh" // this gives us access to the writeout format

#ifndef CONT_POPSTATS_H
#define CONT_POPSTATS_H

using namespace std;

class ContPopStats{

	public:
		ContPopStats(DataAccess *d);
		~ContPopStats();

		/* Run all statistics */
		void prepPopStatsForOutput(int snp, ContPopStatsResults &);

	private:

		DataAccess *data;

		int numTotal;
		double pp, pq, qq;
		double missingQuant, nonMissingQuant;

		void countIndividuals(int snp);
		void performMissingtTest(int snp, ContPopStatsResults &);
		void performHWE(ContPopStatsResults &);

};

#endif
