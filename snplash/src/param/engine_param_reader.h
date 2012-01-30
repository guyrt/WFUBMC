
/*
 *      engine_param_reader.h
 *
 *      Copyright 2009 Richard T. Guy <guyrt7@wfu.edu>
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

#ifndef ENG_PARAM_READER_H
#define ENG_PARAM_READER_H

#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "param_reader.h"

using namespace std;

class EngineParamReader {

	public :

		enum BaggingTypes { ALL, LEAVES, SUBTREE };

		EngineParamReader();

		void read_parameters(vector<string> *);

		/// Getters most of these should be const.
		int get_cross_validation() const {return number_of_crossvalidations;}
		int get_number_of_bags() const {return number_of_bags;}
		ParamReader::EngineTypes get_engine_type() const {return engine_type;}
		int get_threshold_of_bags() const {return threshold_of_bags;}
		double get_division_threshold() const {return division_threshold;}
		int get_number_of_nodes() const {return number_of_nodes;}
		BaggingTypes get_bag_type() const {return bag_type;}
		bool get_dprime_smartpairs() const {return dprime_smartpairs;}
		int get_dprime_fmt() const {return dprime_fmt;}
		int get_dprime_window() const {return dprime_window;}
		
		bool get_snpgwa_dohap() const {return snpgwa_dohaptest;}
		bool get_output_val() const {return output_val;}
		bool get_output_geno() const {return output_geno;}
		bool get_output_haplo() const {return output_haplo;}
		bool get_output_hwe() const {return output_hwe;}
		
		int get_haplo_thresh() const {return haplo_thresh;}
		
		bool get_dandelion_pprob() const {return dandelion_pprob;}
		const int get_dandelion_window() const {return dandelion_window;}
		
		const double getRegressionConditionNumberThreshold() const {return regression_condition_threshold;}
		
	protected :

		ParamReader::EngineTypes engine_type;

		// cv
		int number_of_crossvalidations;
		
		// bagging
		int number_of_bags;
		int threshold_of_bags;
		BaggingTypes bag_type;
		
		// ADTree
		int number_of_nodes;
		double division_threshold;
		
		// dprime
		bool dprime_smartpairs;
		int dprime_fmt;
		int dprime_window;
		
		//snpgwa and qsnpgwa
		bool snpgwa_dohaptest;
		/* formatting variable - 1 is just one file.  2 is all. */
		bool output_val; // if true, make output file with test stat values.

		bool output_geno, output_haplo, output_hwe;

		// Dandelion
		bool dandelion_pprob;
		int haplo_thresh; /// Used by zaykin (via snpgwa, dandelion) to set threshold of 
						  /// acceptance of haplotype for testing.
		int dandelion_window;
		
		/// Stats engines
		// This is a 1-norm threshold.
		double regression_condition_threshold;

		// Methods
		bool check_necessary();
		bool check_conflicting();
		bool check_engine_param_conflict();

};

#endif

