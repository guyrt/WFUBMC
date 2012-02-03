/*
 *      bagging.h
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

#ifndef BAGGING_H
#define BAGGING_H

/**
 * Bagging wrapper.  Utilizes another engine to analyze bootstrap samples. 
 * 
 * Author: Richard T. Guy
 */

#include "../combinable.h"
#include "../adtree/adtree.h"
#include "../../param/engine_param_reader.h"
#include "../../param/param_reader.h"
#include <fstream>


class Bagging : public Classifies {
	
	public : 
	
		explicit Bagging();
		explicit Bagging(DataAccess *);
		virtual ~Bagging();

		/// From Engine.h	
		virtual void init();
        virtual void preProcess();
        virtual void process();
        virtual void enslave(EngineParamReader *);
        virtual void test();
        
        /// Classification
        virtual void classify(int a[4]); // Run test on the actual data.
        virtual void classify(int a[4], SnpData *); // Run test on passed data.
		virtual int classify(vector<short> &); // Run test on single individual.
	
	private : 
	
		bool haveOwner; // Primarily used to deal with garbage collection.
		int order_in_bag;
		vector<ADTree> engines;
		
		ParamReader::EngineTypes engine_type;
		
		EngineParamReader *bag_param;
		
		bool evaluate_trees();
		bool evaluate_trees2();
		void print_and_process(int, AD_Rule &a);
		void print_and_process(int bags, vector<long> SNPs);
		void delete_my_innards();
		
		ofstream outstream;
	
};


#endif
