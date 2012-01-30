/*
 *      cross_val.h
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
#ifndef CROSSVAL_H
#define CROSSVAL_H

#include "../engine.h"
#include "../adtree/adtree.h"
#include "../bagging/bagging.h"
//#include <omp.h> // So we can have OpenMP threads.
#include "../../param/param_reader.h"
#include <fstream>

/*
 * Cross-validation is a technique for judging the accuracy of a prediction rather than
 * for returning a single model (though with a change or two you could extract either, in
 * fact, MDR does.)  Here, the cross-validation number of engines are formed and used
 * as slaves.   
 * 
 */

class CrossValidation : public Classifies {
	
    public : 
    
	explicit CrossValidation();

	~CrossValidation(){};
	
	/// From engine
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
	bool haveOwner;
	int order_in_bag;
	EngineParamReader *cross_param;
	vector<bool> classification;
	vector<int> out_of_bag;
	vector<Classifies *> engines; // Holds the cast of engines.
	ofstream outstream;
	
	/// Result vector
	int results[4];
	
};

#endif
