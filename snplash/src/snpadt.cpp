/*
 *      snpadt.cpp
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
// IDEA: add an include file that incorporates all of these.
#include "reader/linkage_reader.h"
#include "engine/machineLearning/adtree.h"
#include "engine/machineLearning/bagging.h"
#include "engine/machineLearning/cross_val.h"
#include "engine/utils/statistics.h"
#include "engine/qsnpgwa/qsnpgwa.hh"
#include "engine/snpgwa/snpgwa.h"
#include "engine/dandelion/dandelion.hh"
#include "engine/intertwolog/intertwolog.hh"
#include "engine/ld/ld.h"
#include "snplashConfig.h"
#include <iostream>
#include <string>
#include <typeinfo>
//#include "logger/log.hh"

using namespace std;

int main(int argc, char * argv[]){
	
	cout << endl <<  "Snplash, a statistical genetics tool from the Wake Forest University Health Sciences." << endl;
	cout << "Written by Richard T. Guy, Matt Stiegert, Josh D. Grab, and Carl D. Langefeld." << endl;
	cout << "This version compiled on " << __TIMESTAMP__ << " as version " << SNPLASH_VERSION << endl << endl << endl;

	// initialize the global random generator.
	srand ( time(NULL) );
	
	ParamReader *params = ParamReader::Instance();
	
	
	bool _continue = true;
	_continue = params->process_parameters(argc, argv);

	if(!_continue){
		delete params;
		return 0;
	}

	if(params->get_engine_types() == ParamReader::ADTREE){
	
		ADTree *adt_engine = new ADTree();
		
		adt_engine->init();
		adt_engine->preProcess();
		adt_engine->process();
		adt_engine->print_tree_to_file();
		delete adt_engine;
		
	}else if(params->get_engine_types() == ParamReader::BAGGING){
		Bagging *bag_engine = new Bagging();
		bag_engine->init();
		bag_engine->preProcess();
		bag_engine->process();
		delete bag_engine;
	}else if(params->get_engine_types() == ParamReader::CROSSVAL){
		
		CrossValidation *c_engine = new CrossValidation();
		c_engine->init();
		c_engine->preProcess();
		c_engine->process();
		int *results = new int[4];
		c_engine->classify(results);
		
		cout << "Results: " << endl;
		cout << results[0] << " " << results[1] << endl;
		cout << results[2] << " " << results[3] << endl;
		delete[] results;
		delete c_engine;
	}else if(params->get_engine_types() == ParamReader::DPRIME){
		
		LinkageDisequilibrium *d = new LinkageDisequilibrium();
		d->init();
		d->preProcess();
		d->process();
		delete d;
	}else if(params->get_engine_types() == ParamReader::SNPGWA){

			Snpgwa *d = new Snpgwa();
			d->init();
			d->preProcess();
			d->process();
			delete d;
	}else if(params->get_engine_types() == ParamReader::QSNPGWA){
		
		QSnpgwa *q = new QSnpgwa();
		q->init();
		q->preProcess();
		q->process();
		delete q;
	}else if(params->get_engine_types() == ParamReader::DANDELION){
		
		Dandelion *q = new Dandelion();
		q->init();
		q->preProcess();
		q->process();
		delete q;
	}else if(params->get_engine_types() == ParamReader::INTERTWOLOG){
		
		InterTwoLog *q = new InterTwoLog();
		q->init();
		q->preProcess();
		q->process();
		delete q;
	}else{
		cerr << "Engine type unknown." << endl;
	}

	Logger::Instance()->close();

	// Do not delete the new called above.  Once it is passed to an engine, that engine 
	// is responsible for cleanup.  See ~Bagging.
	return 0;
}
