#ifndef ENGINE_H
#define ENGINE_H
/**
 * 
 * Abstract class for an engine.  An engine is a large algorithm class
 * like SNPGWA, INTERTWOLOG, ADTREE, or ADTREEBAG
 * 
 */
#include "../param/param_reader.h"
#include "../param/engine_param_reader.h"
#include "../reader/reader.h"
#include "../reader/linkage_reader.h"
#include "data_plugin.h"
#include <math.h>

class Engine {
	
public: 
	virtual void init() = 0;		// Do setup including *** data reading ***
	virtual void preProcess() = 0;  // Do any preprocessing that isn't done in constructor.
	virtual void process() = 0;		// This is the general driver for the algorithm.
	virtual void enslave(EngineParamReader *) = 0; // Allows another engine to control this one.
	virtual void test() = 0;		// Must include a tester which can do anything that we want.
	
	virtual ~Engine(){};
	
	// In general, a method must include each of these.  Reader must
	// be instantiated with a non-virtual reader.
protected:
	Reader *reader;
	ParamReader *param_reader; // Hold a pointer to the singleton parameter reader.
	DataAccess *data;	
	
};

#endif
