#ifndef CLASSIFY_H
#define CLASSIFY_H
/**
 * 
 * Abstract class for an engine that can be used in
 * a bag or under cross-validation.  
 * 
 * The point here is to include the capacity for classifying instances.
 * 
 * An engine is a large algorithm class
 * like SNPGWA, INTERTWOLOG, ADTREE, or ADTREEBAG
 * 
 */
#include "snp_data.hh"
#include "engine.h"

class Classifies : public Engine {
	
public: 

	virtual void classify(int a[4]) = 0; // Run test on the actual data.
    virtual void classify(int a[4], SnpData *) = 0; // Run test on passed data.
    virtual int classify(vector<short> &) = 0; // Run test on single instance.
	
};

#endif

