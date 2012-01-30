#ifndef ARFF_READER_H
#define ARFF_READER_H

/**
 *  This is a reader that will read linkage formatted files.
 *  Expects there to be a linkage formatted parameter file.
 */

#include "reader.h"

using namespace std;

class ArffReader : public Reader {

public:
	ArffReader();
	virtual void process(SnpData *, ParamReader *);

private:
	
	

};

#endif

