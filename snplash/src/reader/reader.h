#ifndef READER_H
#define READER_H

/**
 * Abstract reader class that contains all of the information 
 * that will be necessary for an actual reader.
 *  
 */
#include "../engine/snp_data.hh"
#include "../param/param_reader.h"
#include <string>
#include <ctype.h> // So we can get isspace.
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <limits>

using namespace std;

class Reader { 

public:
	virtual ~Reader();
	virtual void process(SnpData *, ParamReader *) = 0;
	bool getLine(ifstream *, vector<string> *, ParamReader *, bool check);
};

#endif
