#ifndef LINKAGE_READER_H
#define LINKAGE_READER_H

/**
 *  This is a reader that will read a set of three files:
 * 		geno
 * 		phen
 * 		map
 */

#include "reader.h"
#include <map>		// Used to order individuals.

#define DBG_PROGRESS 0

class LinkageReader : public Reader {

public:
	LinkageReader();
	~LinkageReader();
	virtual void process(SnpData *, ParamReader *);

private:

	map<string, int> order_in_file;
	void getGenotype(SnpData *, ParamReader *);
	void getPhenotype(SnpData *, ParamReader *);
	void getMapFile(SnpData *, ParamReader *);
	
	void verify_phen_header(int, vector<int> *, vector<string> *, ParamReader *);

};

#endif
