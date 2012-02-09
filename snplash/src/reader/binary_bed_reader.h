#ifndef BINARY_BED_READER_H
#define BINARY_BED_READER_H

/**
 *  This is a reader that will read linkage formatted files.
 *  Expects there to be a linkage formatted parameter file.
 */

#include "reader.h"
#include <map>		// Used to order individuals.

using namespace std;

class BinaryBedReader : public Reader {

public:
	BinaryBedReader();
	~BinaryBedReader();
	virtual void process(SnpData *, ParamReader *);

private:

	map<string, int> order_in_file;
	void getBedFile(SnpData *, ParamReader *);
	void getPhenotype(SnpData *, ParamReader *);
	void getMapFile(SnpData *, ParamReader *);
	
	void verify_phen_header(int, vector<int> *, vector<string> *, ParamReader *);

	// Translate the Plink encoding to our encoding.
	inline 
	int binToCode(int t){
		// Note: I am doing an implicit bit reversal.
		// See http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml for encoding.
		if (t == 0){
			// homozygous major/major
			return 1;
		}else if (t == 2){
			// hetero
			return 2;
		}else if (t == 1){
			// missing
			return 0;
		}else if (t == 3){
			// homo minor/minor
			return 4;
		}else{
			// PROBLEM
			cerr << "Corrupt bed file. Can not read. Aborting." << endl;
			exit(0);
		}
	}

};

#endif
