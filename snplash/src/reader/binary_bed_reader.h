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



/****
May need to include some of these (from Tommy, /home/guyrt7/School/Parallel/floyd.h)
to do reading of binary files with fread
probably just stdio and stdlib
(This shows entire text of floyd.h)

#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

//#define DBGIN
//#define DBGFWC
//#define DBGFWR
//#define DBGOUT

int colSize;
int rowSize;
int rootSize;
int rootProcs;

int id;
int r_id, c_id;
int r_limit_lo, r_limit_hi; // holds limits of computation.
int c_limit_lo, c_limit_hi; // inclusive on low, not on hi.
int nprocs;

int block_c, block_r; // Used to determine block sizes.  Contain proc's elements of blockList.

int *blockList;
float *data;

void readMatrix(char *);
void writeMatrix(char *);
void waitForMatrix();
void sendMatrixForWriting();

void performFW();

****/
