#include "arff_reader.h"

ArffReader::ArffReader(){
	
}

/**
 * Process the arff file, placing it in snp data.  The parameter
 * reader passed in will contain the location of the file, ect.
 */
void ArffReader::process(SnpData *data, ParamReader *params){
	
	ifstream infile;
	infile.open(params->get_arff_file().c_str(), ifstream::in);
	if(! infile.is_open()){
		cerr << "Unable to open file " << params->get_arff_file() << endl;
		exit(1);
	}
	
	// Pass this ifstream to the map method then to the data method.
	
}
