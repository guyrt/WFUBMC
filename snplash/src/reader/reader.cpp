#include "reader.h"

Reader::~Reader(){}

/** getLine()
 *
 * Returns a single line of arbitrary length from the file, splitting
 * on spaces to create an array.
 *
 * Clears the array before starting.
 *
 * @param infile An incoming file stream
 * @param line The address of the vector that should hold the tokenized line.
 * @param params A parameter object
 * @param check_used Query this location agains the parameter object's use/not use information
 * @return Result of file->eof()
 *
 */
bool Reader::getLine(ifstream *infile, vector<string> *line, ParamReader *params, bool check_used){

	string this_piece;
	
	getline(*(infile), this_piece);

	// Break the string into pieces.
	bool spaces = true; // true if last was space or at start.
	int start = 0;
	int finish = 0;

	line->clear();
	line->reserve(1000);
	int element = 1;
	for(unsigned int i=0; i < this_piece.length(); i++){
		if(isspace(this_piece.at(i))){
			if(!spaces){
				// We have finished a chunk.
				if(!check_used || params->use_this_column(element++)){ // put in if not checking or if checking and used.
					line->push_back(this_piece.substr(start, finish-start));
				}
			}
			spaces = true;
		}else{

			if(spaces){
				// Start of a new chunk.  Set start and finish.
				start = i;
				finish = i+1;
			}else{
				finish++;
			}
			spaces = false;
		}
	}
	// If we ended with spaces == false then there was a last elt to push.
	if(!spaces){
		if(!check_used || params->use_this_column(element++)){ // put in if not checking or if checking and used.
			line->push_back(this_piece.substr(start, finish-start));
		}
	}

	return infile->eof();
}

