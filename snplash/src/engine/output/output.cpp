#include "output.h"

Output::Output(){
	n_line = p_line = 0;
	outfile = "";
	
}

bool Output::init(string fileName){
	outstream.open(fileName.c_str());
	return outstream.is_open();
}

/* Check that no lines remain in buffer and then print. */	
void Output::close(){
	if(lines.size() > 0){
		cerr << "Caution: calling close on an asynchronous output buffer"
			<< " with " << lines.size() << " unwritten lines."  << endl;
		cerr << "Line " << p_line << " was never received." << endl;
	}
	outstream.close();
	lines.clear();
}

void Output::write_header(string head){
	outstream << head;
}
/*
 * This is called to print a line.
 * 1) Include the line in the map 'lines' at specified order.
 * 2) If the order matches p_line then print it and all others while
 * they exist.
 */
void Output::write_line(string line, int order){
	#pragma omp critical
	{
		lines[order] = line;
		if(order == p_line){
			outstream << lines[p_line];
			lines.erase(p_line);
			p_line++;
			while(lines.count(p_line) > 0){
				outstream << lines[p_line];
				lines.erase(p_line);
				p_line++;
			}
		}else if(p_line > order){
			cerr << "Problem in output class: requested print of a line out of order.  Got " << p_line << " but expected >= " << order << endl;
		}
	}
}
