/*
 *      output.h
 *      
 *      Copyright 2010 Richard T. Guy <guyrt@guyrt-lappy>
 *      
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 2 of the License, or
 *      (at your option) any later version.
 *      
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *      
 *      You should have received a copy of the GNU General Public License
 *      along with this program; if not, write to the Free Software
 *      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *      MA 02110-1301, USA.
 */

/**
 * An output engine used for storing results from asychronous computation
 * in an ordered way.  Lines are assumed to be ordered, which gives 
 */

#ifndef OUTPUT_H
#define OUTPUT_H

#include <iostream>
#include <string>
#include <map>		// to maintain the lines we have to print
#include <fstream> // for file output.
using namespace std;

class Output {
	
	public:
	
		Output();
	 
		bool init(string fileName);
		void close();
		void flush();
		void write_header(string head);
		void write_line(string line, int order);	
		
	private:
	
		ofstream outstream;
		string outfile;
		
		/* Next int contains the next line missing */
		int n_line;
		/* Next int contains the next line to print */
		int p_line;
		
		map<int, string> lines;
};

#endif
