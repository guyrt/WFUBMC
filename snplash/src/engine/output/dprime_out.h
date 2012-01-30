/*
 *      dprime_out.h
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

#ifndef DPRIME_OUT_H
#define DPRIME_OUT_H

/**
 * Output class for dprime.  Uses an output class.  
 */
#include "output.h"
#include "../../param/param_reader.h"
#include <sstream> // Used to create and manage the string.
#include <map>
#include "../utils/stringutils.h"

using namespace std;
struct LinkageMeasures{
	int		index1;
	int		index2;
	string	name1;
	int 	position1;
	string	chr1;
	string  name2;
	int 	position2;
	string	chr2;
	double	dPrime;
	double	dee;
	double	delta;
	double	rsquare;
};

class LinkageOutput {
	
	friend class LinkageDisequilibrium;
	

	
	public:
		LinkageOutput();
		~LinkageOutput();
		
		bool init(int, ParamReader *, int maxMapSize); // use default file name.
		bool init(int, string filename, ParamReader *, int maxMapSize); // use passed in file name.
		
		void close();
	
	protected:
		Output out;
		/* Output types: 3 [default] -> row is a SNP 
		 * 				 2 -> matrix of values. */
		int outputType;
		void printLine(LinkageMeasures m, int order);
		
		/* Used in output format 2 */
		map<int, map<int, LinkageMeasures> > fmt2_storage;
		void create_fmt2_output();
		
		int mapSize;
		int beginSNP;
	
};
#endif
