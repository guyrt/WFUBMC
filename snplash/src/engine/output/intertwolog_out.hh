//      intertwolog_out.hh
//      
//      Copyright 2010 Richard T. Guy <guyrt7@wfu.edu>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#ifndef INTTWOLOG_OUT_H
#define INTTWOLOG_OUT_H

/**
 * Output class for dprime.  Uses an output class.  
 */
#include "output.h"
#include "../../param/param_reader.h"
#include "../../param/engine_param_reader.h"

using namespace std;
struct InterTwoLogMeasures{
	int		index1;
	string name1;
	string name2;
	int		index2;
	double 	pVal;
	double 	beta;
	double 	SE;
};
class InterTwoLogOutput {
	
	friend class InterTwoLog;
	
	public:
		InterTwoLogOutput();
		~InterTwoLogOutput();
		
		bool init(ParamReader *, EngineParamReader *, int maxMapSize, string message);
		void close();
	
	protected:
		Output out;

		void printLine(const InterTwoLogMeasures &itlo, int order);
		
		
		void writeMainHeader(ParamReader *param);
		void writeMainLegend();
		void writeLogHead(ParamReader *param);
		
		int mapSize;

};
#endif
