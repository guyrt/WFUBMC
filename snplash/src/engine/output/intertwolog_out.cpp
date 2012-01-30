//      intertwolog_out.cpp
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


#include "intertwolog_out.hh"

InterTwoLogOutput::InterTwoLogOutput(){
}

InterTwoLogOutput::~InterTwoLogOutput(){
}

/*
 * Initialize all files and write their headers.
 * First, set some variables determining whether some things are even computed.
 *
 * Use the param object for that.
 *
 * @param param 		Parameter store
 * @param eparams		Engine parameter store.
 * @param maxMapSize	The max num characters in a map.
 * @param message		A string message to put in the output.
 * @return bool Successful init.
 */
bool InterTwoLogOutput::init(ParamReader *param, EngineParamReader *eparams, int maxMapSize, string message){

	bool ret = out.init(param->get_out_file());

	mapSize = maxMapSize;
	if(mapSize < 4) mapSize = 4;

	if(ret){
		writeMainHeader(param);
		out.write_header(message);
		writeMainLegend();
	}

	writeLogHead(param);

	return ret;
}

void InterTwoLogOutput::close(){
	out.close();
}

/*************************************************************************
 *
 * Main writing method.
 *
 *************************************************************************/
void InterTwoLogOutput::printLine(const InterTwoLogMeasures &itlo, int order){
	
	stringstream ss;
	ss << strnutils::spaced_number(itlo.index1, 8,0);
	ss << strnutils::spaced_string(itlo.name1, mapSize,2);
	ss << strnutils::spaced_number(itlo.index2, 8,2);
	ss << strnutils::spaced_string(itlo.name2, mapSize,2);
	ss << strnutils::spaced_number(itlo.beta, 10, 7,2);
	ss << strnutils::spaced_number(itlo.pVal, 10, 7,2);
	ss << strnutils::spaced_number(itlo.SE, 8, 6,2);
	ss << endl;
	
	out.write_line(ss.str(), order);
}



/*************************************************************************
 * 
 * Auxilary writers
 * 
 ************************************************************************/
void InterTwoLogOutput::writeMainHeader(ParamReader *param){
	
	out.write_header("**************************************************************************************\n");
	out.write_header("INTERTWOLOG, a part of SNPLASH Version ");
	out.write_header(SNPLASH_VERSION);
	out.write_header("\nDeveloped by: Linda E. Green, Stephanie R. Beck, Ethan M. Lange, \nMatt L. Stiegert, and Carl D. Langefeld.  Incorporated into SNPadt by Richard T. Guy\n");
	out.write_header("INPUT/OUTPUT PARAMETERS AND OPTIONS\n");
	out.write_header("***********************************\n");
	out.write_header("GENOTYPE FILE:          ");
	out.write_header(param->get_linkage_geno_file());
	out.write_header("\nPHENOTYPE FILE:         ");
	out.write_header(param->get_linkage_pheno_file());
	out.write_header("\nOUTPUT FILE:            ");
	out.write_header(param->get_out_file());

	out.write_header("\n\nTRAIT NAME:             ");
	string temp = param->get_trait();
	if(temp.compare("none") == 0){
		out.write_header("default");
	}else{
		out.write_header(temp);
	}

	out.write_header("\nCOVARIATE NAME(S):      ");
	vector<string> tempCov = param->get_covariates();
	if(tempCov.size() > 0){
		for(unsigned int i=0;i < tempCov.size();++i){
			out.write_header(tempCov.at(i));
			out.write_header("  ");
		}
	}else{
		out.write_header("none");
	}
	out.write_header("\n");

	out.write_header("\nOUTPUT NAMES:         ");
	out.write_header(param->get_out_file());
	
	out.write_header("\n********************************************************\n\n");
	
}


void InterTwoLogOutput::writeMainLegend(){
	stringstream ss;
	ss << "            ";
	for(int i=0;i < mapSize;++i) ss << " ";
	ss << "            ";
	for(int i=0;i < mapSize;++i) ss << " ";
	ss << "                        Interact." << endl;
	
	ss << "Marker 1  ";
	for(int i=0;i < mapSize;++i) ss << " ";
	ss << "  Marker 2  ";
	for(int i=0;i < mapSize;++i) ss << " ";
	ss << "  Interact.   Interact.   Standard" << endl;
	
	
	ss << "Index     Name";
	for(int i=0;i < mapSize-4;++i) ss << " ";
	ss << "  Index     Name";
	for(int i=0;i < mapSize-4;++i) ss << " ";
	ss << "  Beta        P Value     Error" << endl;
	
	ss << "--------  ";
	for(int i=0;i < mapSize;++i) ss << "-";
	ss << "  --------  ";
	for(int i=0;i < mapSize;++i) ss << "-";
	ss << "  ----------  ----------  --------" << endl; 
	out.write_header(ss.str());
}

/**
 * Write a header to the logging file.
 */
void InterTwoLogOutput::writeLogHead(ParamReader *param){
	
	stringstream ss;
	ss << "Log file for the SNPGWA engine of SNPADT." << endl << endl;
	
	time_t rawtime;
	struct tm *timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	
	ss << "Run at local time " << asctime(timeinfo) << endl;
	
	ss << "GENOTYPE FILE: " << param->get_linkage_geno_file() << endl;
	ss << "PHENOTYPE FILE: " << param->get_linkage_pheno_file() << endl;
	ss << "OUTPUT FILE:  " << param->get_out_file() << endl;
	ss << "MAP FILE:  " << param->get_linkage_map_file() << endl;
	
	
	ss << "TRAIT NAME:  ";
	string temp = param->get_trait();
	if(temp.compare("none") == 0){
		ss << "default" << endl;
	}else{
		ss << temp << endl;
	}

	ss << "COVARIATE NAME(S): " << endl;
	vector<string> tempCov = param->get_covariates();
	if(tempCov.size() > 0){
		for(unsigned int i=0;i < tempCov.size();++i){
			ss << tempCov.at(i) << "  ";
		}
		ss << endl;
	}else{
		ss << "none" << endl;
	}
	
	Logger::Instance()->writeLine(ss.str());
}
