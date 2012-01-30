/*
 *      untitled
 *
 *      Copyright 2010 Richard T. Guy <guyrt7@wfu.edu>
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

#include "dandelion_out.hh"

DandelionOutput::DandelionOutput(){
	personIdSpace = 8;
}

DandelionOutput::~DandelionOutput(){

}

/**
 * Set personID.
 * 
 * @param int Size of person Id length.
 */
void DandelionOutput::setMaxPersonId(int maxPersonSize){
	personIdSpace = maxPersonSize;
}

/*
 * Open the output files and write the headers.
 */
bool DandelionOutput::init(ParamReader *param, EngineParamReader *engine_param, int numSnps, string message){

	bool ret = outMain.init(param->get_out_file());
	outMain.write_header(message);
	if(ret) writeHeader(numSnps);


	if(engine_param->get_dandelion_pprob()){
		ret = ret && outPProb.init(param->get_out_file() + ".pprob");
	}

	return ret;
}

/*
 * Write the header given the number of SNPs expected.
 */
void DandelionOutput::writeHeader(int numSnps){

	int numSpaces = 3*numSnps;
	if(numSpaces < 9) numSpaces = 9; // must be certain length.

	for(int i=0; i < numSpaces; i++){
		outMain.write_header(" ");
	}
	outMain.write_header("      <--EM Algorithm Est Haplotype Freqs--->    <---Haplo. Stat & PVal--->\n");

	for(int i=0; i < numSpaces-9; i++){
		outMain.write_header(" ");
	}
	outMain.write_header("Haplotype          Cases       Controls       Combined     ChiSq Stat        P-Value     Odds Ratio         95% CI\n");

	for(int i=0; i < numSpaces; i++){
		outMain.write_header("-");
	}
	outMain.write_header("     ----------     ----------     ----------    -----------     ----------     ----------   ------------\n");
}

/**
 * Write a header for a single series.
 * 
 * @param vector<DandelionSnpInfo> Set of SNPs.
 * 
 * Should look like:
 * 
 * Set of <num> SNPs
 *     rs123  ma  mi  chr  pos
 *     ....
 */
void DandelionOutput::writeHeadSet(const vector<DandelionSnpInfo> &snps){
	
	stringstream ss;
	ss << endl << endl << endl << "----------------------------------" << endl;
	ss << "Set of " << snps.size() << " SNPs" << endl;
	for(unsigned int i=0; i < snps.size(); ++i){
		ss << "    " << snps.at(i).name << "  " << snps.at(i).majAllele << "  " << snps.at(i).minAllele << "  " << snps.at(i).chr << "  " << snps.at(i).position << endl;
	}
	ss << "----------------------------------" << endl;
	outMain.write_header(ss.str());
}

/*
 * Write a formatted line given the input data.
 */
void DandelionOutput::writeLine(const DandelionHaploInfo &d){

	stringstream ss;

	if(d.alleles.size() == 1){
		ss << "      ";
	}else if(d.alleles.size() == 2){
		ss << "   ";
	}

	for(unsigned int i=0;i < d.alleles.size(); i++){

		ss << "  " << d.alleles.at(i);

	}

	ss << strnutils::spaced_number(d.freqCs, 10, 8, 5);
	ss << strnutils::spaced_number(d.freqCn, 10, 8, 5);
	ss << strnutils::spaced_number(d.freqCb, 10, 8, 5);

	if(d.z > 0){
		ss << strnutils::spaced_number(d.z, 10, 8, 5);
		ss << strnutils::spaced_number(d.p, 10, 8, 5);
		
		ss << strnutils::spaced_number(d.OR, 10, 4, 5);
		ss << "   (";
		ss << strnutils::spaced_number(d.LCI, 4, 2, 0);
		ss << ",";
		ss << strnutils::spaced_number(d.UCI, 4, 2, 1);
		ss << ")";
	}else{
		ss << "    -----------     ----------     ----------";
	}
	ss << endl;
	outMain.write_header(ss.str());
}

/*
 * Write a line detailing the statistics from Zaykins' global test.
 */
void DandelionOutput::writeStatisticsLine(double pval, double chiSq, int DF){
	stringstream ss;
	ss << endl << endl << endl;
	ss << "<----Likelihood Ratio Statistic---->" << endl;
	ss << "     Global  Degrees Of            " << endl;
	ss << "   Statistic     Freedom     P Value" << endl;
	ss << "------------  ----------  ----------" << endl;
    ss << strnutils::spaced_number(chiSq, 12, 8,0);
    ss << strnutils::spaced_number(DF, 11, 0, 2);
    ss << strnutils::spaced_number(pval, 10, 8, 2) << endl;
	outMain.write_header(ss.str());
}

/**
 * Write a line of output to the pprob file.
 * 
 * @param DandelionPProbInfo Holds the information to be written.
 */
void DandelionOutput::writePProbLine(const DandelionPProbInfo &p){

	stringstream ss;
//	ss << strnutils::spaced_number(p.personNum, personIdSpace,0);
	ss << strnutils::spaced_string(p.personID, personIdSpace, 1);
	ss << strnutils::spaced_number(p.affectionStatus, 9, 0, 2);
	ss << "   ";
	for(unsigned int i=0;i < p.leftHap.size(); i++){
		ss << "  " << p.leftHap.at(i);
	}
	ss << "   ";
	for(unsigned int i=0;i < p.rightHap.size(); i++){
		ss << "  " << p.rightHap.at(i);
	}
	ss << strnutils::spaced_number(p.prob, 11, 4, 4);
	ss << endl;
	outPProb.write_header(ss.str());

}


// Shut down and flush the out objects.
void DandelionOutput::close(){
	outMain.close();
}
