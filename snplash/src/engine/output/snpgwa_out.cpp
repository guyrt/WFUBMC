//      snpgwa_out.cpp
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


#include "snpgwa_out.hh"

SnpgwaOutput::SnpgwaOutput(){
	writeGenoFiles = writeHWEFiles = writeHaploFiles = false;
	writeMainFileMap = false;
	writeMainFileHap = writeRefFile = true;
	writeValFile = false;

	maxMapSize = 0;
}

SnpgwaOutput::~SnpgwaOutput(){
	close();
}

/**
 * Initialize all files and write their headers.
 * First, set some variables determining whether some things are even computed.
 *
 * Use the param object for that.
 *
 * @param param 		Parameter store
 * @param eparams		Engine parameter store.
 * @param maxMapSize	The maximum size of a map.
 * @param numSNPs		The number of SNPs to be counted.
 * @param message		A string message to include between header and legend on the
 * 						main screen.  Will include information on counts, ect.
 *
 * @return bool Successful init.
 */
bool SnpgwaOutput::init(ParamReader *param, EngineParamReader *eparams, int maxMapSize, int numSNPs, string message){

	if(maxMapSize < 6){
		this->maxMapSize = 6;
	}else{
		this->maxMapSize = maxMapSize;
	}

	totalNumSNPs = numSNPs;

	if(param->get_linkage_map_file().compare("none") != 0) writeMainFileMap = true;

	writeHWEFiles = eparams->get_snpgwa_dohap();

	writeGenoFiles = eparams->get_output_geno();
	writeHWEFiles = eparams->get_output_hwe();
	writeHaploFiles = eparams->get_output_haplo();

	writeValFile = eparams->get_output_val();

	bool ret = outMain.init(param->get_out_file());

	string t = param->get_out_file();

	if(writeGenoFiles){
		ret = ret && outGeno1.init(t + ".geno1");
		ret = ret && outGeno2.init(t + ".geno2");
		ret = ret && outGeno3.init(t + ".geno3");
	}

	if(writeHWEFiles){
		ret = ret && outHWEcase.init(t + ".hwecase");
		ret = ret && outHWEcntrl.init(t + ".hwecntrl");
		ret = ret && outHWEcomb.init(t + ".hwe");
	}

	if(writeHaploFiles){
		ret = ret && outGeno1.init(t + ".haplo1");
		ret = ret && outGeno2.init(t + ".haplo2");
		ret = ret && outGeno3.init(t + ".haplo3");
	}

	if(writeRefFile) ret = ret && outRefAllele.init(t + ".ref");

	if(writeValFile) ret = ret && outVal.init(t + ".statvals");

	if(ret){
		writeMainHeader(outMain, param);
		outMain.write_header(message);
		writeMainLegend(outMain);

		if(writeGenoFiles){
			writeMainHeader(outGeno1, param);
			writeMainHeader(outGeno2, param);
			writeMainHeader(outGeno3, param);
			writeGenoLegend(outGeno1, outGeno2, outGeno3);
		}
		if(writeHWEFiles){
			writeMainHeader(outHWEcase, param);
			writeMainHeader(outHWEcntrl, param);
			writeMainHeader(outHWEcomb, param);
			writeHWELegend(outHWEcase, outHWEcntrl, outHWEcomb);
		}
		if(writeHaploFiles){
			writeMainHeader(outHap1, param);
			writeMainHeader(outHap2, param);
			writeMainHeader(outHap3, param);
			writeHaploLegend(outHap1, outHap2, outHap3);
		}
		if(writeValFile){
			writeMainHeader(outVal, param);
			writeValLegend(outVal);
		}
	}

	writeLogHead(param);

	return ret;
}

void SnpgwaOutput::close(){
	outMain.close();
	outGeno1.close();
	outGeno2.close();
	outGeno3.close();
	outHWEcase.close();
	outHWEcntrl.close();
	outHWEcomb.close();
	outRefAllele.close();
	outVal.close();
	outHap1.close();
	outHap3.close();
	outHap2.close();
}

/*************************************************************************
 *
 * Main writing methods.
 *
 *************************************************************************/
/**
 * Main writing method.
 *
 * The structs are created and filled in snpgwa.cpp.
 * The index is for parallel creation of an in-order output file.
 * It must be 0 based.
 *
 * The main file is written here.  Others are written via calls from here.
 */
void SnpgwaOutput::writeLine(int idx, SnpInfo &s, const PopStatsResults &p,
	const HaploStatsResults &h, const GenoStatsResults &g)
{
	stringstream ss;

	s.print(ss, maxMapSize);

	ss << strnutils::spaced_number(p.caseCount, 8,1);
	ss << strnutils::spaced_number(p.cntrlCount, 8, 1);

	ss << strnutils::spaced_number(p.caseRefFreq, 8, 4, 4);
	ss << strnutils::spaced_number(p.cntrlRefFreq, 8, 4, 1);

	ss << strnutils::spaced_number(100*p.pMissingCombined, 8, 2, 4);
	ss << strnutils::spaced_number(100*p.pMissingCase, 8, 2, 1);
	ss << strnutils::spaced_number(100*p.pMissingCntrl, 8, 2, 1);
	ss << strnutils::spaced_number(p.missingPVal, 13, 10, 3);
	ss << strnutils::spaced_number(p.missingOR, 13, 10, 3);

	ss << "  " << s.majAllele << "  " << s.minAllele << "   " << s.refAllele;

	ss << strnutils::spaced_number(p.cntrlPP, 5, 2);
	ss << strnutils::spaced_number(p.casePP, 5, 1);
	if(p.expcmbdPP < 0.01){
		ss << strnutils::spaced_number(0, 8, 2, 1);
	}else{
		ss << strnutils::spaced_number(p.expcmbdPP, 8, 2, 1);	
	}
	ss << strnutils::spaced_number(p.cntrlPQ, 5, 1);
	ss << strnutils::spaced_number(p.casePQ, 5, 1);
	if(p.expcmbdPQ < 0.01){
		ss << strnutils::spaced_number(0, 8, 2, 1);
	}else{
		ss << strnutils::spaced_number(p.expcmbdPQ, 8, 2, 1);	
	}
	ss << strnutils::spaced_number(p.cntrlQQ, 5, 1);
	ss << strnutils::spaced_number(p.caseQQ, 5, 1);
	if(p.expcmbdQQ < 0.01){
		ss << strnutils::spaced_number(0.0, 8, 2, 1);
	}else{
		ss << strnutils::spaced_number(p.expcmbdQQ, 8, 2, 1);	
	}
	ss << strnutils::spaced_number(p.cmbdPVal, 12,10,1);
	ss << strnutils::spaced_number(p.cmbdExactPVal, 12,10,1);
	ss << strnutils::spaced_number(p.caseExactPVal, 12,10,1);
	ss << strnutils::spaced_number(p.cntrlExactPVal, 12,10,1);

	// insert geno information.
	ss << strnutils::spaced_number(g.twodegPVal, 12,10,4);
	ss << strnutils::spaced_number(g.domPVal, 12,10,1);
	ss << strnutils::spaced_number(g.domOR, 7,4,1);
	ss << strnutils::spaced_number(g.domLCI, 7,4,1);
	ss << strnutils::spaced_number(g.domUCI, 7,4,1);
	ss << strnutils::spaced_number(g.domSens, 6,4,1);
	ss << strnutils::spaced_number(g.domSpec, 6,4,1);
	ss << strnutils::spaced_number(g.domCStat, 6,4,1);

	ss << strnutils::spaced_number(g.addPVal, 12,10,1);
	ss << strnutils::spaced_number(g.addOR, 7,4,1);
	ss << strnutils::spaced_number(g.addLCI, 7,4,1);
	ss << strnutils::spaced_number(g.addUCI, 7,4,1);

	ss << strnutils::spaced_number(g.addSensNNRN, 6,4,1);
	ss << strnutils::spaced_number(g.addSpecNNRN, 6,4,1);
	ss << strnutils::spaced_number(g.addSensNNRR, 6,4,1);
	ss << strnutils::spaced_number(g.addSpecNNRR, 6,4,1);
	ss << strnutils::spaced_number(g.addSensNRRR, 6,4,1);
	ss << strnutils::spaced_number(g.addSpecNRRR, 6,4,1);
	ss << strnutils::spaced_number(g.addCStat, 6,4,1);

	ss << strnutils::spaced_number(g.recPVal, 12,10,1);
	ss << strnutils::spaced_number(g.recOR, 7,4,1);
	ss << strnutils::spaced_number(g.recLCI, 7,4,1);
	ss << strnutils::spaced_number(g.recUCI, 7,4,1);
	ss << strnutils::spaced_number(g.recSens, 6,4,1);
	ss << strnutils::spaced_number(g.recSpec, 6,4,1);
	ss << strnutils::spaced_number(g.recCStat, 6,4,1);

	ss << strnutils::spaced_number(g.lofPVal, 12,10,1);

	// LD info.

	ss << strnutils::spaced_number(h.dprime, 12,10,4);
	ss << strnutils::spaced_number(h.rsquare, 12, 10, 1);

	ss << strnutils::spaced_number(h.allelicPval, 12, 10, 4);

	ss << strnutils::spaced_number(h.twoMarkerPval, 12, 10, 4);
	double td;
	for(int i=0;i<4;i++){
		td = h.twoMarkerCaseFreq.at(i);
		if(td < 0.0001){
			ss << strnutils::spaced_number(0.0,6,4,1);
		}else{
			ss << strnutils::spaced_number(h.twoMarkerCaseFreq.at(i),6,4,1);
		}
	}
	for(int i=0;i<4;i++){
		td = h.twoMarkerCntrlFreq.at(i);
		if(td < 0.0001){
			ss << strnutils::spaced_number(0.0,6,4,1);
		}else{
			ss << strnutils::spaced_number(h.twoMarkerCntrlFreq.at(i),6,4,1);
		}
	}
	ss << strnutils::spaced_number(h.threeMarkerPval, 12, 10, 4);
	for(int i=0;i<8;i++){
		td = h.threeMarkerCaseFreq.at(i);
		if(td < 0.0001){
			ss  << strnutils::spaced_number(0.0,7,4,1);
		}else{
			ss  << strnutils::spaced_number(h.threeMarkerCaseFreq.at(i),7,4,1);
		}
	}
	for(int i=0;i<8;i++){
		td = h.threeMarkerCntrlFreq.at(i);
		if(td < 0.0001){
			ss << strnutils::spaced_number(0.0,7,4,1);
		}else{
			ss  << strnutils::spaced_number(h.threeMarkerCntrlFreq.at(i),7,4,1);
		}
	}

	ss << endl;

	outMain.write_line(ss.str(), idx);

	if(writeValFile) writeValLine(idx, s, p, h, g);
	if(writeHWEFiles) writeHWELine(idx, s, p);
	if(writeGenoFiles) writeGenoLine(idx, s, p, g, h);
	if(writeRefFile){
		stringstream ssr;
		ssr << s.chr << " ";
		if(s.name.size() > 0){
			ssr << strnutils::spaced_string(s.name, 10);
		}else{
			ssr << strnutils::spaced_number(s.index, 6);
		}
		ssr << " 0   " << s.position << "  " << s.refAllele << endl;
		outRefAllele.write_line(ssr.str(), idx);
	}

}
/**
 * Write output line to each of the three HWE files.
 * 
 * @param idx Index in final file (for synchronous update)
 * @param s Contains SNP info to write.
 * @param p Contains population stats results to write.
 */
void SnpgwaOutput::writeHWELine(int idx, const SnpInfo &s, const PopStatsResults &p){

	stringstream sscase, sscnt, ss;

	if(s.name.size() > 0){
		ss << strnutils::spaced_string(s.name, 10);
		sscase << strnutils::spaced_string(s.name, 10);
		sscnt << strnutils::spaced_string(s.name, 10);
	}else{
		ss << strnutils::spaced_number(s.index, 10);
		sscase << strnutils::spaced_number(s.index, 10);
		sscnt << strnutils::spaced_number(s.index, 10);
	}

	ss << strnutils::spaced_number(p.caseCount+p.cntrlCount, 9,1);
	sscase << strnutils::spaced_number(p.caseCount, 9,1);
	sscnt << strnutils::spaced_number(p.cntrlCount, 9,1);

	double freq = (static_cast<double>(p.caseCount) * p.caseRefFreq + static_cast<double>(p.cntrlCount) * p.cntrlRefFreq)/(p.caseCount+p.cntrlCount);
	ss << strnutils::spaced_number(freq, 11, 6, 4);
	sscase << strnutils::spaced_number(p.caseRefFreq, 11, 6, 4);
	sscnt << strnutils::spaced_number(p.cntrlRefFreq, 11, 6, 4);

	ss << strnutils::spaced_number(p.pMissingCombined, 9, 2, 4);
	sscase << strnutils::spaced_number(p.pMissingCase, 9, 2, 4);
	sscnt << strnutils::spaced_number(p.pMissingCntrl, 9, 2, 4);

	ss << "  " << s.majAllele << "  " << s.minAllele << "   " << s.refAllele;

	ss << strnutils::spaced_number(p.cntrlPP+p.casePP, 5,2);
	sscnt << strnutils::spaced_number(p.cntrlPP, 5,2);
	sscase << strnutils::spaced_number(p.casePP, 5,2);

	ss << strnutils::spaced_number(p.expcmbdPP, 8,2,1);
	sscnt << strnutils::spaced_number(p.expcntrlPP, 8,2,1);
	sscase << strnutils::spaced_number(p.expcasePP, 8,2,1);

	ss << strnutils::spaced_number(p.cntrlPQ+p.casePQ, 5,1);
	sscnt << strnutils::spaced_number(p.cntrlPQ, 5,1);
	sscase << strnutils::spaced_number(p.casePQ, 5,1);

	ss << strnutils::spaced_number(p.expcmbdPQ, 8,2,1);
	sscnt << strnutils::spaced_number(p.expcntrlPQ, 8,2,1);
	sscase << strnutils::spaced_number(p.expcasePQ, 8,2,1);

	ss << strnutils::spaced_number(p.cntrlQQ+p.caseQQ, 5,1);
	sscnt << strnutils::spaced_number(p.cntrlQQ, 5,1);
	sscase << strnutils::spaced_number(p.caseQQ, 5,1);

	ss << strnutils::spaced_number(p.expcmbdQQ, 8,2,1);
	sscnt << strnutils::spaced_number(p.expcntrlQQ, 8,2,1);
	sscase << strnutils::spaced_number(p.expcaseQQ, 8,2,1);

	ss << strnutils::spaced_number(p.cmbdPVal, 12,10,1);
	sscase << strnutils::spaced_number(p.casePVal, 12,10,1);
	sscnt << strnutils::spaced_number(p.cntrlPVal, 12,10,1);

	ss << strnutils::spaced_number(p.cmbdExactPVal, 12, 10, 1);
	sscase << strnutils::spaced_number(p.caseExactPVal, 12, 10, 1);
	sscnt << strnutils::spaced_number(p.cntrlExactPVal, 12, 10, 1);

	ss << endl;
	sscase << endl;
	sscnt << endl;

	outHWEcase.write_line(sscase.str(), idx);
	outHWEcntrl.write_line(sscnt.str(), idx);
	outHWEcomb.write_line(ss.str(), idx);
}

/**
 * 
 */
void SnpgwaOutput::writeHaploLine(int idx, const SnpInfo &s, const HaploStatsResults &p){
	
	stringstream ss1, ss2, ss3;
	
	
	
}

/**
 * Write output line for the value file.  This is a separate file that contains
 * test statistics rather than their corresponding p-values.
 * 
 * @param ids Index in final file (for asynchronous update)
 * @param s Contains SNP info to write.
 * @param p Contains population stats results to write.
 * @param h Contains haplotype stats results to write.
 * @param g Contains genotypic stats results to write.
 */
void SnpgwaOutput::writeValLine(int ids, const SnpInfo &s, const PopStatsResults &p, const HaploStatsResults &h, const GenoStatsResults &g){

	stringstream ss;

	if(s.name.size() > 0){
		ss << strnutils::spaced_string(s.name, 10);
	}else{
		ss << strnutils::spaced_number(s.index, 10);
	}

	ss << strnutils::spaced_number(p.cmbdTestStat, 14,12,1);
	ss << strnutils::spaced_number(p.caseTestStat, 14,12,1);
	ss << strnutils::spaced_number(p.cntrlTestStat, 14,12,1);

	ss << strnutils::spaced_number(g.twodegTestStat, 14,12,4);
	ss << strnutils::spaced_number(g.domTestStat, 14,12,1);
	ss << strnutils::spaced_number(g.addTestStat, 14,12,1);
	ss << strnutils::spaced_number(g.recTestStat, 14,12,1);
	ss << strnutils::spaced_number(g.lofTestStat, 14,12,1);

	ss << strnutils::spaced_number(h.allelicChiS, 14,12,4);
	ss << strnutils::spaced_number(h.allelicDF, 4,0,1);

	ss << strnutils::spaced_number(h.twoMarkerChiS, 14,12,4);
	ss << strnutils::spaced_number(h.twoMarkerDF, 4,0,1);
	ss << strnutils::spaced_number(h.threeMarkerChiS, 14,12,4);
	ss << strnutils::spaced_number(h.threeMarkerDF, 4,0,1);

	ss << endl;

	outVal.write_line(ss.str(), ids);

}

/**
 * Write lines for geno files.
 * 
 * @param idx Index in final file (for asynchronous update)
 * @param s Contains SNP info to write.
 * @param p Contains population stats results to write.
 * @param h Contains haplotype stats results to write.
 * @param g Contains genotypic stats results to write.
 */
void SnpgwaOutput::writeGenoLine(int idx, const SnpInfo &s, const PopStatsResults &p, const GenoStatsResults &g, const HaploStatsResults &h){

	stringstream ss1, ss2, ss3;

	if(s.name.size() > 0){
		ss1 << strnutils::spaced_string(s.name, 10);
		ss2 << strnutils::spaced_string(s.name, 10);
		ss3 << strnutils::spaced_string(s.name, 10);
	}else{
		ss1 << strnutils::spaced_number(s.index, 10);
		ss2 << strnutils::spaced_number(s.index, 10);
		ss3 << strnutils::spaced_number(s.index, 10);
	}

	ss1 << strnutils::spaced_number(p.caseRefFreq, 5, 3, 1);
	ss1 << strnutils::spaced_number(p.cntrlRefFreq, 5, 3, 1);
	ss1 << strnutils::spaced_number(p.cmbdExactPVal, 8,6,1);
	ss1 << strnutils::spaced_number(h.dprime, 5,3,1);
	ss1 << strnutils::spaced_number(h.rsquare, 5,3,1);

	ss1 << "  " << s.majAllele << "  " << s.minAllele << "   " << s.refAllele;

	ss2 << strnutils::spaced_number(p.caseRefFreq, 5, 3, 1);
	ss2 << strnutils::spaced_number(p.cntrlRefFreq, 5, 3, 1);
	ss2 << strnutils::spaced_number(p.cmbdExactPVal, 8,6,1);
	ss2 << strnutils::spaced_number(h.dprime, 5,3,1);
	ss2 << strnutils::spaced_number(h.rsquare, 5,3,1);

	ss2 << "  " << s.majAllele << "  " << s.minAllele << "   " << s.refAllele;

	ss3 << strnutils::spaced_number(p.caseRefFreq, 5, 3, 1);
	ss3 << strnutils::spaced_number(p.cntrlRefFreq, 5, 3, 1);
	ss3 << strnutils::spaced_number(p.cmbdExactPVal, 8,6,1);
	ss3 << strnutils::spaced_number(h.dprime, 5,3,1);
	ss3 << strnutils::spaced_number(h.rsquare, 5,3,1);

	ss3 << "  " << s.majAllele << "  " << s.minAllele << "   " << s.refAllele;

	// diverges here.
	ss1 << strnutils::spaced_number(g.twodegPVal, 12, 10, 1);
	ss1 << strnutils::spaced_number(g.domPVal, 12, 10, 1);
	ss1 << strnutils::spaced_number(g.addPVal, 12, 10, 1);
	ss1 << strnutils::spaced_number(g.recPVal, 12, 10, 1);
	ss1 << strnutils::spaced_number(g.lofPVal, 12, 10, 1);

	ss2 << strnutils::spaced_number(g.domOR, 7, 2, 1);
	ss2 << strnutils::spaced_number(g.domLCI, 7, 2, 1);
	ss2 << strnutils::spaced_number(g.domUCI, 7, 2, 1);
	ss2 << strnutils::spaced_number(g.addOR, 7, 2, 1);
	ss2 << strnutils::spaced_number(g.addLCI, 7, 2, 1);
	ss2 << strnutils::spaced_number(g.addUCI, 7, 2, 1);
	ss2 << strnutils::spaced_number(g.recOR, 7, 2, 1);
	ss2 << strnutils::spaced_number(g.recLCI, 7, 2, 1);
	ss2 << strnutils::spaced_number(g.recUCI, 7, 2, 1);

	ss3 << strnutils::spaced_number(g.domSpec, 4, 2, 1);
	ss3 << strnutils::spaced_number(g.domSens, 4, 2, 1);
	ss3 << strnutils::spaced_number(g.domCStat, 4, 2, 1);
	ss3 << strnutils::spaced_number(g.addSensNNRN, 4,2,1);
	ss3 << strnutils::spaced_number(g.addSpecNNRN, 4,2,1);
	ss3 << strnutils::spaced_number(g.addSensNNRR, 4,2,1);
	ss3 << strnutils::spaced_number(g.addSpecNNRR, 4,2,1);
	ss3 << strnutils::spaced_number(g.addSensNRRR, 4,2,1);
	ss3 << strnutils::spaced_number(g.addSpecNRRR, 4,2,1);
	ss3 << strnutils::spaced_number(g.addCStat, 4,2,1);
	ss3 << strnutils::spaced_number(g.recSpec, 4, 2, 1);
	ss3 << strnutils::spaced_number(g.recSens, 4, 2, 1);
	ss3 << strnutils::spaced_number(g.recCStat, 4, 2, 1);

}

/**************************************************************************
 * 																		  *
 * Header writers														  *
 * 																		  *
 **************************************************************************/
void SnpgwaOutput::writeMainHeader(Output &out, ParamReader *param){

	out.write_header("**************************************************************************************\n");
	out.write_header("SNPGWA, a part of SNPLASH Version ");
	out.write_header(SNPLASH_VERSION);
	out.write_header("\nDeveloped by: Linda E. Green, Stephanie R. Beck, Ethan M. Lange, \nMatt L. Stiegert, Richard T. Guy, and Carl D. Langefeld.  Incorporated into SNPadt by Richard T. Guy\n");
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

	if(writeValFile) out.write_header("   " + param->get_out_file() + ".statval");
	if(writeRefFile) out.write_header("   " + param->get_out_file() + ".ref");
	out.write_header("\n");

	int spaces=1; // number of extra lines to include.
	if(writeHWEFiles){
		out.write_header(param->get_out_file() + ".hwe      " + param->get_out_file() + ".hwecase  " + param->get_out_file() + ".hwecntrl\n");
	}else{
		spaces++;
	}
	if(writeGenoFiles){
		out.write_header(param->get_out_file() + ".geno1     " + param->get_out_file() + ".geno2     " + param->get_out_file() + ".geno3\n");
	}else{
		spaces++;
	}
	if(writeHaploFiles){
		out.write_header(param->get_out_file() + ".haplo1    " + param->get_out_file() + ".haplo2    " + param->get_out_file() + ".haplo3\n");
	}else{
		spaces++;
	}

	for(int i=0;i < spaces; i++) out.write_header("\n");

	out.write_header("START/END MARKER:              ");
	if(param->get_end() == INT_MAX){
		out.write_header("NOT CHOSEN\n");
	}else{
		stringstream ss;
		ss << param->get_begin() << " " << param->get_end() << endl;
		out.write_header(ss.str());
	}
	stringstream ss;

	ss << "   Begin Marker:               " << param->get_begin() << endl <<
		"   End Marker:                 " << param->get_begin() + totalNumSNPs -1 << endl;
	out.write_header(ss.str());
	out.write_header("\nOTHER OPTIONS\n**************************************************************************************\n");

	out.write_header("MAP FILE:                      ");
	if(param->get_linkage_map_file().compare("none") == 0){
		out.write_header("Not Chosen\n\n");
	}else{
		out.write_header(param->get_linkage_map_file() + "\n\n");
	}

	out.write_header("EXCL INDIV W/MISSING ALLELES:  ");
	if(param->get_ign()){
		out.write_header("Chosen\n\n");
	}else{
		out.write_header("Not Chosen\n\n");
	}

	out.write_header("**************************************************************************************\n");

}
/**
 * Write legend for the main lines.
 * The only decisions that have to be made are
 *  whether to print the map file info and whether to print haplo info.
 *
 * The spacing is very precise this method.
 * 
 * @param out Output class object to write to.  Utilizes asynchronized output.
 */
void SnpgwaOutput::writeMainLegend(Output &out){

	string mapFilePart1, mapFilePart2, mapFilePart3, mapFilePart4;
	string haploFilePart1, haploFilePart2, haploFilePart3, haploFilePart4;

	if(writeMainFileMap){

		stringstream m1, m2, m3, m4;
		m1 << "    ";
		m2 << "    ";
		m3 << "Chr ";
		m4 << "--- ";

		for(int i = 0; i < maxMapSize; i++){
			m1 << " ";
			m2 << " ";
			m4 << "-";
		}

		m3 << "Marker";
		for(int i = 0; i < maxMapSize-6; i++){
			m3 << " ";
		}

		m1 << "                       ";
		m2 << "                       ";
		m3 << " Position   Difference ";
		m4 << " ---------- ---------- ";

		mapFilePart1 = m1.str();
		mapFilePart2 = m2.str();
		mapFilePart3 = m3.str();
		mapFilePart4 = m4.str();

	}else{

		mapFilePart2 = mapFilePart1 = "           ";
		mapFilePart3 = "Marker     ";
		mapFilePart4 = "---------- ";
	}

	if(writeMainFileHap){
		haploFilePart1 = "    <------------------Two Marker Haplotype Analysis------------------->    <-----------------------------------------------------Three Marker Haplotype Analysis------------------------------------------------------>\n";
		haploFilePart2 = "                 <----------Cases----------> <--------Controls--------->                 <-----------------------------Cases---------------------------> <---------------------------Controls-------------------------->\n";
		haploFilePart3 = "    LRS_PValue   FreqPP FreqPQ FreqQP FreqQQ FreqPP FreqPQ FreqQP FreqQQ    LRS_PValue   FreqPPP FreqPPQ FreqPQP FreqPQQ FreqQPP FreqQPQ FreqQQP FreqQQQ FreqPPP FreqPPQ FreqPQP FreqPQQ FreqQPP FreqQPQ FreqQQP FreqQQQ\n";
		haploFilePart4 = "    ------------ ------ ------ ------ ------ ------ ------ ------ ------    ------------ ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------\n";
	}else{
		haploFilePart1 = haploFilePart2 = haploFilePart3 = haploFilePart4 = "\n";
	}

	out.write_header(mapFilePart1);
	out.write_header("                                                                                                                <---------------------------------------------Hardy-Weinberg Analysis-------------------------------------------->    <------------------------------------------------------------------------------------------------------Genotypic Association------------------------------------------------------------------------------------------------------>                                             ");
	out.write_header(haploFilePart1);

	out.write_header(mapFilePart2);
	out.write_header("   Individuals        Ref Allele Freq     <------------------  Percent Missing  ------------------->             Cntl  Case           Cntl  Case           Cntl  Case                       Combined     Case         Control         <-2 Deg Fr-> <------------------Dominant Test------------------------> <----------Additive Test-------------(NN v RN)-----(NN v RR)-----(NR v RR)----------> <--------------------Recessive Test---------------------> <-Lack Fit->                                 Allelic Test");
	out.write_header(haploFilePart2);

	out.write_header(mapFilePart3);
	out.write_header("Cases    Controls    Cases    Controls    Combined Cases    Controls   Pvalue          Odds Ratio    P  Q Ref     PP    PP   ENumPP    PQ    PQ   ENumPQ    QQ    QQ   ENumQQ  X^2_PValue   Prob_HWE     Prob_HWE     Prob_HWE        PV_2DF       PV_Dom       OR      LCI     UCI     Sens   Spec   C-St   PV_Add       OR      LCI     UCI     Sens   Spec   Sens   Spec   Sens   Spec   C-St   PV_Rec       OR      LCI     UCI     Sens   Spec   C-St   PV_LOF          D_Prime      R_Squared          PValue   ");
	out.write_header(haploFilePart3);

	out.write_header(mapFilePart4);
	out.write_header("-------- --------    -------- --------    -------- -------- --------   -------------   ------------- -- -- ---  ----- ----- -------- ----- ----- -------- ----- ----- -------- ------------ ------------ ------------ ------------    ------------ ------------ ------- ------- ------- ------ ------ ------ ------------ ------- ------- ------- ------ ------ ------ ------ ------ ------ ------ ------------ ------- ------- ------- ------ ------ ------ ------------    ------------ ------------    ------------");
	out.write_header(haploFilePart4);
}
/* Write all three legends for the hwe outputs. */
void SnpgwaOutput::writeHWELegend(Output &outCS, Output &outCN, Output &outCB){

	string cs1, cn1, cb1;

	cn1 = "            Number        Control       Control              <------------------------Hardy-Weinberg Analysis--------------------->\n";
	cs1 = "            Number          Case         Case                <------------------------Hardy-Weinberg Analysis--------------------->\n";
	cb1 = "            Number        Combined      Combined             <------------------------Hardy-Weinberg Analysis--------------------->\n";

	outCS.write_header(cs1);
	outCN.write_header(cn1);
	outCB.write_header(cb1);

	cn1 = "            Control      Reference      Percent                                                            Control      Control \n";
	cs1 = "             Case        Reference      Percent                                                              Case         Case  \n";
	cb1 = "            Case/Ctrl     Reference      Percent                                                            Combined     Combined\n";

	outCS.write_header(cs1);
	outCN.write_header(cn1);
	outCB.write_header(cb1);

	cn1 = "Marker      Indivs      Allele Freq     Missing    P  Q Ref  NumPP  ENumPP  NumPQ  ENumPQ  NumQQ  ENumQQ    X^2_PVal     Prob HWE\n";

	outCS.write_header(cn1);
	outCN.write_header(cn1);
	outCB.write_header(cn1);

	cn1 = "---------- ---------    -----------    ---------  -- -- ---  ----- -------- ----- -------- ----- -------- ------------ ------------\n";

	outCS.write_header(cn1);
	outCN.write_header(cn1);
	outCB.write_header(cn1);
}

void SnpgwaOutput::writeValLegend(Output &out){

	out.write_header("			 <---------Hardy-Weinberg Analysis---------->    <-------------------------Genotypic Association-------------------------->                             Two Marker            Three Marker\n");
	out.write_header("             Case/Cntl       Case          Control         <--2 Deg Fr--> <--Dom Test--> <--Add Test--> <--Rec Test--> <--Lack Fit-->     Allelic Test             Haplotype              Haplotype\n");
	out.write_header("Marker       X^2_Value      X^2_Value      X^2_Value          X^2Value       X^2Value       X^2Value       X^2Value       X^2Value          X^2Value    DegF       X^2Value    DegF       X^2Value    DegF\n");
	out.write_header("---------- -------------- -------------- --------------    -------------- -------------- -------------- -------------- --------------    -------------- ----    -------------- ----    -------------- ----\n");

}

/**
 * Write three legends.
 */
void SnpgwaOutput::writeGenoLegend(Output &o1, Output &o2, Output &o3){

	o1.write_header("            Reference                                 <--------------------Genotypic Association--------------------->\n");
	o1.write_header("           Allele Freq                                <-2 Deg Fr-> <-Dom Test-> <-Add Test-> <-Rec Test-> <-Lack Fit->\n");
	o1.write_header("Marker     Cases Ctrls Prob HWE D Pri R Squ  P  Q Ref PV 2DF       PV Dom       PV Add       PV Rec       PV LOF      \n");
	o1.write_header("---------- ----- ----- -------- ----- ----- -- -- --- ------------ ------------ ------------ ------------ ------------\n");

	o2.write_header("            Reference                                 <---------------------------Genotypic Test---------------------------->\n");
	o2.write_header("           Allele Freq                                <----Dominant Test----> <----Additive Test----> <---Recessive Test---->\n");
	o2.write_header("Marker     Cases Ctrls Prob HWE D Pri R Squ  P  Q Ref OR      LCI     UCI     OR      LCI     UCI     OR      LCI     UCI    \n");
	o2.write_header("---------- ----- ----- -------- ----- ----- -- -- --- ------- ------- ------- ------- ------- ------- ------- ------- -------\n");

	o3.write_header("            Reference                                 <--Dom Test--> <---------Additive Test----------> <--Rec Test-->\n");
	o3.write_header("           Allele Freq                                               (NN v RN) (NN v RR) (NR v RR)                    \n");
	o3.write_header("Marker     Cases Ctrls Prob HWE D Pri R Squ  P  Q Ref Spec Sens C-St Spec Sens Spec Sens Spec Sens C-St Spec Sens C-St\n");
	o3.write_header("---------- ----- ----- -------- ----- ----- -- -- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----\n");

}

/**
 * Write three legends
 */
void SnpgwaOutput::writeHaploLegend(Output &o1, Output &o2, Output &o3){

	o1.write_header("                        Reference\n");
	o1.write_header("                       Allele Freq                                        Allelic Test Two Marker   Three Marker\n");
	o1.write_header("Number  Chr Marker     Cases Ctrls    Prob HWE    D Pri R Squ   P  Q Ref     PValue    LRS PValue   LRS PValue  \n");
	o1.write_header("------- --- ---------- ----- -----    --------    ----- -----  -- -- ---  ------------ ------------ ------------\n");

	o2.write_header("                        Reference                          <--------Two Marker Haplotype Analysis-------->\n");
	o2.write_header("                       Allele Freq                         <--------Cases--------> <------Controls------->\n");
	o2.write_header("Number  Chr Marker     Cases Ctrls    Prob HWE   P  Q Ref  FrqPP FrqPQ FrqQP FrqQQ FrqPP FrqPQ FrqQP FrqQQ\n");
	o2.write_header("------- --- ---------- ----- -----    --------  -- -- ---  ----- ----- ----- ----- ----- ----- ----- -----\n");

	o3.write_header("                        Reference                          <-------------------------------Three Marker Haplotype Analysis------------------------------->\n");
	o3.write_header("                       Allele Freq                         <---------------------Cases-------------------> <-------------------Controls------------------>\n");
	o3.write_header("Number  Chr Marker     Cases Ctrls    Prob HWE   P  Q Ref  FrPPP FrPPQ FrPQP FrPQQ FrQPP FrQPQ FrQQP FrQQQ FrPPP FrPPQ FrPQP FrPQQ FrQPP FrQPQ FrQQP FrQQQ\n");
	o3.write_header("------- --- ---------- ----- -----    --------  -- -- ---  ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n");
}

/**
 * Write a header to the logging file.
 */
void SnpgwaOutput::writeLogHead(ParamReader *param){
	
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
	ss << endl << endl;
	
	Logger::Instance()->writeLine(ss.str());
}
