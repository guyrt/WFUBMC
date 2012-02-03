//      qsnpgwa_out.cpp
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

#include "qsnpgwa_out.hh"
#include "../../snplashConfig.h"

QSnpgwaOutput::QSnpgwaOutput(){
	writeGenoFiles = writeHWEFiles = writeRefFile = false;
	writeMainFileMap = false;
	writeMainFileHap = true;
}

QSnpgwaOutput::~QSnpgwaOutput(){
}

void QSnpgwaOutput::close(){
	outMain.close();
	outGeno1.close();
	outGeno2.close();
	outGeno3.close();
	outHWE.close();
	outRefAllele.close();
}

/*
 * Initialize all files and write their headers.
 * First, set some variables determining whether some things are even computed.
 * 
 * Use the param object for that.
 */
bool QSnpgwaOutput::init(ParamReader *param, EngineParamReader *eparams, int maxMapSize, string message){

	mapSize = maxMapSize;
	if(mapSize < 6) mapSize = 6;

	if(param->get_linkage_map_file().compare("none") != 0) writeMainFileMap = true;
	
	writeHWEFiles = eparams->get_output_haplo();
	writeGenoFiles = eparams->get_output_geno();

	bool ret = outMain.init(param->get_out_file());

	string t = param->get_out_file();

	if(writeGenoFiles){
		ret = ret && outGeno1.init(t + ".geno1");
		ret = ret && outGeno2.init(t + ".geno2");
		ret = ret && outGeno3.init(t + ".geno3");
	}

	if(writeHWEFiles){
		ret = ret && outHWE.init(t + ".hwe");
	}

	if(writeRefFile) ret = ret && outRefAllele.init(t + ".ref");

	if(ret){
		writeMainHeader(outMain, param);
		outMain.write_header(message);
		writeMainLegend(outMain);

		if(writeGenoFiles){
			writeMainHeader(outGeno1, param);
			outGeno1.write_header("specific");
			writeMainHeader(outGeno2, param);
			outGeno2.write_header("specific");
			writeMainHeader(outGeno3, param);
			outGeno3.write_header("specific");
		}
		if(writeHWEFiles){
			writeMainHeader(outHWE, param);
			writeHWELegend(outHWE);
		}
		if(writeRefFile){
			writeMainHeader(outRefAllele, param);
			outRefAllele.write_header("specific");
		}
	}

	return ret;
}

/**
 * Write the header for the main output.
 */
void QSnpgwaOutput::writeMainHeader(Output &out,ParamReader *param){
	
	out.write_header("**************************************************************************************\n");
	out.write_header("Qsnpgwa, a part of SNPLAST Version ");
	out.write_header(SNPLASH_VERSION);
	out.write_header("\nDeveloped by: Richard T. Guy, Matt L. Stiegert, Carl D. Langefeld, and Joshua D. Grab\n");
	out.write_header("Incorporated into SNPadt by Richard T. Guy\n");
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

	out.write_header("\nCOVARIATE NAME(S):             ");
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
}

/**
 * Write the main output legend.
 */
void QSnpgwaOutput::writeMainLegend(Output &out){
	
	out.write_header("\n\n   "); 
	for(int i=0;i < mapSize;i++)
		out.write_header(" ");
	out.write_header("                          ");
	out.write_header("Total   Minor                                      <------------------------Hardy-Weinberg Analysis------------------------>    <------------------------------------------------------------------------------------------------------Quantitative Genotypic Association------------------------------------------------------------------------------------------------------>\n");
	
	// line 2
	out.write_header("   ");
	for(int i=0;i < mapSize;i++)
		out.write_header(" ");
	out.write_header("                          ");
	out.write_header("Number  Allele                                                                                                                  <-2 Deg Fr-> <---------Dominant Test----------> <---------Additive Test----------> <---------Recessive Test---------> <-Lack Fit-> <---------------------------------------Means and Standard Deviations--------------------------------------->\n");	
	
	// line 3
	out.write_header("Chr Marker");
	for(int i=0;i < mapSize-6;i++)
		out.write_header(" ");
	out.write_header(" Position   Difference   ");
	out.write_header("Indivs  Freq    %_Miss %_Miss PVal    P  Q Ref(a)  NumPP   ENumPP  NumPQ   ENumPQ  NumQQ   ENumQQ  X^2_PValue   Prob_HWE        PV_2DF       PV_Dom       Dom_Beta   Dom_SE     PV_Add       Add_Beta   Add_SE     PV_Rec       Rec_Beta   Rec_SE     PV_LOF       Mean_AA    SD_AA      Mean_Aa    SD_Aa      Mean_aa    SD_aa      Mean_AA&Aa SD_AA&Aa   Mean_Aa&aa SD_Aa&aa      D_Prime R_Sqr\n");
	
	out.write_header("--- ------");
	for(int i=0;i < mapSize-6;i++)
		out.write_header("-");
	out.write_header(" ---------- ----------   ");	
	out.write_header("------  ------  ------ ------------  -- -- ------  ------- ------- ------- ------- ------- ------- ------------ ------------    ------------ ------------ ---------- ---------- ------------ ---------- ---------- ------------ ---------- ---------- ------------ ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------    ------- -------\n");
	
	
}

/**
 * Write legend for HWE
 */
void QSnpgwaOutput::writeHWELegend(Output &out){
	
out.write_header("           Total  Minor                                   <---------------Hardy-Weinberg Analysis---------------->\n");
out.write_header("           Number Allele                                                                                          \n");
out.write_header("Marker     Indivs Freq   %_Miss %_Miss PVal   P  Q Ref(a) NumPP ENumPP NumPQ ENumPQ NumQQ ENumQQ X^2_PVal Prob_HWE\n");
out.write_header("---------- ------ ------ ------ ------------ -- -- ------ ----- ------ ----- ------ ----- ------ -------- --------\n");
	
}

/**
 * Write a line.
 */
void QSnpgwaOutput::writeLine(int idx, SnpInfo &q, const ContPopStatsResults &cp
	, const ContGenoStatsResults &cg){
		
	stringstream ss;
	// print starting information.
	q.print(ss, mapSize);
	
	ss << strnutils::spaced_number(cp.totalIndiv, 6,3);
	ss << strnutils::spaced_number(cp.maf,6,4,2);
	
	ss << strnutils::spaced_number(cp.perMissing,6,2,2);
	ss << strnutils::spaced_number(cp.perMissingPVal,12,10,1);
	
	ss << "  " << q.majAllele << "  " << q.minAllele << "   " << q.refAllele;
	
	// HWE
	
	ss << strnutils::spaced_number(cp.numPP,7,6);
	ss << strnutils::spaced_number(cp.expPP,7,1,1);
	ss << strnutils::spaced_number(cp.numPQ,7,1);
	ss << strnutils::spaced_number(cp.expPQ,7,1,1);
	ss << strnutils::spaced_number(cp.numQQ,7,1);
	ss << strnutils::spaced_number(cp.expQQ,7,1,1);
	ss << strnutils::spaced_number(cp.chiSqPval,12,10,1);
	ss << strnutils::spaced_number(cp.pHWE,12,10,1);
	
	// Regression
	
	ss << strnutils::spaced_number(cg.twodegfree_pval,12,10,4);
	
	ss << strnutils::spaced_number(cg.dom_pval,12,10,1);
	ss << strnutils::spaced_number(cg.dom_beta,10,5,1);
	ss << strnutils::spaced_number(cg.dom_se,10,5,1);
	
	ss << strnutils::spaced_number(cg.add_pval,12,10,1);
	ss << strnutils::spaced_number(cg.add_beta,10,5,1);
	ss << strnutils::spaced_number(cg.add_se,10,5,1);
	
	ss << strnutils::spaced_number(cg.rec_pval,12,10,1);
	ss << strnutils::spaced_number(cg.rec_beta,10,5,1);
	ss << strnutils::spaced_number(cg.rec_se,10,5,1);
	
	ss << strnutils::spaced_number(cg.lof_pval, 12, 10, 1);

	// Moments
	
	ss << strnutils::spaced_number(cg.meanAA,10,4,1);
	ss << strnutils::spaced_number(cg.sdAA,10,4,1);
	ss << strnutils::spaced_number(cg.meanAa,10,4,1);
	ss << strnutils::spaced_number(cg.sdAa,10,4,1);
	ss << strnutils::spaced_number(cg.meanaa,10,4,1);
	ss << strnutils::spaced_number(cg.sdaa,10,4,1);
	ss << strnutils::spaced_number(cg.meanAA_Aa,10,4,1);
	ss << strnutils::spaced_number(cg.sdAA_Aa,10,4,1);
	ss << strnutils::spaced_number(cg.meanAa_aa,10,4,1);
	ss << strnutils::spaced_number(cg.sdAa_aa,10,4,1);
	
	// LD
	ss << strnutils::spaced_number(cg.dprime, 7, 5, 4);
	ss << strnutils::spaced_number(cg.rsquare, 7, 5, 1);
	
	ss << endl;
	outMain.write_line(ss.str(), idx);
	
	if(writeHWEFiles) writeHWELine(idx, q, cp);
}

void QSnpgwaOutput::writeHWELine(int idx, const SnpInfo &s, const ContPopStatsResults &p){
	
	stringstream ss;
	
	if(s.name.size() > 0){
		ss << strnutils::spaced_string(s.name, 10);
	}else{
		ss << strnutils::spaced_number(s.index, 10);
	}
	
	ss << strnutils::spaced_number(p.totalIndiv, 9,1);
	ss << strnutils::spaced_number(p.maf, 11, 6, 4);
	
	ss << strnutils::spaced_number(100*p.perMissing, 9, 2, 4);
	ss << strnutils::spaced_number(p.perMissingPVal, 12, 10, 4);
	
	ss << " " << s.majAllele;
	ss << " " << s.minAllele;
	
	ss << "     " << s.refAllele;
	
	ss << strnutils::spaced_number(p.numPP, 5,2);
	ss << strnutils::spaced_number(p.expPP, 5,2);
	
	ss << strnutils::spaced_number(p.numPQ, 5,2);
	ss << strnutils::spaced_number(p.expPQ, 5,2);
	
	ss << strnutils::spaced_number(p.numQQ, 5,2);
	ss << strnutils::spaced_number(p.expQQ, 5,2);
	
	ss << strnutils::spaced_number(p.chiSqPval, 8,6);
	ss << strnutils::spaced_number(p.pHWE, 8,6);
}
