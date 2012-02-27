//      qsnpgwa.cpp
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


#include "qsnpgwa.hh"

QSnpgwa::QSnpgwa(){
	param_reader = ParamReader::Instance();
	this->data = new DataAccess;
	this->data->init(NULL);
	this->snp_param = new EngineParamReader;
}


QSnpgwa::~QSnpgwa(){
	delete_my_innards();
}

void QSnpgwa::delete_my_innards(){
	delete this->param_reader;
	this->param_reader = NULL;
	delete this->data;
	data = NULL;
	delete this->snp_param;
	snp_param = NULL;
}

/*
 * Get data and populate the parameter holder.
 */
void QSnpgwa::init(){
	snp_param->read_parameters(param_reader->get_engine_specific_params());

	initializeReader(); // defined in engine.h

	if(param_reader->get_linkage_map_file().compare("none") == 0){
		cerr << "QSNPGWA requires that you enter a map file.  Aborting." << endl;
		exit(0);
	}

	Logger::Instance()->init(param_reader->get_out_file() + ".log");

	reader->process(data->getDataObject(), param_reader);
	
	numInitSNPs = data->geno_size();
	numInitPhen = data->pheno_size();

	delete reader; // No longer needed.
	data->getDataObject()->remove_phenotype(numeric_limits<double>::max());
	data->getDataObject()->remove_covariate(numeric_limits<double>::max());
	
	data->getDataObject()->prep_data(param_reader);
}

void QSnpgwa::preProcess(){

	stringstream ss;
	ss << "INDIVIDUALS READ FROM THE INPUT FILE:            " << numInitPhen << endl;
	int indivLeft = data->pheno_size();
	ss << "INDIVIDUALS LEFT FROM INPUT FILE:                " << indivLeft << endl;

	cout << ss.str();

	if(!out.init(param_reader, snp_param, data->max_map_size(),ss.str())){
		cerr << "Error opening files.  Aborting." << endl;
		exit(0);
	}
}

void QSnpgwa::process(){

	int sz = data->geno_size();
	
	// This is the primary loop.
	for(int i=0;i < sz;i++){
		
		#if DB_V_QSNP
			cout << i << " starting" << endl;
		#endif
		
		SnpInfo s;
		ContPopStatsResults p;
		ContGenoStatsResults g;
		s.index = i+param_reader->get_begin();

		#if DB_V_QSNP
			cout << i << " start setup" << endl;
		#endif
		data->get_map_info(i, s.chr, s.name, s.position);
		#if DB_V_QSNP
			cout << i << " end setup" << endl;
		#endif
		
		
		data->get_allele_codes(i, s.majAllele, s.minAllele, s.refAllele);
		if(s.majAllele == ' '){
			s.majAllele = '.';
		}
		if(s.minAllele == ' '){
			s.minAllele = '.';
		}
		if(s.refAllele == ' '){
			s.refAllele = '.';
		}
		
		// calc difference in genetic position.
		if(i < sz-1){

			string t, chr;
			int pos = 0;

			data->get_map_info(i+1, chr, t, pos);
			if(chr.compare(s.chr) == 0){
				s.diff = pos-s.position;
			}else{
				s.diff = -999;
			}

		}else{
			s.diff = -999;
		}

		// This if statement passes only if this is a useable snp.
		if(data->getDataObject()->isUsable(i)){
			
			ContPopStats pop_calc(data);
			ContGenoStats gen_calc(data);
	
			LinkageDisequilibrium ld(data);
			ld.enslave(snp_param); // VERY IMPORTANT.  We have to do this or the ld engine will think
									// that it owns the data and might will erase it.
									
			LinkageMeasures lr;
			
			#if DB_V_QSNP
			cout << i << " start pop_calc" << endl;
			#endif
			
			pop_calc.prepPopStatsForOutput(i,p);
			
			#if DB_V_QSNP
			cout << i << " end pop_calc" << endl;
			cout << i << " start gen_calc" << endl;
			#endif
			gen_calc.prepGenoStatsForOutput(i,g);
			
			#if DB_V_QSNP
			cout << i << " end gen_calc" << endl;
			cout << i << " start ld calc" << endl;
			#endif		
	
			if(i + 1 < data->geno_size()){
				ld.dprimeOnPair(i, i+1, lr);
				g.rsquare = lr.rsquare;
				g.dprime = lr.dPrime;
			}else{
				g.rsquare = g.dprime = -1.0;
			}
			
			#if DB_V_QSNP
			cout << i << " end ld calc" << endl;
			#endif
	
		}else{
			initToZero(p, g);
		}
		#if DB_V_QSNP
			cout << i << " write start" << endl;
		#endif
		out.writeLine(i, s,p,g);
		#if DB_V_QSNP
			cout << i << " write end" << endl;
		#endif
	}
	out.close();
}

// zero out.
void QSnpgwa::initToZero(ContPopStatsResults &p, ContGenoStatsResults &ge){
	p.totalIndiv = 0;
	p.maf = 0;
	p.perMissing = 0;
	p.perMissingPVal = 2.0;
	
	p.numPP = 0;
	p.numPQ = 0;
	p.numQQ = 0;
	
	p.expPP = 0;
	p.expPQ = 0;
	p.expQQ = 0;
	
	p.chiSqPval = -1;
	p.pHWE = 0;
	
	ge.twodegfree_pval = 2;
	
	ge.add_pval = 2;
	ge.add_beta = -1;
	ge.add_se = -1;  
	
	ge.dom_pval = 2;
	ge.dom_beta = -1;
	ge.dom_se = -1;
	
	ge.rec_pval = 2;
	ge.rec_beta = -1;
	ge.rec_se = -1;
	
	ge.meanAA = 0;
	ge.sdAA = 0;
	ge.meanAa = 0;
	ge.sdAa = 0;
	ge.meanaa = 0;
	ge.sdaa = 0;
	ge.meanAA_Aa = 0;
	ge.sdAA_Aa = 0;
	ge.meanAa_aa = 0;
	ge.sdAa_aa = 0;
	
	ge.rsquare = -1;
	ge.dprime = -1;
}

void QSnpgwa::enslave(EngineParamReader *){
	cerr << "Enslavement not supported for SNPGWA." << endl;
}

void QSnpgwa::test(){

}
