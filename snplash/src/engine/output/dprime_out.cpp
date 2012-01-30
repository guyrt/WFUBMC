#include "dprime_out.h"


LinkageOutput::LinkageOutput(){
	outputType = 0;
	beginSNP = 0;
}
LinkageOutput::~LinkageOutput(){
	
}
bool LinkageOutput::init(int outputType, ParamReader *param, int maxMapSize){
	return init(outputType, param->get_out_file(), param, maxMapSize);
}
bool LinkageOutput::init(int outputType, string fileName, ParamReader *param, int maxMapSize){
	
	if(maxMapSize < 8)
		maxMapSize = 8;
	mapSize = maxMapSize;
	
	beginSNP = param->get_begin();
	
	this->outputType = outputType;
	bool ret = out.init(param->get_out_file());
	if(ret){
		
		out.write_header("**************************************************************************************\n");
		out.write_header("Drime, a part of SNPLASH Version ");
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
		
		out.write_header("\n\n*********************************************************\n\n");
        
        switch(this->outputType){
			case 2:
				// Do nothing.			
			break;
			case 3:
				for(int i=0;i < 2*mapSize; i++)
					out.write_header(" ");
				out.write_header("                Multi-\n");
				
				out.write_header("Marker 1");
				for(int i=0;i < mapSize-8; i++)
					out.write_header(" ");
				out.write_header("  ");
				out.write_header("Marker 2");
				for(int i=0;i < mapSize-8; i++)
					out.write_header(" ");	
				out.write_header("              allelic     Biallelic   Biallelic\n");
				out.write_header("  Name  ");
				for(int i=0;i < mapSize-8; i++)
					out.write_header(" ");
				out.write_header("  ");
				out.write_header("  Name  ");
				for(int i=0;i < mapSize-8; i++)
					out.write_header(" ");
				out.write_header("  D           D Prime     R Squared   Delta\n");
				for(int i=0;i < mapSize; i++)
					out.write_header("-");
				out.write_header("  ");
				for(int i=0;i < mapSize; i++)
					out.write_header("-");
				out.write_header("  ----------  ----------  ----------- ----------\n");
			break;
			default :
				cerr << "You shouldn't see this." << endl;
			break;
		}
        
        
       
	}
	return ret;
}

/*
 * Build the line requested and call output.
 */
void LinkageOutput::printLine(LinkageMeasures m, int order){
	if(outputType == 3){
		
		stringstream ss;
		ss << strnutils::spaced_string(m.name1,8) <<  strnutils::spaced_string(m.name2,8,2);
		ss << strnutils::spaced_number(m.dee,10,7,2) << strnutils::spaced_number(m.dPrime,10,8,2);
		ss << strnutils::spaced_number(m.rsquare,10,8,2) << strnutils::spaced_number(m.delta,10,7,2) << endl;
		out.write_line(ss.str(), order);
	
	}else if(outputType == 2){
		fmt2_storage[m.index1][m.index2] = m;
	}else{
		cerr << "Dprime output error: output type " << outputType << " unknown." << endl;
	}
}

/*
 * Build the output matrix.
 * 
 */
void LinkageOutput::create_fmt2_output(){
	
	int i_max = fmt2_storage.size()+1+beginSNP;
	
	
	stringstream ss;
	ss << "Marker-Marker D" << endl;
	ss << "      " ;
	for(int i=beginSNP;i<i_max;i++){
		ss << strnutils::spaced_number(i,8,1);
	}
	ss << endl;
	ss << "      " ;
	for(int i=beginSNP;i<i_max;i++){
		ss << " --------";
	}
	ss << endl;
	
	for(int i=beginSNP;i<i_max;i++){
		ss << strnutils::spaced_number(i,5,0);
		ss << " " ;
		for(int j=beginSNP;j<=i;j++)
			ss << "        .";
		for(int j=i+1;j<i_max;j++)
			ss << strnutils::spaced_number(fmt2_storage[i][j].dee,8,5,1);
		ss << endl;
	}
	ss << endl;
	
	// Next
	ss << "Marker-Marker Multiallelic D' (bounded between 0 and 1)" << endl;
	ss << "      " ;
	for(int i=beginSNP;i<i_max;i++){
		ss << strnutils::spaced_number(i,8,1);
	}
	ss << endl;
	ss << "      " ;
	for(int i=beginSNP;i<i_max;i++){
		ss << " --------";
	}
	ss << endl;
	
	for(int i=beginSNP;i<i_max;i++){
		ss << strnutils::spaced_number(i,5);
		ss << " " ;
		for(int j=beginSNP;j<=i;j++)
			ss << "        .";
		for(int j=i+1;j<i_max;j++)
			ss << strnutils::spaced_number(fmt2_storage[i][j].dPrime,8,5,1);
		ss << endl;
	}
	ss << endl;
	
	// Next
	ss << "Marker-Marker r^2" << endl;
	ss << "      " ;
	for(int i=beginSNP;i<i_max;i++){
		ss << strnutils::spaced_number(i,8,1);
	}
	ss << endl;
	ss << "      " ;
	for(int i=beginSNP;i<i_max;i++){
		ss << " --------";
	}
	ss << endl;
	
	for(int i=beginSNP;i<i_max;i++){
		ss << strnutils::spaced_number(i,5,0);
		ss << " " ;
		for(int j=beginSNP;j<=i;j++)
			ss << "        .";
		for(int j=i+1;j<i_max;j++)
			ss << strnutils::spaced_number(fmt2_storage[i][j].rsquare,8,5,1);
		ss << endl;
	}
	ss << endl;
	
	// Next
	ss << "Marker-Marker Delta" << endl;
	ss << "      " ;
	for(int i=beginSNP;i<i_max;i++){
		ss << strnutils::spaced_number(i,8,1);
	}
	ss << endl;
	ss << "      " ;
	for(int i=beginSNP;i<i_max;i++){
		ss << " --------";
	}
	ss << endl;
	
	for(int i=beginSNP;i<i_max;i++){
		ss << strnutils::spaced_number(i,5,0);
		ss << " " ;
		for(int j=beginSNP;j<=i;j++)
			ss << "        .";
		for(int j=i+1;j<i_max;j++)
			ss << strnutils::spaced_number(fmt2_storage[i][j].delta,8,5,1);
		ss << endl;
	}
	ss << endl;
	out.write_header(ss.str());
}


/*
 * Flush and close the output.
 */
void LinkageOutput::close(){
	
	if(outputType == 2){
		create_fmt2_output();
	}
	
	out.close();
}
