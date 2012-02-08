#include "ld.h"

LinkageDisequilibrium::LinkageDisequilibrium(){
	this->param_reader = ParamReader::Instance();
	this->data = new DataAccess;
	this->data->init(NULL);
	this->ld_param = new EngineParamReader;
	this->order_in_bag = 0;
	this->window = 0;
	this->haveOwner = false;
	diagonal = true;
}

LinkageDisequilibrium::LinkageDisequilibrium(DataAccess *data){
	this->param_reader = ParamReader::Instance();
	this->data = data;
	this->ld_param = new EngineParamReader;
	this->order_in_bag = 0;
	this->window = 0;
	this->haveOwner = false;
	diagonal = true;
}

LinkageDisequilibrium::~LinkageDisequilibrium(){
	if(!haveOwner) delete_my_innards();
}

void LinkageDisequilibrium::delete_my_innards(){
	
	this->param_reader = NULL;
	delete this->data;
	data = NULL;
	delete this->ld_param;
	ld_param = NULL;
}

void LinkageDisequilibrium::init(){
	ld_param->read_parameters(param_reader->get_engine_specific_params());

	if(param_reader->get_linkage_map_file().compare("none") == 0){
		cerr << "DPRIME requires that you enter a map file.  Aborting." << endl;
		//exit(0);
	}

	if(ld_param->get_dprime_smartpairs() && param_reader->get_linkage_map_file().compare("none") == 0){
		cerr << "Warning: Dprime parameter mismatch.  Smart pairs called, but no map file supplied." << endl;
		cerr << "All pairs will be computed, but the overhead will be slightly higher." << endl;
	}

	initializeReader(); // defined in engine.h

	reader->process(data->getDataObject(), param_reader);

	numInitSNPs = data->geno_size();
	numInitPhen = data->pheno_size();

	data->getDataObject()->remove_phenotype(numeric_limits<double>::max());
	data->getDataObject()->remove_covariate(numeric_limits<double>::max());
	data->getDataObject()->remove_phenotype(0);	
	data->getDataObject()->prep_data(param_reader);
	/*try{
		data->make_categorical(1,2);
	}catch(DataException d){
		cerr << "Data Exception with message: " << endl << d.message << endl;
		exit(0);
	}catch(...){
		cerr << "Exception making categorical." << endl;
	}*/
	delete reader; // No longer needed.
	
	numFinalPhen = data->pheno_size();
	
	
}

/*
 * Set up output.
 */
void LinkageDisequilibrium::preProcess(){
	
	int numCase = 0;
	for(int i=0;i < data->pheno_size(); ++i)
		if(equal(data->get_phenotype(i), 2)) numCase++;

	
	stringstream ss;
	
	if(!haveOwner){
		ss << "INDIVIDUALS READ FROM THE INPUT FILE:            " << numInitPhen << endl;
		ss << "INDIVIDUALS DELETED                              " << numInitPhen - numFinalPhen << endl;
		ss << "INDIVIDUALS LEFT FROM THE INPUT FILE:            " << numFinalPhen << "  (" << numCase << " cases and " << numFinalPhen - numCase << " controls)" << endl;
		
		cout << ss.str();
	}else{
		ss << "Dprime called as a subsidiary of another program." << endl;
	}
	
	if (haveOwner){
		if(! output.init(2,param_reader->get_out_file() + ".dprime" , param_reader,data->max_map_size())){
			cerr << "Bad setup.  Aborting." << endl;
			exit(0);
		}
	}else{
		if(! output.init(ld_param->get_dprime_fmt(),param_reader,data->max_map_size())){
			cerr << "Bad setup.  Aborting." << endl;
			exit(0);
		}	
	}
	
}

/**
 * Compute DPrime for all pairs of SNPs.
 * 
 * If the dp-smartpairs option was chosen, then only compute if the map chr matches.
 */
void LinkageDisequilibrium::process(){

	int run_size = data->geno_size();	
	// Set up window with window size.
	window = ld_param->get_dprime_window();
	if(window < 0) window = run_size + 1;
	int order = 0;
	int ceil;
	if(ld_param->get_dprime_smartpairs()){
		
		vector<int> run_pairs;
		for(int i=0;i<run_size;i++){
			// compute the pairs we will run.
			run_pairs.clear();
			ceil = run_size < i+window ? run_size : window+i;
			for(int j=i+1;j<ceil;j++){
				if(data->get_chrom(i) == data->get_chrom(j))	
					run_pairs.push_back(j);
			}

			#if RUN_IN_PARALLEL_LD
			#pragma omp parallel 
			{
			#endif
				int rp_sz = run_pairs.size();
				#if RUN_IN_PARALLEL_LD
				#pragma omp for schedule(static,1)
				#endif
				for(int ju=0;ju<rp_sz;ju++){
					LinkageMeasures l;
					if(data->getDataObject()->isUsable(i) && data->getDataObject()->isUsable(ju)){
						dprimeOnPair(i,ju, l);
					}else{
						l.dPrime = 0;
						l.dee = 0;
						l.delta = 0;
						l.rsquare = 0;
					}
					l.index1 = i+param_reader->get_begin();
					l.index2 = ju+param_reader->get_begin();
					
					data->get_map_info(l.index1-1, l.chr1, l.name1, l.position1);
					data->get_map_info(l.index2-1, l.chr2, l.name2, l.position2);
					#pragma omp critical
					{
						output.printLine(l, order+ju);
					}
				}
			#if RUN_IN_PARALLEL_LD
			}
			#endif
			order += run_pairs.size();
		}	
		
	}else{
		
		for(int i=0;i<run_size;i++){
			ceil = run_size < i+window ? run_size : i+window;
			#if RUN_IN_PARALLEL_LD
			#pragma omp parallel shared(ceil, order) 
			{
				#pragma omp for schedule(static)
			#endif
				for(int j=1+i;j<ceil;j++){
					LinkageMeasures l;
					if(data->getDataObject()->isUsable(i) && data->getDataObject()->isUsable(j)){
						dprimeOnPair(i,j, l);
					}else{
						l.dPrime = 0;
						l.dee = 0;
						l.delta = 0;
						l.rsquare = 0;
					}
					l.index1 = i+param_reader->get_begin();
					l.index2 = j+param_reader->get_begin();
					data->get_map_info(i, l.chr1, l.name1, l.position1);
					data->get_map_info(j, l.chr2, l.name2, l.position2);
					output.printLine(l, order+j-i-1);
				}
			#if RUN_IN_PARALLEL_LD
			}
			#endif
			order += ceil-i-1;
		}
		output.close();
	}
	
}

/**
 * Enslave this LD engine. For an LD engine, enslavement means that the output
 * file name has to change and the format becomes 2.
 * 
 */
void LinkageDisequilibrium::enslave(EngineParamReader *e){
	haveOwner = true;
}

void LinkageDisequilibrium::test(){
	// Keep this blank for now.
}


///////////////////// COMPUTATION ///////////////////////
/**
 * Perform all linkage disequilibrium measures on two SNPs
 * 
 * @param s1 SNP 1
 * @param s2 SNP 2
 * @return lr A filled linkage measures struct.
 * 
 * @bug The error spits to cerr not logger.
 */
bool LinkageDisequilibrium::dprimeOnPair(int s1, int s2, LinkageMeasures &lr){
	
	if(s1 == s2) return false;
	
	if(s1 >= data->geno_size() || s2 >= data->geno_size()){
		lr.dPrime = -1;
		lr.dee = -1;
		lr.rsquare = -1;
		lr.delta = -1;
		return false;
	}
	
	vector<int> numAlleles;
	vector<double> unknownProb;
	EM emAlgorithm;
	
	vector<short> v1, v2;
	for(int i=0; i<data->pheno_size(); i++ ){
		if(data->get_data(i)->at(s1) != 0 && data->get_data(i)->at(s2) != 0 && data->get_phenotype(i) != 0){
			// Data to be used only if not missing in both and not in phenotype
			v1.push_back( data->get_data(i)->at(s1) );
			v2.push_back( data->get_data(i)->at(s2) );
		}
	}
	
	if(v1.empty()){
		lr.dPrime = lr.dee = lr.rsquare = lr.delta = -1;
		return true;
	}
	
	vector<vector<short> > t;
	t.push_back(v1);
	t.push_back(v2);
	try{
		emAlgorithm.setup(t);
		emAlgorithm.run();
	}catch(...){
		// TODO: make this a better error.
		
		cerr << "Error running EM algorithm on SNPs " << s1+1 << " and " << s2+1 << "." << endl;
		lr.dPrime = -1;
		lr.dee = -1;
		lr.rsquare = -1;
		lr.delta = -1;
		return false;
	}
	double a[2][3];

	emAlgorithm.getAlleleFreqs(a);
	numAlleles = emAlgorithm.getNumberAlleles();
	
	// Analyze results of em.
	lr.dPrime = computeDPrime(emAlgorithm,a);
	lr.dee = 0, lr.rsquare = 0; lr.delta = 0;
	if((numAlleles[0] == 2) && (numAlleles[1] == 2)){
		lr.dee = compDee(emAlgorithm);
		lr.rsquare = compRSquare(lr.dee, a);
		lr.delta = compDelta(lr.dee,emAlgorithm);
	}

	return true;
}

/**
 * Compute DPrime, from which all other measures are computed.
 *
 * @param aFreq contains allele freqs in form A,a,B,b
 * @param emAlg Contains an EM object that includes data from previous run.
 * @return Dprime measurement.
 */
double LinkageDisequilibrium::computeDPrime(EM &emAlg, double aFreq[2][3]){
	int iallele, jallele;
	double total = 0;
	vector<int> numAlleles = emAlg.getNumberAlleles();
	vector<double> hapProbs = emAlg.getHapProbs();

	for(iallele = 1; iallele <= numAlleles[0]; iallele++){
		if(fabs(aFreq[0][iallele] - 1) < EPS()){
			return(-99.0);
		}
	}
	for(jallele = 1; jallele <= numAlleles[1]; jallele++){
		if(fabs(aFreq[1][jallele] - 1) < EPS()){
			return(-99.0);
		}
	}
	for(iallele = 1; iallele <= numAlleles[0]; iallele++){
		for(jallele = 1; jallele <= numAlleles[1]; jallele++){
			total += aFreq[0][iallele]*aFreq[1][jallele]*dPrime(iallele, jallele, aFreq, numAlleles, hapProbs);
		}
	}
	return(total);
}

/**
 * Called by computeDPrime.
 */
double LinkageDisequilibrium::dPrime(int iallele, int jallele, double aFreq[2][3], const vector<int> &numAlleles, const vector<double> &hapProbs)
{
	int tempVec[2] = {iallele, jallele};
	int hapIndex = 0;
	hapIndex = (numAlleles[1] + 1)*tempVec[0] + tempVec[1];

	double p = aFreq[0][iallele];
	double q = aFreq[1][jallele];

	double d11 = hapProbs[hapIndex] - p*q;

	double dmax = min(p*(1.0 - q), q*(1.0 - p));
	double dmin = max(-((1.0 - p)*(1.0 - q)), -(p*q));

	if((fabs(dmax) < EPS()) || (fabs(dmin) < EPS())){
		d11 = 0;
	}
	else if(d11 > 0){
		d11 /= dmax;
	}
	else if(d11 < 0){
		d11 /= dmin;
	}
	return(d11);
}

/**
 * Compute D
 */
double LinkageDisequilibrium::compDee(EM &e)
{
   double dee;
   vector<double> emFreqs = e.getEMFreqs();
   dee = emFreqs.at(0) * emFreqs.at(3) - emFreqs.at(1) * emFreqs.at(2);
   return dee;
}
/*
 * Compute r^2
 */
double LinkageDisequilibrium::compRSquare(double dee, double aFreq[2][3])
{
   double r, rSq;
   r = dee/sqrt(aFreq[0][1] * aFreq[0][2] * aFreq[1][1] * aFreq[1][2]);
   rSq = r*r;
   return rSq;
}
/*
 * Compute delta
 */
double LinkageDisequilibrium::compDelta(double dee, EM &e)
{
	
	vector<double> emFreqs = e.getEMFreqs();
	double delta = dee / ((emFreqs.at(0) + emFreqs.at(2)) * emFreqs.at(3));
	return delta;
	
}
