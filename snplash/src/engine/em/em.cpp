#include "em.h"
using namespace std;

EM::EM(): numHaps(0), numSplits(4), setFor(0){}

EM::~EM(){
	for(unsigned int i=0;i < incompleteTable.size();i++)
		delete incompleteTable.at(i);
}

/**
 * Perform the initilization steps to run on an arbitrary number of
 * individuals.
 *
 * @param &data Each inner vector is a SNP.
 */
void EM::setup(const vector<vector<short> > &data){
	setFor = data.size();

	numberAlleles = data.size();

	inTotal = data;

	numAlleles.clear();
	for(unsigned int i=0; i < data.size(); i++)
		numAlleles.push_back(uniqueElements(data.at(i)));
//		numAlleles.push_back(2);

	numHaps = 1;
	for(unsigned int i=0; i < numAlleles.size(); i++)
		numHaps *= (numAlleles.at(i) + 1);

	unknownProb.clear();
	unknownProb.reserve(numHaps);
	tempCount.clear();
	tempCount.reserve(numHaps);
	emFrequencies.clear();
	emFrequencies.reserve(numHaps);
	for(int i=0;i < numHaps;i++){
		unknownProb.push_back(0);
		tempCount.push_back(0);
		emFrequencies.push_back(0);
	}

	numSplits = 1 << data.size();//pow(2, data.size());

	haplotype *nextHaplotype;
	for(int ihap = 0; ihap < numHaps; ihap++){
		nextHaplotype = new haplotype(ihap, numAlleles);
		incompleteTable.push_back(nextHaplotype);
	}

	/* Initialize haplotype1 and haplotype2. */
	haplotype1.clear();
	haplotype t;
	haplotype1.resize(data.at(0).size(),t);
	haplotype2.resize(data.at(0).size(),t);
	/* These require that we break our storage back apart. */
	vector<int> t1, t2;

	for(unsigned int i = 0; i < data.at(0).size(); i++){
		for(unsigned int j = 0; j < data.size(); j++){

			switch(data.at(j).at(i)){
				case 0:
					t1.push_back(0);
					t2.push_back(0);
				break;
				case 1:
					t1.push_back(1);
					t2.push_back(1);
				break;
				case 2:
					t1.push_back(1);
					t2.push_back(2);
				break;
				case 3:
					t1.push_back(2);
					t2.push_back(1);
				break;
				case 4:
					t1.push_back(2);
					t2.push_back(2);
				break;
				default:

				break;

			}
		}
		haplotype1.at(i).init(t1);
		haplotype2.at(i).init(t2);
		haplotype1.at(i).setHaplotypeIndex(numAlleles);
		haplotype2.at(i).setHaplotypeIndex(numAlleles);
		t1.clear();
		t2.clear();
	}
}


/**
 * Execute the EM algorithm for the individuals that have been loaded.
 */
bool EM::run(){

	#if DBG_EM
	cout << "EM run start" << endl;
	#endif
	if(setFor < 1) throw EMAlgorithmNoSetup();

	initEMAlgorithm();

	int ihap;
	double diff, epsilon = 0.000001;
	double *oldProb = NULL;
	try{
		oldProb = new double[numHaps];
	}
	catch(bad_alloc &memoryAllocationException){
		delete [] oldProb;
		cout << "Exception occurred: "
			<< memoryAllocationException.what() << endl;
		return false;
	}

	do{

		for(ihap = 0; ihap < numHaps; ++ihap){
			oldProb[ihap] = unknownProb.at(ihap);
		}
		computeIncompleteFreqs();
		cycleThroughPeople();
		distributeCounts();
		computeCompleteFreqs();

		diff = 0;
		for(ihap = 0; ihap < numHaps; ++ihap){
			diff += fabs(oldProb[ihap] - unknownProb.at(ihap));
		}

	}
	while(diff > epsilon);
	delete [] oldProb;

	#if DBG_EM
	cout << "EM run finalize 1" << endl;
	#endif
	computeIncompleteFreqs(); // so that loglike test can be run efficiently
	#if DBG_EM
	cout << "EM run finalize 2" << endl;
	#endif
	reweightKnownHaps();
	#if DBG_EM
	cout << "EM run finalize 3" << endl;
	#endif

	for(ihap = 0; ihap < numHaps; ihap++){
		emFrequencies.at(ihap) = unknownProb.at(ihap);
	}

	#if DBG_EM
	cout << "EM run end" << endl;
	#endif

	return true;
}

/*
 * Return the allele frequencies computed by the EM algorithm.
 * They are placed in the double array passed in.
 *
 * @input double** aFreq[2][2]
 * @input vector<int> numAlleles gets the interal number alleles array.
 * @input vector<double> hapProbs gets vector of haplo probabilities directly.
 */
void EM::result(double **aFreq, vector<int> &numAlleles, vector<double> &hapProbs){

	#if DBG_EM
	cout << "EM result start" << endl;
	#endif

	if(setFor != 2) throw EMAlgorithmRunResultsMismatch();

	int tempVec[2];
	int imark;
	tempVec[0] = tempVec[1] = 0;
	aFreq[0][0] = aFreq[0][1] = aFreq[1][0] = aFreq[1][1] = 0;

	for(imark = 0; imark < 2; imark++){
		int hapIndex, iallele;
		for(iallele = 1; iallele <= numAlleles.at(imark); iallele++){
			tempVec[imark] = iallele;
			hapIndex = (numAlleles.at(1) + 1)*tempVec[0] + tempVec[1];
			aFreq[imark][iallele] = unknownProb.at(hapIndex);
		}
		tempVec[imark] = 0;
	}
	numAlleles = this->numAlleles;
	hapProbs.clear();
	for(int i=0;i<numHaps;i++){
		hapProbs.push_back(unknownProb.at(i));
	}

	#if DBG_EM
	cout << "EM result end" << endl;
	#endif

}

void EM::getAlleleFreqs(double aFreq[2][3]){
	int tempVec[2];
	int imark;
	tempVec[0] = tempVec[1] = 0;
	aFreq[0][0] = aFreq[0][1] = aFreq[1][0] = aFreq[1][1] = 0;
	aFreq[1][2] = aFreq[0][2] = 0;
	for(imark = 0; imark < 2; imark++){
		int hapIndex, iallele;
		for(iallele = 1; iallele <= numAlleles.at(imark); iallele++){
			tempVec[imark] = iallele;
			hapIndex = (numAlleles.at(1) + 1)*tempVec[0] + tempVec[1];
			aFreq[imark][iallele] = unknownProb.at(hapIndex);
		}
		tempVec[imark] = 0;
	}
}

vector<int> EM::getNumberAlleles(){
	return numAlleles;
}

vector<double> EM::getHapProbs(){
	vector<double> t;
	for(int i=0;i<numHaps;i++){
		t.push_back(unknownProb.at(i));
	}
	return t;
}

/**
 * Return the haplotype frequencies for each possible haplotype.
 *
 */
vector<double> EM::getEMFreqs(){

	#if DBG_EM
	cout << "EM get freqs start" << endl;
	cout << "Sizes: " << endl;
	cout << "emFreq: " << emFrequencies.size() << endl;
	cout << "numHaps: " << numHaps << endl;
	cout << "incompleteTable: " << incompleteTable.size() << endl;
	cout << "setFor: " << setFor << endl;
	#endif

	vector<double> hapFreq;
	for(int i=0;i < numSplits;i++) hapFreq.push_back(0);

	haplotype *nextHaplotype;

	for(int ihap = 0; ihap < numHaps; ihap++){
		nextHaplotype = incompleteTable.at(ihap);
		if(nextHaplotype->getNumDescendents() == 0){
			vector<int> alleles;
			for(int i=0;i < setFor; i++){ // setFor == number of snps.
				#if DBG_EM
				cout << "EM inner loop start" << endl;
				#endif
				alleles.push_back(nextHaplotype->getAllele(i));
				#if DBG_EM
				cout << "EM inner loop end" << endl;
				#endif	

			}
			#if DBG_EM
			cout << "EM haptoindex: " << hapToIndex(alleles) << endl;
			cout << "alleles: " << alleles.size() << " - ";
			for(unsigned int i=0; i < alleles.size(); i++)
				cout << alleles.at(i) << " ";
			cout << endl;
			#endif
			int hapIdx = hapToIndex(alleles);
			if(hapIdx >= 0){
				hapFreq.at(hapIdx) = emFrequencies.at(ihap);
			}else{
				cout << "throwing" << endl;
				throw EMAlgorithmFailureException();
			}
		}
	}

	#if DBG_EM
	cout << "EM get freqs start" << endl;
	#endif

	return hapFreq;
}

/*
 * Return index of a haplotype assuming it uses a binary number system
 * with 1=>0 and 2=>1.
 */
int EM::hapToIndex(const vector<int> &haps){

	int ret = 0;
	for(unsigned int i=0;i<haps.size();i++){
		ret *= 2;
		ret += haps.at(i) - 1;
	}
	return ret;
}

/**
 * Calculate all haplotype probs for an individual.
 */
vector<EMPersonalProbsResults> EM::getPersonalProbabilities(){

	vector<EMPersonalProbsResults> ret;

	int *newhap, *partnerhap;
	newhap = new int[numSplits];
	partnerhap = new int[numSplits];
	vector<bool> tmp(numHaps, 0);
	vector<vector<bool> > used(numHaps, tmp);
	personalProb = vecops::getDblVec(numHaps, numHaps);

	for(int iindiv = 0; iindiv < static_cast<int>(inTotal.at(0).size()); iindiv++){

		for(int ihap = 0; ihap < numHaps; ihap++){
			for(int jhap = 0; jhap < numHaps; jhap++){
				personalProb[ihap][jhap] = 0;
				used[ihap][jhap] = false;
			}
		}

		double denom = 0;
		for(int isplit = 0; isplit < numSplits/2; isplit++){
			getRelHaps(iindiv, isplit, newhap[isplit], partnerhap[isplit]);
//			newhap[isplit] = getRelHap1(iindiv, isplit, numAlleles);
//			partnerhap[isplit] = getRelHap2(iindiv, isplit, numAlleles);
			denom += emFrequencies[newhap[isplit]]*emFrequencies[partnerhap[isplit]];
		}

		for(int isplit = 0; isplit < numSplits/2; isplit++){
			personalProb[newhap[isplit]][partnerhap[isplit]] += emFrequencies[newhap[isplit]]*emFrequencies[partnerhap[isplit]]/denom;
		}

		for(int isplit = 0; isplit < numSplits/2; isplit++){
			if(!used[newhap[isplit]][partnerhap[isplit]]){
				resolve(personalProb, newhap[isplit], partnerhap[isplit]);
			}
			used[newhap[isplit]][partnerhap[isplit]] = true;
		}

		int lefthap, righthap;
		for(lefthap = 0; lefthap < numHaps - 1; lefthap++){
			for(righthap = lefthap + 1; righthap < numHaps; righthap++){
				personalProb[lefthap][righthap] += personalProb[righthap][lefthap];
				personalProb[righthap][lefthap] = 0;
			}
		}

		for(lefthap = 0; lefthap < numHaps; lefthap++){
			for(righthap = lefthap; righthap < numHaps; righthap++){
				int numleftdesc, numrightdesc;
				haplotype *leftHaploPtr, *rightHaploPtr;
				leftHaploPtr = incompleteTable[lefthap];
				numleftdesc = leftHaploPtr->getNumDescendents();
				rightHaploPtr = incompleteTable[righthap];
				numrightdesc = rightHaploPtr->getNumDescendents();
				if(((numrightdesc == 0) && (numleftdesc == 0)) && (personalProb[lefthap][righthap] > 0)){
					EMPersonalProbsResults e;
					e.personId = iindiv;

					vector<int> t(setFor,0);
					for(int i=0; i < setFor; i++){
						t[i] = index2Allele(i, lefthap);
					}
					e.leftHap = hapToIndex(t);

					for(int i=0; i < setFor; i++){
						t[i] = index2Allele(i, righthap);
					}
					e.rightHap = hapToIndex(t);
					e.prob = personalProb[lefthap][righthap];
					ret.push_back(e);
				}
			}
		}
	}

	return ret;
}

/*
 * Stopgap hack to translate the undocumented internal EM storage pattern.
 */
int EM::index2Allele(int imark, int haplotypeIndex)
{
   int jmark;
   int product = 1;
   int alleleNumber;

   for (jmark = static_cast<int>(numAlleles.size()) - 1 ; jmark > imark ; jmark--){
      product *= (numAlleles[jmark] + 1);
   }
   alleleNumber = (haplotypeIndex/product) % (numAlleles[imark] + 1);
   return(alleleNumber);
}

/*
 * No comments in original code.  I have NO IDEA what this does.
 */
void EM::resolve(vector<vector<double> > &personalProb, int lefthap, int righthap)
{
	double epsilon = 0.000001;
	double prob = personalProb[lefthap][righthap];
	if(prob <= epsilon){
		return;
	}

	haplotype *leftHaploPtr;
	haplotype *rightHaploPtr;
	int idesc, numleftdesc, numrightdesc, nextdeschap;

	leftHaploPtr = incompleteTable[lefthap];
	numleftdesc = leftHaploPtr->getNumDescendents();
	rightHaploPtr = incompleteTable[righthap];
	numrightdesc = rightHaploPtr->getNumDescendents();
	if(numleftdesc > 0){
		for(idesc = 0; idesc < numleftdesc; idesc++){
			nextdeschap = leftHaploPtr->getDescendent(idesc);
			personalProb[nextdeschap][righthap] += emFrequencies[lefthap] ? personalProb[lefthap][righthap]*emFrequencies[nextdeschap]/	emFrequencies[lefthap] : 0 ;
			resolve(personalProb, nextdeschap, righthap);
		}
	}
	else if(numrightdesc > 0){
		for(idesc = 0; idesc < numrightdesc; idesc++){
			nextdeschap = rightHaploPtr->getDescendent(idesc);
			personalProb[lefthap][nextdeschap] += emFrequencies[righthap] ? personalProb[lefthap][righthap]*emFrequencies[nextdeschap] / emFrequencies[righthap] : 0 ;
			resolve(personalProb, lefthap, nextdeschap);
		}
	}
}






/*
 * Count unique elements in array.
 * ASSUMPTION: elements have values between 0 and 4.
 * We are looking for presence of a 2,3,or 4.  If we see one
 * return 2.  Otherwise, return 1.
 *
 */
int EM::uniqueElements(const vector<short> &v){

	bool has_a_one = false;
	for(unsigned int i=0; i < v.size(); ++i){
		if(v.at(i) >= 2){
			return 2;
		}else if(v.at(i) == 1){
			has_a_one = true;
		}
	}
	if(has_a_one)	return 1;
	return 0;
}
/****************************************************************************/
/**
 * Perform initilization for the EM by giving all haplotypes
 * equal probability
 */
void EM::initEMAlgorithm()
{
   int ihap;
   haplotype *nextHaplo;

   for(ihap = numHaps - 1; ihap >= 0; ihap--){
      nextHaplo = incompleteTable[ihap];
      if(nextHaplo->getNumDescendents() == 0){
         unknownProb[ihap] = 1/static_cast<double>(numHaps);
      }
      else{
         unknownProb[ihap] = 0;
      }
   }
}
void EM::computeIncompleteFreqs()
{
	haplotype *nextHaplo;
	int numdesc, idesc, ihap;
	int nextdeschap;

	for(ihap = numHaps - 1; ihap >= 0; ihap--){
		nextHaplo = incompleteTable[ihap];
		numdesc = nextHaplo->getNumDescendents();
		if(numdesc != 0){
			for(idesc = 0; idesc < numdesc; ++ idesc){
				nextdeschap = nextHaplo->getDescendent(idesc);
				unknownProb[ihap] += unknownProb[nextdeschap];
			}
		}
	}
}

void EM::computeCompleteFreqs()
{
	haplotype *nextHaplo;
	int ihap, numdesc;

	for(ihap = numHaps - 1; ihap >= 0; --ihap){
		nextHaplo = incompleteTable[ihap];
		numdesc = nextHaplo->getNumDescendents();
		if(numdesc == 0){
			unknownProb[ihap] = tempCount[ihap]/(2*static_cast<int>(inTotal.at(0).size()));
		}
		else{
			unknownProb[ihap] = 0;
		}
	}
}

void EM::cycleThroughPeople()
{
	int ihap;
	for(ihap = 0; ihap < numHaps; ihap++){
		tempCount[ihap] = 0.0;
	}

	double denom;
	int isplit;
	int *newhap = new int[numSplits];
	int *partnerhap = new int[numSplits];
	for(unsigned int iindiv = 0; iindiv < inTotal.at(0).size(); iindiv++){
		denom = 0;
		for(isplit = 0; isplit < numSplits/2; isplit++){
			getRelHaps(iindiv, isplit, newhap[isplit], partnerhap[isplit]);
//			newhap[isplit] = getRelHap1(iindiv, isplit, numAlleles);
//			partnerhap[isplit] = getRelHap2(iindiv, isplit, numAlleles);
			denom += unknownProb[newhap[isplit]]*unknownProb[partnerhap[isplit]];
		}

		for(isplit = 0; isplit < numSplits/2; ++isplit){
			tempCount[newhap[isplit]] += (unknownProb[newhap[isplit]]*unknownProb[partnerhap[isplit]])/denom;
			tempCount[partnerhap[isplit]] += (unknownProb[newhap[isplit]]*unknownProb[partnerhap[isplit]])/denom;
		}
	}

	delete [] newhap;
	delete [] partnerhap;
}

void EM::distributeCounts()
{
	int ihap, numdesc, idesc, nextdeschap;
	haplotype *nextHaplo;

	for(ihap = 0; ihap < numHaps; ++ihap){
		if(tempCount[ihap] != 0){
			nextHaplo = incompleteTable[ihap];
			numdesc = nextHaplo->getNumDescendents();
			for(idesc = 0; idesc < numdesc; idesc++){
				nextdeschap = nextHaplo->getDescendent(idesc);
				tempCount[nextdeschap] += tempCount[ihap]*unknownProb[nextdeschap]/unknownProb[ihap];
			}
		}
	}
}

void EM::reweightKnownHaps()
{
	haplotype *nextHaplo;
	int ihap;
	double checkSum = 0, epsilon = 0.000001;
	for(ihap = 0; ihap < numHaps; ++ihap){
		nextHaplo = incompleteTable[ihap];
		if(nextHaplo->getNumDescendents() == 0){
			if(unknownProb[ihap] < epsilon){
				unknownProb[ihap] = 0;
			}
			checkSum += unknownProb[ihap];
		}
	}
	for(ihap = 0; ihap < numHaps; ++ihap){
		unknownProb[ihap] /= checkSum;
	}
}

void EM::getRelHaps(int indiv, int split, int &newHapNum1, int &newHapNum2){
	newHapNum1 = newHapNum2 = 0;
	int mask, bit;
	int weight1, weight2;
	for (unsigned int imark=0; imark < numAlleles.size(); ++imark){
		// Grab weights
		weight1 = haplotype1.at(indiv).index2Allele(imark, numAlleles);
		weight2 = haplotype2.at(indiv).index2Allele(imark, numAlleles);
		
		// Calculate the correct bit.
		mask = 1<<(numAlleles.size() - imark - 1);//static_cast<int>(pow(2.0, (static_cast<int>(numAlleles.size()) - imark - 1)));
		bit = (split & mask);
		if (bit > 0){
			// equiv of old getRelHap1.
			newHapNum1 = (numAlleles[imark] + 1)*newHapNum1 + weight1;
			newHapNum2 = (numAlleles[imark] + 1)*newHapNum2 + weight2;
		}else{
			newHapNum1 = (numAlleles[imark] + 1)*newHapNum1 + weight2;
			newHapNum2 = (numAlleles[imark] + 1)*newHapNum2 + weight1;
		}
	}
}


int EM::getRelHap1(int indiv, int split, vector<int> &numAlleles) const
{
	int newHapNum = 0, mask;
	int bit, whichAllele;

	int imark;
	for(imark = 0; imark < static_cast<int>(numAlleles.size()); imark++){
		
		// The final product, bit, is equivalent to asking whether the single digit in
		// split is 1.
		
		mask = 1<<(numAlleles.size() - imark - 1);//static_cast<int>(pow(2.0, (static_cast<int>(numAlleles.size()) - imark - 1)));
		
		bit = (split & mask);
		whichAllele = bit > 0 ? haplotype1.at(indiv).index2Allele(imark, numAlleles) : haplotype2.at(indiv).index2Allele(imark, numAlleles);
		newHapNum = (numAlleles[imark] + 1)*newHapNum + (whichAllele);
	}
	return newHapNum;
}
/* Only difference from one above is haplotype1,2 are switched on whichAllele = line. */
int EM::getRelHap2(int indiv, int split, vector<int> &numAlleles) const
{
	int newHapNum = 0, mask;
	int bit, whichAllele;

	int imark;
	for(imark = 0; imark < static_cast<int>(numAlleles.size()); imark++){
		mask = 1 <<(numAlleles.size() - imark - 1);//static_cast<int>(pow(2.0, (static_cast<int>(numAlleles.size()) - imark - 1)));
		bit = (split & mask)/mask;
		whichAllele = bit ? haplotype2.at(indiv).index2Allele(imark, numAlleles) : haplotype1.at(indiv).index2Allele(imark, numAlleles);
		newHapNum = (numAlleles[imark] + 1)*newHapNum + (whichAllele);

	}
	return newHapNum;
}
