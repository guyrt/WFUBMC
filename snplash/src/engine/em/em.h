/*
 *      em.h
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

#ifndef EXPMAX_H
#define EXPMAX_H

/**
 * Expectation maximization algorithm.
 */

#include <vector>
#include <stdlib.h>
#include "../utils/statistics.h"
#include "haplotype.h"
#include "../utils/vecops.hh"

#define DBG_EM 0



struct EMPersonalProbsResults {
	int personId;
	int leftHap;
	int rightHap;
	double prob;
};

class EM {
	
	friend class HaploStats;
	
	public : 
	
		EM();
		~EM();
		
		void setup(const vector<vector<short> > &);
		
		bool run();
		// get the results in a [2][2] matrix. 	Correspond to allele freqs. A,a,B,b
		void result(double**, vector<int> &, vector<double> &); 
		
		// Getters
		void getAlleleFreqs(double a[2][3]);
		vector<int> getNumberAlleles();
		vector<double> getHapProbs();
		vector<double> getEMFreqs();
		
		vector<EMPersonalProbsResults> getPersonalProbabilities();
		
	protected : 
		vector<vector<short> > inTotal; // Input vectors aggregated.
		vector<int> numAlleles;
		
		int numberAlleles;
		
		int numHaps, numSplits;
		vector<double> unknownProb, tempCount, emFrequencies;// Used in core of EM algo.
		vector<haplotype *> incompleteTable; // same
		/* next two are initialized with the two alleles that come in from in1/2. */
		/* We break the <short> storage scheme back into two vectors of length two per person. */
		vector<haplotype> haplotype1, haplotype2; // index is by individual.
		
		int uniqueElements(const vector<short> &);
		void initEMAlgorithm();
		void computeIncompleteFreqs();
		void computeCompleteFreqs();
		void cycleThroughPeople();
		void distributeCounts();
		void reweightKnownHaps();
		
		void getRelHaps(int indiv, int split, int &newHapNum1, int &newHapNum2);
		
		/// @Depricated. Use the combined form.
		int getRelHap1(int indiv, int split, vector<int> &numAlleles) const;
		int getRelHap2(int indiv, int split, vector<int> &numAlleles) const;
		/* No bleeping idea what this thing does -- RTG*/
		void resolve(vector<vector<double> > &, int lefthap, int righthap);
		
		int setFor;
		
		int hapToIndex(const vector<int> &);
		int index2Allele(int imark, int haplotypeIndex); // translate the getRelHap format.
		
		// Used in the get personal probabilities method.
		// Declared outside the stack for size reasons.
		vector<vector<double> > personalProb; 
};

#endif
