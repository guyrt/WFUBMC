/*
 *      adtree.h
 *      
 *      Copyright 2009 Richard T. Guy <guyrt7@wfu.edu>
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


#ifndef AD_TREE_H
#define AD_TREE_H

/**
 *
 *
 * 
 *
 */
#include "../combinable.h"
#include "ad_rule.h"
#include "ad_tree_data.h"
#include "../../param/engine_param_reader.h" 
#include <limits.h>

class ADTree : public Classifies { // Which means we're also an engine.

    public:
        explicit ADTree();
        explicit ADTree(DataAccess *);
      	~ADTree();
        
        /// Must declare these as part of Engine abstract class.
        virtual void init();
        virtual void preProcess();
        virtual void process();
        virtual void enslave(EngineParamReader *);
        virtual void test();
        
        void report(map<string, int> &, map<string, AD_Rule> &); // Return information about the tree.
        void report_leaves(map<string, int> &, map<string, AD_Rule> &); // Return information about the leaves.

        AD_Data get_tree(){return tree;}
        string print_tree(){ return tree.print(data); }
		
	void print_tree_to_file();
        
        /// Classification
        virtual void classify(int a[4]); // Run test on the actual data.
        virtual void classify(int a[4], SnpData *); // Run test on passed data.
        virtual int classify(vector<short> &); // Run test on single individual.

    protected:

	bool haveOwner;
	double fudge;
	
	vector<double> weight_vec;
        
        EngineParamReader *ad_param;
        int order_in_bag;
        
        void weights(double * , int, long, Condition::comparison, short); // Calc pos + neg weights for a given attribute and condition.
        void weights2(double *, int, long); // Calc pos and neg weights for each type of score.
        void weights(double *);
        
        double score_categorical(int , long , Condition::comparison &, short &);
        double score_ordinal(int , long, short &, vector<double> &weight_vec);
        
        double minimize(int *, long * , Condition::comparison * , short *);
        double minimize_ordinal(int &, long &, Condition::comparison &, short &);
        bool keep_going(); // Will keep track of stopping conditions on the tree.
	void delete_my_innards();
	
	/// Used in stopping criteria.
	vector<double> percentTrue;

	AD_Data tree; // Holds the tree built by this ad_tree engine.  Can be large, so stuck down here.

	void processNoCov(); // Process without covariates.


};

#endif
