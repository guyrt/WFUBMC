/*
 *      chisquare.h
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

/**
 * @class Statistics
 * @author Unknown, Aggregated by Richard T. Guy
 * 
 * Static class contains statistical tools.
 * 
 * Current tool set:
 * 
 * 	Chi2Prob
 *  
 */

#ifndef STATS_H
#define STATS_H

#include "exceptions.h"
#include <math.h>
#include <iostream>

using namespace std;

class Statistics {
	
	public :
	
		/**
		 * Convert a chi^2 value with a given number of degrees
		 * of freedom to a p-value.  
		 * 
		 * Throws
		 * 		EX_INVALID_CHI_SQUARE
		 * 		EX_GAMMA_FXN_FAILURE
		 * 
		 */
		static double chi2prob(double chi2val, double df); // Calculate degrees of freedom.
		static double normalPValue(double value);
		static double tdist(double t, int df);
		static double fdist(double x, int mval, int nval);
		
		/// Helpers for chi2prob
		static double gammq(double, double);
		static void gser(double *, double, double, double *);
		static void gcf(double *, double, double, double *);
		static double gammln(double);
		static double gamma_log(double x);
		
		static double normal_01_cdf ( double x );
		
		static double ITMAX(){return 1000000.0;}
		static double EPS(){return 3.0e-7;}
		
		/// Helpers for Tdist
		static double beta_inc ( double a, double b, double x );
		static double beta(double x, double y);
		
		/// General utilities
		static double d_epsilon();
		
		/// Test code
		static void runAllTests();
		static void test_beta_inc();
};
#endif
