//      linear_regression.hh
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

/**
 * Linear regression algorithm.
 * 
 * Performs linear regression using the Ordinary Least Squares approach:
 * 
 * 		beta = inv(X'X) * X' * y
 * 
 */

#ifndef LINEAR_REGRESSION_H
#define LINEAR_REGRESSION_H

#include "../engine.h"


using namespace std;

struct LinRegStats{
		int i;
};

class LinearRegression{

		public:
			LinearRegression();
		
			vector<double> leastSquares(const vector<vector<double> > &data, const vector<double> &response);
			/* Calculate sum squared errors 
			 * and sum squared regression
			 * for a given set of betas */
			void sumSquaredStats(const vector<vector<double> > &data, const vector<double> &betas, 
						const vector<double> &response, double mean, double &sse, double &ssr);
			
			/*
			 * Calc residuals for given data and betas.
			 */
			vector<double> residuals(const vector<vector<double> > &data, const vector<double> &betas, 
							const vector<double> &response);
			
			/*
			 * Perform an anova on response variable.
			 */
			double anova(const vector<double> &response, const vector<double> &division, double means[4]);
			
			void test();
		
		private:
			static const double LINEAR_REGRESSION_CONDITION_NUMBER_LIMIT = 10e6;
	
};

#endif
