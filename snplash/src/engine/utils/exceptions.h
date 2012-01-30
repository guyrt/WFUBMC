/*
 *      exceptions.h
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

/*
 * Exceptions are going to be a set of classes.
 */

#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>

using namespace std;

class ADTException {
	public:
	virtual ~ADTException(){};	
};

class DataException : public ADTException { 
	public :
		virtual ~DataException(){};
		string message;
};

/**
 * @class RunOrderException
 * 
 * Used when something called out of order.
 * 
 */
class RunOrderException : public ADTException {
	public:
		string message;
};

class StatsException : public ADTException { };

class InvalidChiSquareException : public StatsException {};

class GammaFxnFailureException : public StatsException {};

class FDistributionException : public StatsException {};

class EMAlgorithmNoSetup : public ADTException {};

class EMAlgorithmRunResultsMismatch : public ADTException {};

class EMAlgorithmFailureException : public ADTException {};

class LogisticRegressionException : public ADTException { };

class NewtonRaphsonFailureEx : public LogisticRegressionException { };

class NewtonRaphsonIterationEx : public LogisticRegressionException { };

class ConditionNumberEx : public LogisticRegressionException {
	public:
	double conditionNumber; // in 1-norm
	ConditionNumberEx(double c){
		conditionNumber = c;
	}
};

class InvalidNewRaphInputEx : public LogisticRegressionException { };

class MatrixNotSquareEx : public LogisticRegressionException { };

class DeterminantCalculationEx : public LogisticRegressionException { };

class SingularMatrixEx : public LogisticRegressionException {};

class QSnpgwaException : public ADTException { };

class LinearRegressionException : public ADTException {  };

#endif


