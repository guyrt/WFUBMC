#include "float_ops.hh"

/**
 * Compare two floating point numbers.
 * 
 * @param double Number 1
 * @param double Number 2
 * @return bool True if within EPS_EQUALS of each other.
 */
bool equal(double i, double d){


	if((i - d) > EPS_EQUALS){
		return false;
	}else if((i - d) < -EPS_EQUALS){
		return false;
	}
	return true;
}
