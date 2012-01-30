#include "vecops.hh"


/**
 * getVec.  Get a double vector of the specified size.
 * 
 * @param size1 Size in first component.
 * @param size2 Size in second component.
 * @return double vector of size (size1, size2)
 */
std::vector<std::vector<double> > vecops::getDblVec(int size1, int size2){
	
	std::vector<double> tmp(size2, 0);
	std::vector<std::vector<double> > ret(size1, tmp);
	return ret;
	
}
