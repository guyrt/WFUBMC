#include "stringutils.h"

/**
 * Correctly format left justified string with lead and trailing spaces.
 * 
 */
std::string strnutils::spaced_string(std::string s, int numspaces, int buffer){
	#if DB_V_STRNUTILS
	std::cout << "s " << s << std::endl;
	#endif
	std::stringstream ss;
	for(int i=0;i<buffer;i++) ss << " ";

	int left = numspaces-s.length();
	ss << s;
	for(int i=0;i<left;i++) ss << " " ;
	#if DB_V_STRNUTILS
	std::cout << "done" << std::endl;
	#endif
	return ss.str();
}

/**
 * Create a string with the 'number' embedded in 'numspaces' spaces.
 * 
 */
std::string strnutils::spaced_number(int number, int numspaces, int buffer){
	#if DB_V_STRNUTILS
	std::cout << "n1 " << number << std::endl;
	#endif
	std::stringstream ss;
	if(number < 0) buffer--;
	for(int i=0;i<buffer;i++) ss << " ";
	
	
	if(number > pow(10, numspaces)){
		// make if an inf.
		for(int i=0;i < numspaces; i++) ss << "9";
		return ss.str();
	}
	
	int sz = int_log(number)+1;
	for(int i=sz;i<numspaces;i++)
		ss << " ";
	ss << number;
	#if DB_V_STRNUTILS
	std::cout << "done" << std::endl;
	#endif
	return ss.str();
	
}

/**
 * Create a string with the 'number' embedded in 'numspaces' spaces.
 * 
 * Doubles can satisfy several criteria:
 * 	very closeto zero, requiring scientific notation.
 *  have digits on either side of the decimal.
 *  
 * 	@param number To be placed
 *  @param numspaces Number of slots alloted
 *  @param precision The desired numerical precision to be printed.
 *  @param buffer Number of spaces to print ahead of the number.
 *  @return std::string A string with the correct embedded number.
 * 
 * 	@bug Sometimes prints with too many spaces if a large exponent is involved.
 */
std::string strnutils::spaced_number(double number, int numspaces, int precision, int buffer){
	
	std::stringstream ss;
	for(int i=0;i<buffer;i++) ss << " ";
	
	double lim = pow(.1,precision);
	
	#if DB_V_STRNUTILS
	std::cout << "n2 " << number << " spaces " << numspaces << " prec " << precision << " lim " << lim <<  std::endl;
	#endif
	if(number < lim && number > 0){
		
		int exp = small_log(number);
		exp += 2; // for 'e-'
		
		int leading = numspaces - exp - 1;
		if(number < 0) leading -= 1;
		
		if(numspaces - exp - 2 > 0){
			if(numspaces - exp - 2 > precision){
				ss.precision(precision);
				leading -= precision;
			}else{
				leading = 0;
				ss.precision(numspaces-exp-2);
			}
		}else{
			// ran out of room.
			ss.precision(0);
		}
		
		
		for(int i=0;i<leading;i++) ss << " ";
		ss << std::scientific << number;
	}else if(number > pow(10, numspaces)){
		// make it an inf.
		for(int i=0;i < numspaces; i++) ss << "9";
		
	}else if(-1*number > pow(10, numspaces)){
		// make it an inf.
		ss << "-";
		for(int i=0;i < numspaces-1; i++) ss << "9";
		
	}else{
		ss << std::fixed;
		int leading = numspaces - precision - 1;
		leading -= (int_log(number) + 1);
		if(number < 0) leading -= 1;
		for(int i=0;i<leading;i++) ss << " ";
		if(leading < 0) precision += leading;
		ss.precision(precision);
		ss << number;
		
	}
	#if DB_V_STRNUTILS
	std::cout << "done" << std::endl;
	#endif
	return ss.str();
}


// Return the integer logarithm, so the log rounded down.  
// Possible to use a very cheap algorithm.
int strnutils::int_log(int i){
	if(i < 0) i *= -1;
	int r = 10;
	int c = 0;
	while(i >= r){
		r *= 10;
		c++;
	}
	return c;
}
int strnutils::int_log(double i){
	
	if(i < 0) i *= -1;
	double r = 10;
	int c = 0;
	while(i > r){
		r *= 10;
		c++;
	}
	return c;
}
/* Return number of digits in exponent. */
int strnutils::small_log(double i){
	
	if(i < 0) i *= -1;
	
	if(i >= 0.1){
		return -1;
	}else if(i >= 0.01){
		return 0;
	}else if(i >= 10e-10){
		return 1;
	}else if(i >= 10e-100){
		return 2;
	}
	return 3;
}

/**
 * Test method.
 * For now, record the results in cout.
 */
void strnutils::test(){
	
	std::cout << "Log: " << int_log(7.08507e+09) << std::endl;
	
	
	std::cout << "Printing 0.000" << std::endl;
	std::cout << "Log: " << int_log(0.000) << std::endl;
	std::cout << " --------" << std::endl << spaced_number(0.0, 8, 2, 1) << std::endl;
	std::cout << " --------" << std::endl << spaced_number(0.0, 8, 3, 1) << std::endl;
	std::cout << " --------" << std::endl << spaced_number(0.0, 8, 4, 1) << std::endl;
	std::cout << " --------" << std::endl << spaced_number(0.0, 8, 5, 1) << std::endl;

	std::cout << "Printing 23.000" << std::endl;
	std::cout << "Log: " << int_log(23.000) << std::endl;
	std::cout << " --------" << std::endl << spaced_number(23.0, 8, 2, 1) << std::endl;
	std::cout << " --------" << std::endl << spaced_number(23.0, 8, 3, 1) << std::endl;
	std::cout << " --------" << std::endl << spaced_number(23.0, 8, 4, 1) << std::endl;
	std::cout << " --------" << std::endl << spaced_number(23.0, 8, 5, 1) << std::endl;
	
	std::cout << "Printing 234.000" << std::endl;
	std::cout << "Log: " << int_log(234.000) << std::endl;
	std::cout << " --------" << std::endl << spaced_number(234.0, 8, 2, 1) << std::endl;
	std::cout << " --------" << std::endl << spaced_number(234.0, 8, 3, 1) << std::endl;
	std::cout << " --------" << std::endl << spaced_number(234.0, 8, 4, 1) << std::endl;
	std::cout << " --------" << std::endl << spaced_number(234.0, 8, 5, 1) << std::endl;
	
	std::cout << "Printing .000000000000234" << std::endl;
	std::cout << "Log: " << small_log(.000000000000234) << std::endl;
	std::cout << " ------" << std::endl << "|" << spaced_number(.000000000000234, 6, 2, 0) << "|" << std::endl;
	std::cout << " -------" << std::endl << "|" << spaced_number(.000000000000234, 7, 3, 0) << "|" << std::endl;
	std::cout << " --------" << std::endl << "|" << spaced_number(.000000000000234, 8, 4, 0) << "|" << std::endl;
	std::cout << " --------" << std::endl << "|" << spaced_number(.000000000000234, 8, 5, 0) << "|" << std::endl;
}
