#ifndef RANDW_H
#define RANDW_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
/*
 * Implementation of the Wichmann-Hill random generator to 2005 specs.
 *
 * See Wichmann, Hill, "Generating Good psuedo-random numbers." from 2005.
 *
 * Author: Richard T. Guy
 *
 * September, 2008
 *
 * Usage:
 *   call initilize to set up the seeds.
 *   call get to get a random number.
 */
using namespace std;

class RandWH {

	public:
		RandWH();  // Standard constructor

		void init(long seeds[]);
		void init(char *);
		void state(char *);
		double get();
		void state(long *s);

	private:
		long seed[4];
		static const long p1 = 2147483579;
		static const long p2 = 2147483543;
		static const long p3 = 2147483423;
		static const long p4 = 2147483123;


};

#endif
