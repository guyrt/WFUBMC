#include "randwh.h"

using namespace std;

RandWH::RandWH(){
	seed[0] = rand();
	seed[1] = rand();
	seed[2] = rand();
	seed[3] = rand();
}

/*
* Set seeds to the four values in seeds[].
*/
void RandWH::init(long seeds[])
{
  int j;
  for(j=0; j < 4; j++)
    {
      seed[j] = seeds[j]; // We are filling in our held seeds.
    }
}

/*
* 	Takes a filename.  Puts first four numbers of the last line of the file in as the random seeds.
*/
void RandWH::init(char * filename)
{
	ifstream infile (filename, ios::ate);
	infile.unget();
	infile.unget();
	infile.unget();
	char c;
	while ((c = infile.get()) >= 32)
	{
	infile.unget();
	infile.unget();
	}
	
	
	infile >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	
}

/*
*	 Fill the input array with the 4 values of random seeds.
*/
void RandWH::state(long *a)
{
  for (int i=0; i < 4; i++)
    {
      a[i] = seed[i];
    }
}

/*
* Append the file with a list of the random seeds.
*/
void RandWH::state(char * filename)
{
		ofstream outfile (filename, ios::app);
		outfile << seed[0] << " " << seed[1] << " " << seed[2] << " " << seed[3] << endl;
		outfile.close();
}

/*
 * Update seeds according to W-H algorithm then return a random number.
 */
double RandWH::get()
{  
  long t1, t2;
  double r = 0.0;

  // Seed 0
  t1 = seed[0] % 185127;  // Be careful here: t1 and t2 must be longs
  t1 *= 11600;		  // If we let them be doubles then we don't have our generator.
  t2 = seed[0] / 185127;
  t2 *= 10379;
  t1 -= t2;
  if (t1 < 0) t1 += p1;

  seed[0] = t1;
  r += (static_cast<double>(t1)/p1);
  // cout << r << " ";

  // Seed 1
  t1 = seed[1] % 45688;
  t1 *= 47003;
  t2 = seed[1] / 45688;
  t2 *= 10479;
  t1 -= t2;
  if (t1 < 0) t1 += p2;

  seed[1] = t1;
	
  r += (static_cast<double>(t1)/p2);
	
  // cout << t1 << " " << t1/p2 << " ";
  // Seed 2
  t1 = seed[2] % 93368;
  t1 *= 23000;
  t2 = seed[2] / 93368;
  t2 *= 19423;
  t1 -= t2;
  if (t1 < 0) t1 += p3;

  seed[2] = t1;
	
  r += (static_cast<double>(t1)/p3);
	
  //  cout << r << " ";
	
  // Seed 3
  t1 = seed[3] % 65075;
  t1 *= 33000;
  t2 = seed[3] / 65075;
  t2 *= 8123;
  t1 -= t2;
  if (t1 < 0) t1 += p4;

  seed[3] = t1;
  r += (static_cast<double>(t1)/p4);
  //  cout << r << endl;
  return r - static_cast<long>(r);
	
	
}


