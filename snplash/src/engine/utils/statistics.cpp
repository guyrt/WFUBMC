#include "statistics.h"

/**
 * Compute f-statistic for a given input. 
 * 
 * @author Matt Steiger
 * 
 * @param x The test statistic
 * @param mval Degrees of freedom explained
 * @param nval Degrees of freedom unexplained.
 * @return double The probability under the F(mval,nval) of observing x or less extreme.
 */
double Statistics::fdist(double x, int mval, int nval)
{
   double pvalue;

   if(x <= 0.0){
      pvalue = 0.0;
   }
   else{
      double arg1, arg2, arg3;
      arg1 = 0.5 * static_cast<double>(nval);
      arg2 = 0.5 * static_cast<double>(mval);
      arg3 = static_cast<double>(nval)/(static_cast<double>(nval) + static_cast<double>(mval) * x);
	  pvalue = beta_inc(arg1, arg2, arg3);
   }
   
   return(pvalue);
}

/**
 * Compute p-value for the T distribution
 * 
 * @param t - the T statistic
 * @param df - degrees of freedom
 * @return pvalue - the pvalue.
 */
double Statistics::tdist(double t, int df)
{
   double a2, b2, c2;
   double cdf, degFree;
   degFree = static_cast<double>(df);

   a2 = 0.5 * degFree;
   b2 = 0.5;
   c2 = degFree/(degFree + t*t);

   if(t <= 0.0){
      cdf = 0.5*beta_inc(a2, b2, c2);
   }
   else{
      cdf = 1.0 - 0.5*beta_inc(a2, b2, c2);
   }

   double pvalue;
   if(t <= 0.0){
      pvalue = 2*cdf;
   }
   else{
      pvalue = 2*(1.0 - cdf);
   }
   return(pvalue);
}


/*
 * Calculate pval from normal distribution.
 */
double Statistics::normalPValue(double value)
{
   double normalVal, pVal;
   if(value > 0){
      double negValue;
      negValue = -1.0 * value;
      normalVal = normal_01_cdf(negValue);
   }
   else{
      normalVal = normal_01_cdf(value);
   }
   pVal = 2.0*normalVal;

   return(pVal);
}

/*
 * Compute normal CDF.
 * 
 * written by J. Grab.
 * 
 */
double Statistics::normal_01_cdf ( double x )
{

  double a2 = 0.399903438504;
  double a3 = 5.75885480458;
  double a4 = 29.8213557808;
  double a5 = 2.62433121679;
  double a6 = 48.6959930692;
  double a7 = 5.92885724438;
  double b0 = 0.398942280385;
  double b1 = 3.8052E-08;
  double b2 = 1.00000615302;
  double b3 = 3.98064794E-04;
  double b4 = 1.98615381364;
  double b5 = 0.151679116635;
  double b6 = 5.29330324926;
  double b7 = 4.8385912808;
  double b8 = 15.1508972451;
  double b9 = 0.742380924027;
  double b10 = 30.789933034;
  double b11 = 3.99019417011;
  double cdf;
  double q;
  double y;
//
//  |X| <= 1.28.
//
  if ( fabs ( x ) <= 1.28 )
  {
    double a1 = 0.398942280444;
    y = 0.5 * x * x;

    q = 0.5 - fabs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5
      + a6 / ( y + a7 ) ) ) );
//
//  1.28 < |X| <= 12.7
//
  }
  else if ( fabs ( x ) <= 12.7 )
  {
    y = 0.5 * x * x;

    q = exp ( - y ) * b0 / ( fabs ( x ) - b1
      + b2 / ( fabs ( x ) + b3
      + b4 / ( fabs ( x ) - b5
      + b6 / ( fabs ( x ) + b7
      - b8 / ( fabs ( x ) + b9
      + b10 / ( fabs ( x ) + b11 ) ) ) ) ) );
//
//  12.7 < |X|
//
  }
  else
  {
    q = 0.0;
  }
//
//  Take account of negative X.
//
  if ( x < 0.0 )
  {
    cdf = q;
  }
  else
//***********************************************************************************
// This code can handle positive values of x, but it is best to avoid the           *
// subtraction in the code below that occurs when x is positive, due to concerns    *
// related to floating point arithmetic (subtracting near-zero or near-one values   *
// from 1).  One way to handle this is to pre-multiply your a positive input value  *
// by -1.0, and then account for the sign change OUTSIDE this function.             *
//***********************************************************************************
  {
    cdf = 1.0 - q;
  }

  return cdf;
}

/**
 * Convert a chi^2 value with a given number of degrees of freedom to a p-value.  
 * 
 * @author Matt Steigert
 * @param x2val The chi^2 distributed statistic
 * @param df Degrees of freedom
 * @return Value of chi2(x2val, df)
 * 
 * @throw InvalidChiSquareException() If df < 1 or x2val < 0.0
 */
double Statistics::chi2prob(double x2val, double df)
{
   if(df < 1){
      throw InvalidChiSquareException();
      cerr << "Chi Square Value = " << x2val;
   }
   if(x2val < 0.0){
      throw InvalidChiSquareException();
      cerr << "Degrees of Freedom = " << df;
   }
   
   double pval;
   pval = gammq(df/2.0, x2val/2.0);

   return pval;
}

double Statistics::gammq(double a, double x)
{
   double gamser, gammcf, gln, pval;
   if(x < (a + 1.0)){
      gser(&gamser, a, x, &gln);
      pval = 1.0 - gamser;
   }
   else{
      gcf(&gammcf, a, x, &gln);
      pval = gammcf;
   }
   return pval;
}
void Statistics::gser(double *gamser, double a, double x, double *gln)
{
   int n; 
   double sum, del, ap;

   *gln = gammln(a);
   if(x <= 0.0){
      *gamser=0.0;
      return;
   }
   else
   {
      ap = a;
      del = sum = 1.0/a;
      for(n = 1; n <= Statistics::ITMAX(); n++){
         ap += 1.0;
         del *= x/ap;
         sum += del;
         if(fabs(del) < fabs(sum)*Statistics::EPS()){
            *gamser=sum*exp(-x+a*log(x)-(*gln));
            return;
         }
      }
      throw GammaFxnFailureException();
      return;
   }
}
void Statistics::gcf(double *gammcf, double a, double x, double *gln)
{
   int n;
   double gold = 0.0, g, fac = 1.0, b1 = 1.0;
   double b0= 0.0, anf, ana, an, a1, a0 = 1.0;
   
   *gln = gammln(a);
   a1 = x;
   for(n = 1 ; n <= Statistics::ITMAX() ; n++){
      an = static_cast<double>(n);
      ana = an - a;
      a0 = (a1 + a0 * ana) * fac;
      b0 = (b1 + b0 * ana) * fac;
      anf = an * fac;
      a1 = x * a0 + anf * a1;
      b1 = x * b0 + anf * b1;
      if(a1 > pow(Statistics::EPS(),3)){ // CHANGED BY RTG to avoid f.p. compare.
         fac = 1.0/a1;
         g= b1 * fac;
         if(fabs((g - gold)/g) < Statistics::EPS()){
            *gammcf = exp(-x + a * log(x) - (*gln)) * g;
            return;
         }
         gold = g;
      }
   }
   throw GammaFxnFailureException();
}

/**
 * Compute the gamma function.
 * 
 * @author Matt Steigart ?
 * 
 * @param xx 
 * @return double 
 */
double Statistics::gammln(double xx)
{
   double x, y, tmp, ser;
   static double cof[6]={76.18009172947146, -86.50532032941677,
                         24.01409824083091, -1.231739572450155,
                         0.1208650973866179e-2, -0.5395239384953e-5};
   int j;  
   x = xx;
   y = xx;
  
   tmp = x + 5.5;
   tmp -= (x + 0.5) * log(tmp);
   ser = 1.000000000190015;
   for(j = 0 ; j <= 5 ; j++){
      ser += cof[j]/++y;
   }
   return (-tmp + log(2.5066282746310005 * ser/x));
}

/**
 * Compute natural log of the gamma function.
 * 
 * @param x
 * @return double
 */
double Statistics::gamma_log(double x)
{
  double c[7] = {-1.910444077728E-03, 8.4171387781295E-04, -5.952379913043012E-04, 
                 7.93650793500350248E-04, -2.777777777777681622553E-03, 
                 8.333333333333333331554247E-02, 5.7083835261E-03 };
  double corr;
  double d1 = - 5.772156649015328605195174E-01;
  double d2 =   4.227843350984671393993777E-01;
  double d4 =   1.791759469228055000094023;
  double frtbig = 1.42E+09;
  int i;
  double p1[8] = {
    4.945235359296727046734888, 
    2.018112620856775083915565E+02, 
    2.290838373831346393026739E+03, 
    1.131967205903380828685045E+04, 
    2.855724635671635335736389E+04, 
    3.848496228443793359990269E+04, 
    2.637748787624195437963534E+04, 
    7.225813979700288197698961E+03 };
  double p2[8] = {
    4.974607845568932035012064, 
    5.424138599891070494101986E+02, 
    1.550693864978364947665077E+04, 
    1.847932904445632425417223E+05, 
    1.088204769468828767498470E+06, 
    3.338152967987029735917223E+06, 
    5.106661678927352456275255E+06, 
    3.074109054850539556250927E+06 };
  double p4[8] = {
    1.474502166059939948905062E+04, 
    2.426813369486704502836312E+06, 
    1.214755574045093227939592E+08, 
    2.663432449630976949898078E+09, 
    2.940378956634553899906876E+010,
    1.702665737765398868392998E+011,
    4.926125793377430887588120E+011, 
    5.606251856223951465078242E+011 };
  double pnt68 = 0.6796875;
  double q1[8] = {
    6.748212550303777196073036E+01, 
    1.113332393857199323513008E+03, 
    7.738757056935398733233834E+03, 
    2.763987074403340708898585E+04, 
    5.499310206226157329794414E+04, 
    6.161122180066002127833352E+04, 
    3.635127591501940507276287E+04, 
    8.785536302431013170870835E+03 };
  double q2[8] = {
    1.830328399370592604055942E+02, 
    7.765049321445005871323047E+03, 
    1.331903827966074194402448E+05, 
    1.136705821321969608938755E+06, 
    5.267964117437946917577538E+06, 
    1.346701454311101692290052E+07, 
    1.782736530353274213975932E+07, 
    9.533095591844353613395747E+06 };
  double q4[8] = {
    2.690530175870899333379843E+03, 
    6.393885654300092398984238E+05, 
    4.135599930241388052042842E+07, 
    1.120872109616147941376570E+09, 
    1.488613728678813811542398E+010, 
    1.016803586272438228077304E+011, 
    3.417476345507377132798597E+011, 
    4.463158187419713286462081E+011 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double xbig = 4.08E+36;
  double xden;
  double xm1;
  double xm2;
  double xm4;
  double xnum;
  double xsq;
//
//  Return immediately if the argument is out of range.
//
   if(x <= 0.0 || xbig < x){
      return HUGE_VAL;
   }

   if(x <= d_epsilon()){
      res = -log(x);
   }
   else if(x <= 1.5){
      if(x < pnt68){
         corr = -log(x);
         xm1 = x;
      }
      else{
         corr = 0.0;
         xm1 = (x - 0.5) - 0.5;
      }

      if(x <= 0.5 || pnt68 <= x){
         xden = 1.0;
         xnum = 0.0;

         for(i = 0; i < 8; i++){
            xnum = xnum*xm1 + p1[i];
            xden = xden*xm1 + q1[i];
         }

         res = corr + (xm1*(d1 + xm1*(xnum/xden)));
      }
      else{
         xm2 = ( x - 0.5 ) - 0.5;
         xden = 1.0;
         xnum = 0.0;
         for ( i = 0; i < 8; i++ ){
            xnum = xnum*xm2 + p2[i];
            xden = xden*xm2 + q2[i];
         }
         res = corr + xm2*(d2 + xm2*(xnum/xden));
      }
   }
   else if ( x <= 4.0 )
   {
      xm2 = x - 2.0;
      xden = 1.0;
      xnum = 0.0;
      for(i = 0; i < 8; i++){
         xnum = xnum * xm2 + p2[i];
         xden = xden * xm2 + q2[i];
      }
      res = xm2*(d2 + xm2*( xnum/xden ) );
   }
   else if ( x <= 12.0 )
   {
      xm4 = x - 4.0;
      xden = - 1.0;
      xnum = 0.0;
      for(i = 0; i < 8; i++){
         xnum = xnum * xm4 + p4[i];
         xden = xden * xm4 + q4[i];
      }
      res = d4 + xm4*(xnum/xden);
   }
   else{
      res = 0.0;

      if(x <= frtbig){
         res = c[6];
         xsq = x*x;

         for( i = 0; i < 6; i++ ){
            res = res/xsq + c[i];
         }

      }

      res = res/x;
      corr = log(x);
      res = res + sqrtpi - 0.5*corr;
      res = res + x*(corr - 1.0);

   }

  return res;
}

/**
 * Return value of incomplete beta function.
 * 
 * @author Matt Steigart?
 * 
 * @param a
 * @param b
 * @param x
 * @return double 
 */
double Statistics::beta_inc ( double a, double b, double x )
{
  double cx;
  int i;
  int it;
  int it_max = 100000; // changed from 10,000 by RTG on 10/26/2010
  					   // this was somewhat of a stopgap, but appears to help.
  bool indx;
  int ns;
  double pp;
  double psq;
  double qq;
  double rx;
  double temp;
  double term;
  double tol = 1.0E-07;
  double value;
  double xx;

   if(a <= 0.0){
      throw FDistributionException();
   }

   if(b <= 0.0){
      throw FDistributionException();
   }

   if(x <= 0.0){
      return 0.0;
   }
   else if(1.0 <= x){
      return 1.0;
   }
//
//  Change tail if necessary and determine S.
//
   psq = a + b;
   if(a < (a + b) * x){
      xx = 1.0 - x;
      cx = x;
      pp = b;
      qq = a;
      indx = true;
   }
   else{
      xx = x;
      cx = 1.0 - x;
      pp = a;
      qq = b;
      indx = false;
   }

   term = 1.0;
   i = 1;
   value = 1.0;

   ns = static_cast<int>(qq + cx*(a + b));
//
//  Use Soper's reduction formulas.
//
   rx = xx/cx;

   temp = qq - static_cast<double>(i);
   if(ns == 0){
      rx = xx;
   }

   it = 0;

   while(1){
      it++;

      if(it_max < it){
         throw FDistributionException();
      }

      term = term*temp*rx/(pp + static_cast<double>(i));
      value = value + term;
      temp = fabs(term);

      if((temp <= tol) && (temp <= tol*value)){
         break;
      }

      i = i + 1;
      ns = ns - 1;

      if ( 0 <= ns ){
         temp = qq - static_cast<double>(i);
         if(ns == 0){
            rx = xx;
         }
      }
      else{
         temp = psq;
         psq = psq + 1.0;
      }
   }
//
//  Finish calculation.
//
   value = value*exp(pp*log(xx) + (qq - 1.0)*log(cx))/(beta(a, b)*pp);

   if(indx){
      value = 1.0 - value;
   }
   return value;
}

/**
 * Return value of beta function.
 * 
 * @author Matt Steigart ?
 * 
 * @param x
 * @param y
 * @return double
 */
double Statistics::beta(double x, double y)
{
   if(x <= 0.0 || y <= 0.0){
      throw FDistributionException();
   }

   return(exp(gamma_log(x) + gamma_log(y) - gamma_log(x + y)));
}

/**
 * Return roundoff error: the smallest number of form 1 / 2^k 
 * different from one.
 * 
 * @return r - Roundoff size.
 */
double Statistics::d_epsilon()
{
   double r;
   r = 1.0;

   while(1.0 < (double)(1.0 + r)){
      r = r/2.0;
   }

   return(2.0*r);
}

/**
 * Run full gamut of tests against stats functions.
 * 
 * 
 */
void Statistics::runAllTests(){
	test_beta_inc();	
}
/**
 * Test beta_inc() function by running full range of all parameters and outputing results
 * 
 * beta_inc takes three inputs:
 * a positive
 * b positive
 * x in range [=(0,1)
 * 
 * This was failing once:
 * 1974 0.5 0.999655
 */
void Statistics::test_beta_inc(){
	
	double a = 0.0,b = 0.0,x = 0.0; // define test values.
	double ma, mb, mx; // define max value
	ma = mb = 1;
	mx = 1.0;
	double na, nb, nx; // define granularity.
	na = nb = 100;
	nx = 100;
	for (int i=1; i < na; i++){
		a = ma / na * double(i);
		for (int j = 1; j < nb; j++){
		b = mb / nb * j;
			for(int k=1; k < nx; k++){
				x = mx / nx * k;
				try{
					beta_inc(a,b,x);
//					cout << a << " " << b << " " << x << " " << beta_inc(a,b,x) << endl;
				}catch(...){
					cout << a << " " << b << " " << x << " " << -111.111 << endl;
				}
			}
		}
	}
}
