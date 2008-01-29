/**** FILE: generatorlib.c                                                           ****/
/**** Distributed deviates generators                                                ****/
/****                                                                                ****/
/**** ran0(idum)                         : random deviate generator                  ****/
/**** ran1(idum)                         : random deviate generator with shuffling   ****/
/**** expdev(idum,lambda)                : positive exponential distribution deviate ****/
/**** gasdev(idum,m,sigma)               : gaussian positive distribution deviate    ****/
/**** g2cdev(idum,m1,m2,s1,s2,r,*y1,*y2) : correlated gaussian distributed deviates  ****/
/**** lorz(idum,m,b)                     : Lorentzian positive distribution deviate  ****/
/**** logndev(idum,m,sigma)              : log normal positive distribution deviate  ****/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "generatorlib.h"

#define IA 16807
#define IM 2147483647 
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI 3.141592653589793

float ran0(long *idum)
/* Minimal random number generator of Park & Miller */
/* Returns a uniform random deviates between 0.0 and 1.0 */
/* Set or reset idum to any integer value except the */
/* value MASK to initialize the sequence */
/* idum must not be altered in successive calls in a sequence */
{
 long k;
 float ans;
 *idum ^= MASK;
 k = (*idum)/IQ;
 *idum = IA*(*idum-k*IQ)-IR*k;
 if (*idum < 0) *idum += IM;
 ans = AM*(*idum);
 *idum ^= MASK;
 return ans;
}

float ran1(long *idum)
/* Minimal random number generator of Park & Miller */
/* with Bays-Durham shuffle and added safeguards */
/* Returns a uniform random deviates between in [0.0,1.0[ (exclusive) */
/* Initial call with idum to any negative integer value */
/* idum must not be altered in successive calls in a sequence */
{
 int j;
 long k;
 static long iy=0;
 static long iv[NTAB];
 float temp;
 if (*idum <= 0 || !iy)
 {
  if (-(*idum) < 1) *idum = 1;
  else *idum = -(*idum);
  for (j = NTAB+7; j>=0; j--)
  {
   k = (*idum)/IQ;
   *idum = IA*(*idum-k*IQ)-IR*k;
   if (*idum < 0) *idum += IM;
   if (j < NTAB) iv[j] = *idum;
  }
  iy = iv[0];
 }
 k = (*idum)/IQ;
 *idum = IA*(*idum-k*IQ)-IR*k;
 if (*idum < 0) *idum += IM;
 j = iy/NDIV;
 iy = iv[j];
 iv[j] = *idum;
 if ((temp = AM*iy) > RNMX) return RNMX;
 else return temp;
}

float expdev(long *idum, float lambda)
/* Returns a strictly positive exponential distribution deviate */
/* with parameter lambda using ran1(idum) */
/* Expectation = lambda */
/* Variance = lambba^2 */
/* Standard deviation = lambda */
{
 float dum;
 do
  dum = ran1(idum);
 while (dum == 0.0);
 return -log(dum)*lambda; 
}

float gasdev(long *idum, float m, float sigma)
/* Returns a normal (gaussian) positive distribution deviate */
/* with mean m and variance sigma*sigma using ran1(idum) */
/* Box-Muller method */
/* Expectation = m */
/* Variance = sigma^2 */
/* Standard deviation = sigma */
{
 static int iset=0;
 static float gset;
 float fac,rsq,v1,v2;
 if (*idum < 0) iset =0;
 if (iset == 0)
 {
  do
  {
   v1 = 2.0*ran1(idum)-1.0;
   v2 = 2.0*ran1(idum)-1.0;
   rsq = v1*v1 + v2*v2;
  }
  while (rsq >= 1.0 || rsq == 0.0);
  fac = sqrt(-2.0*log(rsq)/rsq);
  gset = sigma*v1*fac + m;
  iset = 1;
  return (sigma*v2*fac + m);
 }
 else
 {
  iset = 0;
  return gset;
 }
}

void g2cdev(long *idum, float m1, float m2, float s1, float s2, float r, float *y1, float *y2)
/* Returns two normal (gaussian) positive distributed */
/* deviates y1 and y2 with correlation r, respective expectations */
/* m1 and m2 and respective variances s1*s1 and s2*s2 using ran1(idum) */
/* Box-Muller method */
/* Expectatiosn = m1, m2 */ 
/* Variances = s1^2,s2^2 */
/* Standard deviation = s1,s2 */
{
 float g1,g2,fac,rsq,v1,v2;
 do
 {
  v1 = 2.0*ran1(idum)-1.0;
  v2 = 2.0*ran1(idum)-1.0;
  rsq = v1*v1 + v2*v2;
 }
 while (rsq >= 1.0 || rsq == 0.0);
 fac = sqrt(-2.0*log(rsq)/rsq);
 g1 = v1*fac;
 g2 = v2*fac*sqrt(1.0-r*r);
 *y1 = s1*g1 + m1;
 *y2 = (g2+r*g1)*s2 + m2;
 return;
}

float lorz(long *idum, float m, float b)
/* Returns a strictly positive Lorentzian (Cauchy) distribution deviate */
/* with median m and half width b (stricly positive) using ran1(idum) */
/* p(x) = b/pi/((x-m)^2 + b^2)
/* All moments diverge
/* Expectation can be considered by symmetry equal to the median m */
/* Expectation = infinity */
/* Variance = infinity */
/* Standard deviation = infinity */
{
 float dum;
 do
  dum = ran1(idum);
 while (dum == 0.0);
 return (b*tan(PI*(dum-0.5)) + m);
}

float logndev(long *idum, float m, float sigma)
/* Returns a strictly positive log normal distribution deviate */
/* Gaussian : mean m, standard deviation sigma */
/* Uses gasdev and ran1(idum) */
/* p(x) = 1/(x*sigma*sqrt(2*pi))*exp(-(log(x)-m)^2/(2*sigma^2)) */
/* Expectation = exp(m + 0.5*sigma^2) */
/* Variance = expectation*expectation*(exp(sigma^2)-1)
/* Standard deviation = expectation*sqrt(exp(sigma^2)-1) */
{
 float dum;
 dum = gasdev(idum,m,sigma);
 return exp(dum);
}
