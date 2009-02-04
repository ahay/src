/* Monte Carlo integration of cos(2t).P(x1,x2).P(y1,y2) */
/* P is correlated Gaussian joint probability distribution */
/* t is the angle between vector fields (x1,y1) and (x2,y2) */
/*
  Copyright (C) 2009 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>
#include <math.h>

#define IA 16807
#define IM 2147483647 
#define AM (1.0/IM)
#define NTAB 32
#define IQ 127773
#define IR 2836
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI 3.141592653589793

float rand1_gen(int *idum)
{
/* Minimal random number generator of Park & Miller */
/* with Bays-Durham shuffle and added safeguards */
/* Returns a uniform random deviates between in [0.0,1.0[ (exclusive) */
/* Initial call with idum to any negative integer value */
/* idum must not be altered in successive calls in a sequence */

    int j;
    int k;
    static int iy=0;
    static int iv[NTAB];
    float temp;

    if (*idum <= 0 || !iy) {

	if (-(*idum) < 1) {
	    *idum = 1;
	} else {
	    *idum = -(*idum);
	}

	for (j = NTAB+7; j>=0; j--) {
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

    if ((temp = AM*iy) > RNMX) {
	return RNMX;
    } else {
	return temp;
    }

}

static void gauss_joint (float *x, float r, float m1, float s1, float m2, float s2, int *iseed) 
{
/* Generates one pair of correlated Gaussian distributed random deviates */
/* Computed with modified Box Mulller algorithm */

    float rsq,v1,v2,fac,gd1,gd2;

    rsq = 0.0;

    while ((rsq >= 1.0) || (rsq == 0.0)) {
	v1 = 2.*rand1_gen(&iseed) - 1.0;
	v2 = 2.*rand1_gen(&iseed) - 1.0;
	rsq = v1*v1 + v2*v2;
    }

    fac = sqrt(-2.*log(rsq)/rsq);
    gd1 = v1*fac;
    gd2 = v2*fac*sqrtf(1.-r*r);

    x[1] = gd1*s1 + m1;
    x[2] = (gd2 + r*gd1)*s2 + m2;
    
    return;
}


int main(int argc, char* argv[])
{
    int i,k;

    int nr;
    float dr,or,r;         /* correlation coefficient */

    int n;                 /* number of pairs of random deviates */
    int iseed;             /* random generator seed */
    float mm,ss;           /* mean and standard deviation */

    float x1,y1,x2,y2;
    float inrm1,inrm2,cs12;

    float **x,**y;         /* Gaussian distributed correlated random deviates */
    float *mc;

    sf_file in, out;

    sf_init (argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    /* correlation coefficient grid */
    if (!sf_histint(in,"n1",&nr)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dr)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&or)) sf_error("No o1= in input");

    if (!sf_getint("n",&n)) n=100;
    /* number of random deviates pairs */

    if (!sf_getfloat("mm",&mm)) mm=0.0;
    /* mean for deviates */

    if (!sf_getfloat("ss",&ss)) ss=1.0;
    /* standard deviation for deviates */

   if (!sf_getint("iseed",&iseed)) iseed=-33;
    /* random generator seed */
   if (iseed >= 0) sf_error("Need strictly negative iseed");

   /* memory allocations */
    mc = sf_floatalloc(nr);

    /* read input */
    sf_floatread(mc,nr,in);

    /* Correlated Gaussian deviates */
    x = sf_floatalloc2(2,n);
    y = sf_floatalloc2(2,n);

    for (k=0; k<nr; k++) {

	r = or + k*dr;

        /* Generates n pairs of Gaussian distributed correlated deviates */
	for (i=0; i<n; i++) {
	    gauss_joint(x[i],r,mm,ss,mm,ss,&iseed);
	}
	for (i=0; i<n; i++) {
	    gauss_joint(y[i],r,mm,ss,mm,ss,&iseed);
	}

        /* estimates <cos(2t)> = <2*cos^2(t)-1> */
	mc[k] = 0.0;

	for (i=0; i<n; i++) {

	    x1 = x[i][1];
	    y1 = y[i][1];
	    x2 = x[i][2];
	    y2 = y[i][2];

	    inrm1 = 1.0/sqrtf(x1*x1 + y1*y1);
	    inrm2 = 1.0/sqrtf(x2*x2 + y2*y2);

	    cs12 = (x1*x2 + y1*y2)*inrm1*inrm2;
	    mc[k] += 2.*cs12*cs12 - 1.0;

	}

	mc[k] /= n;
    }

    /* output */
    sf_floatwrite(mc,nr,out);

    exit(0);
}
