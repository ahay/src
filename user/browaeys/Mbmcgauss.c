/* Correlated Gaussian joint probability distribution histogram generated with modified Box Mulller algorithm */
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

#include <math.h>
#include <rsf.h>

#define IA 16807
#define IM 2147483647 
#define AM (1.0/IM)
#define NTAB 32
#define IQ 127773
#define IR 2836
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

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

void gauss_joint(int *iseed, float m1, float m2, float s1, float s2, float r, float *y1, float *y2)
/* Returns two Gaussian positive distributed deviates y1 and y2 with correlation r */
/* using a modified Box Mulller algorithm and rand1_gen(idum) */
/* Expectations = m1,m2 */ 
/* Standard deviation = s1,s2 */
{
 float g1,g2,fac,rsq,v1,v2;
 do
 {
  v1 = 2.0*rand1_gen(iseed)-1.0;
  v2 = 2.0*rand1_gen(iseed)-1.0;
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

void histgauss(float **hist, int n, float r, float m1, float m2, float s1, float s2, int nbin, float dbin, float obin, int iseed) 
{
/* Joint Gaussian distribution histogram */

    int i,i1,i2,isd;
    float x1,x2,norm;   
     
    isd = iseed;
    norm = 1.0/(dbin*dbin*n);

    for (i = 0; i < n; i++) {

	gauss_joint(&isd,m1,m2,s1,s2,r,&x1,&x2);

	i1 = floorf((x1-obin)/dbin);
	i2 = floorf((x2-obin)/dbin);

	if ( (i1 >= 0) && (i1 < nbin) && (i2 >= 0) && (i2 < nbin) ) {
	    hist[i1][i2] += norm;
	}

    }

    return;
}

int main(int argc, char* argv[])
{
    int i,j;
    int n,iseed;
    int nbin;

    float dbin,obin;
    float m1,m2,s1,s2,r;
    float **hist;

    sf_axis axis_1,axis_2;
    sf_file out;

    sf_init (argc,argv);

    out = sf_output("out");
    sf_settype(out,SF_FLOAT);

    if (!sf_getint("n",&n)) n=100;
    /* number of random deviates pairs */

    if (!sf_getfloat("m1",&m1)) m1=0.0;
    /* mean for deviate 1 */

    if (!sf_getfloat("m2",&m2)) m2=0.0;
    /* mean for deviate 2 */

    if (!sf_getfloat("s1",&s1)) s1=1.0;
    /* standard deviation for deviate 1 */

    if (!sf_getfloat("s2",&s2)) s2=1.0;
    /* standard deviation for deviate 2 */

    if (!sf_getfloat("r",&r)) r=0.0;
    /* correlation coefficient */

    if (!sf_getint("nbin",&nbin)) nbin=51;
    /* number of bins for histogram */

    if (!sf_getfloat("dbin",&dbin)) dbin=0.1;
    /* histogram bin size */

    if (!sf_getfloat("obin",&obin)) obin=0.0;
    /* histogram origin */

    if (!sf_getint("iseed",&iseed)) iseed=-33;
    /* random generator seed */
    if (iseed >= 0) sf_error("Need strictly negative iseed");

    /* output file parameters */
    axis_1 = sf_maxa(nbin,obin,dbin);
    axis_2 = sf_maxa(nbin,obin,dbin);

    sf_oaxa(out,axis_1,1);
    sf_oaxa(out,axis_2,2);

    /* Histogram */
    hist = sf_floatalloc2(nbin,nbin);

    for (j = 0; j < nbin; j++) {
	for (i = 0; i < nbin; i++) {
	    hist[j][i] = 0.0;
	}
    }
    histgauss(hist,n,r,m1,m2,s1,s2,nbin,dbin,obin,iseed);

    /* output */
    sf_floatwrite(hist[0],nbin*nbin,out);

    exit(0);
}
