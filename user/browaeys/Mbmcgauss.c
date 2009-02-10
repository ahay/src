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
#include "bmgauss.h"

int main(int argc, char* argv[])
{
    int i,k;
    int nbin,nbh;          /* number of bins for histogram */ 
    int n;                 /* number of random deviates pairs */ 
    int iseed;             /* random generator seed */

    int nr;
    float dr,or,r;         /* correlation coefficient */
    float m1,s1;           /* mean and standard deviation */
    float m2,s2;           /* mean and standard deviation */

    float dbin,cs12;
    float inrm1,inrm2;
    float x1,y1,x2,y2;
    float **x,**y;         /* Gaussian distributed correlated random deviates */
    float *mc,**hist;

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

    if (!sf_getint("nbh",&nbh)) nbh=50;
    /* half number of bins for histogram */

    if (!sf_getfloat("r",&r)) r=0.7;
    /* histogram correlation coefficient */

    if (!sf_getfloat("m1",&m1)) m1=0.0;
    /* mean for deviate 1 */

    if (!sf_getfloat("s1",&s1)) s1=2.0;
    /* standard deviation for deviate 1 */

    if (!sf_getfloat("m2",&m2)) m2=2.0;
    /* mean for deviate 2 */

    if (!sf_getfloat("s2",&s2)) s2=1.0;
    /* standard deviation for deviate 2 */

    if (!sf_getfloat("dbin",&dbin)) dbin=0.1;
    /* histogram bin size */

    if (!sf_getint("iseed",&iseed)) iseed=-33;
    /* random generator seed */
    if (iseed >= 0) sf_error("Need strictly negative iseed");

    /* memory allocations */
    mc = sf_floatalloc(nr);

    /* read input */
    sf_floatread(mc,nr,in);

    /* Histogram */
    nbin = 2*nbh + 1;
    hist = sf_floatalloc2(nbin,nbin);
    for (k = 0; k < nbin; k++) {
	for (i = 0; i < nbin; i++) {
	    hist[k][i] = 0.0;
	}
    }
    hist_gauss(hist,n,r,m1,m2,s1,s2,nbh,dbin,&iseed);

    /* Correlated Gaussian deviates */
    x = sf_floatalloc2(2,n);
    y = sf_floatalloc2(2,n);

    for (k = 0; k < nr; k++) {

	r = or + k*dr;

        /* Generates n pairs of Gaussian distributed correlated deviates */
	for (i=0; i<n; i++) {
	    gauss_joint(x[i],r,m1,s1,m1,s1,&iseed);
	}
	for (i=0; i<n; i++) {
	    gauss_joint(y[i],r,m1,s1,m1,s1,&iseed);
	}

        /* Estimates <cos(2t)> = <2*cos^2(t)-1> */
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
