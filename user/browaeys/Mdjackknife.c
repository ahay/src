/* D-delete jackknife estimate of variance. */
/* Leave out d observations with sqrt(n) < d < n-1
   Jackknife samples where d elements removed are
   x(1),...,x(i-1),x(i+1),...,x(n) */
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

static int binomial (int n, int k)
/* Computes the binomial coefficient C(n,k) = n! / k!*(n-k)! */
/* This code is distributed under the GNU LGPL license. Modified from original author John Burkardt */
{
    int i, icnk, mn, mx;

    mn = SF_MIN(k,(n-k));

    if (mn < 0) {

	icnk = 0;

    } else if (mn == 0) {

	icnk = 1;

    } else {

	mx = SF_MAX(k,(n-k));
	icnk = mx + 1;

	for (i = 2; i < (mn+1); i++) {
	    icnk = (icnk*(mx + i))/i;
	}
    }

    return(icnk);
}

static void comb_next (int n, int k, int *a, int *done )
/* Computes combinations of k elements out of n, computed one at a time, in lexicographical order. */
/* This code is distributed under the GNU LGPL license. Modified from original author John Burkardt */
{
    int i;
    int j;

    if (*done) {

	for (i = 0; i < k; i++) {
	    a[i] = i+1;
	}

	if ( 1 < k ) {
	    *done = 0;
	} else {
	    *done = 1;
	}

    } else {

	if (a[k-1] < n) {

	    a[k-1] += 1;
	    return;

	}

	for (i = k; i > 1; i--) {

	    if (a[i-2] < (n-k+i-1)) {

		a[i-2] += 1;
		for (j = i; j < (k+1); j++) {
		    a[j-1] = a[i-2] + j - (i-1);
		}
		return;

	    }
	}

	*done = 1;

    }
    return;
}


int main(int argc, char* argv[])
{
    int i,j,k,ie,nrm;
    int n,d,ns,ss;
    int done;          /* boolean for more combinations to compute */
    int *a;            /* list of elements in the current combination (not needed at startup) */

    float ms, mt, sp, sedj;
    float *data, *xx, *sj;

    sf_file in,out;

    sf_init (argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n)) sf_error("No n1=");
    /* number of data */

    /* memory allocations */
    data = sf_floatalloc(n);

    /* read input data */
    sf_floatread(data,n,in);

    /* d-delete jackknife number */
    d = floor(sqrtf(n));

    /* samples size  */
    ss = n - d;

    /* number of samples (combinations) = n!/(d!(n-d)!) */
    ns = binomial(n,d);
    sf_warning("Number of combinations is  %3d",ns);

    /* current jacknife sample */
    xx = sf_floatalloc(ss);
    for (i = 0; i < ss; i++) {
	xx[i] = 0.0;
    }

    /* statistical estimator jackknife replication */
    sj = sf_floatalloc(ns);
    for (i = 0; i < ns; i++) {
	sj[i] = 0.0;
    }

    i = 0;
    a = sf_intalloc(d);
    done = 1;
    while (done) {

        /* Combination of d elements out of n */
	comb_next(n,d,a,&done);
        /* done = 1 as long as there are more combinations to compute */
        /* done = 0 when the list is exhausted. */
	i += 1;

        /* Extract jackknife samples from data */
	ie = 0;
	for (j = 0; j < n; j++) {
	    nrm = 1;
	    for (k = 0; k < d; k++) {
		if (j==a[k]) nrm = 0;
	    }
	    if (nrm) {
		xx[ie] = data[j];
		ie += 1;
	    }
	}
    
        /* Mean of jackknife sample */
	sp = 0.0;
	for (j = 0; j < ss; j++) {
	    sp += xx[j];
	}
	sp /= ss;

        /* variance of jackknife sample */
	for (j = 0; j < ss; j++) {
	    mt = xx[j]-sp;
	    sj[i] += mt*mt;
	}
	sj[i] /= ss;

    }

    sf_warning("Number of combinations performed is  %3d",n);

    /* Expectation of jackknife replications sj of statistical estimator */
    sp = 0.0;
    for (i = 0; i < ns; i++) {
	sp += sj[i];
    }
    sp /= ns;

    /* Delete-d jackknife estimate of standard error */
    ms = 0.0;
    for (i = 0; i < ns; i++) {
	mt = sj[i]-sp;
	ms += mt*mt;
    }
    sedj = sqrtf(ms*ss/(d*ns));
    
    exit(0);
}


