/* Combinations */
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

static void comb_next (int n, int k, int *a, int *done)
/* Computes combinations of k elements out of n, computed one at a time, in lexicographical order. */
/* This code is distributed under the GNU LGPL license. Modified from original author John Burkardt */
{
    int i;
    int j;

    if (*done) {

	for (i = 0; i < k; i++) a[i] = i+1;
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
		for (j = i; j < (k+1); j++) a[j-1] = a[i-2] + j - (i-1);
		return;
	    }
	}

	*done = 1;

    }
    return;
}

int main(int argc, char* argv[])
{
    int n,k,nc,i,j;

    int done = 1;   /* boolean for more combinations to compute */
    int *a;         /* list of elements in the current combination (not needed at startup) */
    int **c;

    sf_axis acb, ake;
    sf_file in,out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getint("k",&k)) sf_error("Need k=");
    /* combination of k elements */

    /* input file */
    if (!sf_histint(in,"n1",&n)) sf_error("No n1=");

    nc = binomial(n,k);
    sf_warning("Number of combinations is %3d",nc);

    /* output file parameters */
    acb = sf_maxa(nc,0,1);
    sf_oaxa(out,acb,1);

    ake = sf_maxa(k,0,1);
    sf_oaxa(out,ake,2);

    sf_putstring (out,"label1", "combination");
    sf_putstring (out,"label2", "elements");

    /* memory allocations */
    a = sf_intalloc(k);
    c = sf_intalloc2(nc,k);

    j = 0;

    while (done) {
        /* Combination of k elements out of n */
	comb_next(n,k,a,&done);
        /* done = 1 as long as there are more combinations to compute */
        /* done = 0 when the list is exhausted. */
	for (i = 0; i < k; i++) {
	    sf_warning("  %3d",a[i]);
	    c[j][i] = a[i];
	}
	j++;
    }

    /* output */ 
    sf_intwrite (c[0],k*nc,out);

    exit(0);
}

