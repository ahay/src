/* L1-like inversion */
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

#include "threshold.h"

static bool verb;
static int nd, niter;
static float *n, *r;
static char type;

void l1_init(int nd1    /* data size */,
	     int niter1 /* number of iterations */,
	     float perc /* thresholding percent */,
	     float fact /* thresholding factor */,
	     const char *type1 /* thresholding type */,
	     bool verb1  /* verbosity flag */)
/*< initialize >*/
{
    nd=nd1;
    niter=niter1;
    verb=verb1;
    type=type1[0];

    n = sf_floatalloc(nd);
    r = sf_floatalloc(nd);

    switch(type) {
	case 's':
	    sf_sharpen_init(nd,perc,fact);
	    break;
	case 't':
	    threshold_init(nd,fact);
	    break;
	default:
	    sf_error("unknown thresholding \"%s\"",type1);
    }
}

void l1_close(void)
/*< free allocated storage >*/
{
    free(n);
    free(r);

    switch(type) {
	case 's':
	    sf_sharpen_close();
	    break;
	case 't':
	    threshold_close();
	    break;
	default:
	    break;
    }
}

void unil1(const float *d, const float *a, float *alfa)
/*< uni-variate inversion >*/
{
    int i, iter;
    double ad, aa, a0, da;
    double alf;

    alf = 0.;
    for (i=0; i < nd; i++) {
	n[i] = 0.;
    }

    aa = cblas_dsdot( nd, a, 1, a, 1);
    a0 = 1.+aa;
    
    for (iter=0; iter < niter; iter++) {
	/* Solve |d - alpha * a - n|_2 */
	/* -------------------------------------- */
	for (i=0; i < nd; i++) {
	    r[i] = d[i]-n[i]-alf*a[i];
	}
	ad = cblas_dsdot( nd, a, 1, r, 1);
	da = ad/a0;
	alf += da;
	for (i=0; i < nd; i++) {
	    n[i] += r[i] - a[i]*da;
	}
	/* Threshold n */
	/* ----------- */

	switch(type) {
	    case 's':
		sf_sharpen(n);
		sf_weight_apply(nd,n);
		break;
	    case 't':
		threshold_set(n);
		threshold(n);
		break;
	    default:
		break;
	}

	if (verb) sf_warning("%d %g",iter,alf);
    }

    *alfa = alf;
}

void bil1(const float *d, const float *a, const float *b, 
	  float *alfa, float *beta)
/*< bi-variate inversion >*/
{
    int i, iter;
    double ad, bd, aa, bb, a0, b0, da, db, ab, det;
    double alf, bet;

    alf = bet = 0.;
    for (i=0; i < nd; i++) {
	n[i] = 0.;
    }

    aa = cblas_dsdot( nd, a, 1, a, 1);
    bb = cblas_dsdot( nd, b, 1, b, 1);
    ab = cblas_dsdot( nd, a, 1, b, 1);
    a0 = 1.+aa;
    b0 = 1.+bb;
    det = a0*b0-ab*ab;
    
    for (iter=0; iter < niter; iter++) {
	/* Solve |d - alpha * a - beta * b - n|_2 */
	/* -------------------------------------- */
	for (i=0; i < nd; i++) {
	    r[i] = d[i]-n[i]-alf*a[i]-bet*b[i];
	}
	ad = cblas_dsdot( nd, a, 1, r, 1);
	bd = cblas_dsdot( nd, b, 1, r, 1);
	da = (b0*ad-ab*bd)/det;
	db = (a0*bd-ab*ad)/det;
	alf += da;
	bet += db;
	for (i=0; i < nd; i++) {
	    n[i] += r[i] - a[i]*da - b[i]*db;
	}
	/* Threshold n */
	/* ----------- */

	switch(type) {
	    case 's':
		sf_sharpen(n);
		sf_weight_apply(nd,n);
		break;
	    case 't':
		threshold_set(n);
		threshold(n);
		break;
	    default:
		break;
	}

	if (verb) sf_warning("%d %g %g",iter,alf,bet);
    }

    *alfa = alf;
    *beta = bet;
}
