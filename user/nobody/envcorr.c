/* Correlation with the envelope. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include "envcorr.h"
#include "hilbert.h"
#include "divn.h"

static float **hlb, *den, *num, *rat2;
static int n1, ntr, nd;

void envcorr_init(int dim    /* dimensionality */, 
		  int *m     /* data dimensions [dim] */, 
		  int *rect  /* smoothing radius [dim] */, 
		  int niter  /* number of iterations */)
/*< initialize >*/
{
    int i;
    const int order=6;
    const float c=1.;

    n1 = m[0];
    nd = 1;
    for (i=0; i < dim; i++) {
	nd *= m[i];
    }
    ntr = nd/n1;
    
    hlb = sf_floatalloc2(n1,ntr);
    den = sf_floatalloc(nd);
    num = sf_floatalloc(nd);
    rat2 = sf_floatalloc(nd);

    hilbert_init(n1, order, c);
    divn_init(dim, nd, m, rect, niter);
}

void envcorr_close(void)
/*< free allocated storage >*/
{
    free(*hlb);
    free(hlb);
    free(num);
    free(den);
    free(rat2);
    hilbert_free();
    divn_close();
}

void envcorr(float** inp /* input data [ntr][n1] */, float* rat1)
/*< compute correlation >*/
{
    float doth, dout, *data, *hilb;
    int i1, i2, i;

    doth = 0.;
    dout = 0.;
    for (i2=0; i2 < ntr; i2++) {
	data = inp[i2];
	hilb = hlb[i2];

	hilbert(data,hilb);

	for (i1=0; i1 < n1; i1++) {
	    hilb[i1] = hypotf(data[i1],hilb[i1]);
	    dout += data[i1]*data[i1];
	    doth += hilb[i1]*hilb[i1];
	}
    }
    doth = sqrtf(nd/doth);
    dout = sqrtf(nd/dout);

    for (i2=0; i2 < ntr; i2++) {
	data = inp[i2];
	hilb = hlb[i2];

	for (i1=0; i1 < n1; i1++) {
	    i = i2*n1 + i1;
	    den[i] = data[i1]*dout;
	    num[i] = hilb[i1]*dout;
	}
    }

    divn(num,den,rat1);
	
    for (i2=0; i2 < ntr; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    i = i2*n1 + i1;
	    num[i] = inp[0][i1]*doth;
	    den[i] = hlb[0][i1]*doth;
	}
    }

    divn(num,den,rat2);
	
    for (i=0; i < nd; i++) {
	if (rat1[i] > 0.) {
	    if (rat2[i] > 0. || -rat2[i] < rat1[i]) {
		rat1[i] = sqrtf(fabsf(rat1[i]*rat2[i]));
	    } else {
		rat1[i] = -sqrtf(fabsf(rat1[i]*rat2[i]));
	    }
	} else {
	    if (rat2[i] < 0. || rat2[i] < -rat1[i]) {
		rat1[i] = -sqrtf(fabsf(rat1[i]*rat2[i]));
	    } else {
		rat1[i] = sqrtf(fabsf(rat1[i]*rat2[i]));
	    }
	}
    }
}

/* 	$Id: envcorr.c 744 2004-08-17 18:46:07Z fomels $	 */
