/* L1 Smooth division with several components (first attempt) */
/*
  Copyright (C) 2010 University Politecnico di Milano.
  
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

#include "helix.h"
/*^*/

#include "multidivnL1.h"
#include "trianglen.h"
#include "repeat.h"
#include "weight2.h"
#include "weight.h"
#include "error.h"
#include "conjgrad.h"
#include "alloc.h"
#include "helicon.h"
#include "sharpen.h"

static float *p,*nn,*rr,*num_tmp,*rat_tmp;
static int n,n2;
static bool prec,verb;

void sf_multidivnL1_init(int nw       /* number of components */,
		       int ndim     /* number of dimensions */, 
		       int n1        /* data size */,
		       int *ndat    /* data dimensions [ndim] */, 
		       int *nbox    /* smoothing radius [ndim] */,
		       float* den   /* denominator [nw*nd] */,
		       sf_filter aa /* data filter */,
		       float perc   /* percentage for sharpening */,
		       bool verb1    /* verbosity flag */)
/*< initialize >*/
{
    int i;
    verb=verb1;
    n = n1; /* data size*/
    n2 = n*nw; /* model size */

    sf_trianglen_init(ndim, nbox, ndat);
    sf_repeat_init(n,nw,sf_trianglen_lop);

    sf_conjgrad_init(n2, n2, n, n, 1., 1.e-6, verb, true);
    p = sf_floatalloc (n2);
    sf_weight2_init(nw,n,den);

    prec = (bool) (NULL != aa);
    if (prec) sf_helicon_init(aa);

    sf_sharpen_init(n,perc);
     /* initialization of sharpening regularization*/

    num_tmp = sf_floatalloc(n);
    nn = sf_floatalloc(n);
    rr = sf_floatalloc(n);

    rat_tmp = sf_floatalloc(n2);

    /* initialize with zero */
    for (i=0;i<n;i++) num_tmp[i]=nn[i]=rr[i]=0.0;

    for (i=0;i<n2;i++) rat_tmp[i]=p[i]=0.0;

}

void sf_multidivnL1_close (void)
/*< free allocated storage >*/
{
    sf_trianglen_close();
    sf_conjgrad_close();
    free (p);
    sf_weight2_close();
}

void sf_multidivnL1 (float* num  /* numerator [data]*/,
		   float* rat  /* ratio [model]*/,
		   int niter   /* number of POCS iterations */,
		   int liter   /* number of shaping CG iteration */)
/*< smoothly L1 divide num/rat >*/
{
	float en=0.0,eb=0.0;
	int i,iter;

    for (iter=0; iter < niter; iter++) {
    	/* Solve smooth division in L2 */
    	/* --------------------------- */
		eb=en=0.0;
		for (i=0; i < n; i++)
			rr[i] = num_tmp[i] = num[i]-nn[i]; /* dd = d - n */

		for (i=0; i < n2; i++)
			rat_tmp[i] = (-1) * rat[i];  /* -rat */

		sf_weight2_lop(false,true,n2,n,rat_tmp,rr);

    	sf_conjgrad(prec? sf_helicon_lop: NULL,
    			sf_weight2_lop,sf_repeat_lop,p,rat,num_tmp,liter);

		for (i=0; i < n2; i++) {
			rat_tmp[i]= -rat_tmp[i] - rat[i];  /* -d(rat) */
			eb+=rat[i]*rat[i];
		}
		for (i=0; i < n; i++)
			nn[i] += rr[i];

		sf_weight2_lop(false,true,n2,n,rat_tmp,nn); /* n[i] += r - M * drat;*/

		/* apply sharpening regularization*/
    	sf_sharpen(nn);
    	sf_weight_apply(n,nn);

		for (i=0; i < n; i++)
				en+= nn[i]*nn[i];

		if (verb) sf_warning("Sharpening iteration=%d: eb=%g en=%g",iter,eb,en);

    }

}


