/* Smooth division with several components */
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

#include "helix.h"
/*^*/

#include "multidivn.h"
#include "trianglen.h"
#include "repeat.h"
#include "weight2.h"
#include "conjgrad.h"
#include "alloc.h"
#include "helicon.h"

static float *p;
static bool prec;

void sf_multidivn_init(int nw       /* number of components */, 
		       int ndim     /* number of dimensions */, 
		       int n        /* data size */, 
		       int *ndat    /* data dimensions [ndim] */, 
		       int *nbox    /* smoothing radius [ndim] */,
		       float* den   /* denominator [nw*nd] */,
		       sf_filter aa /* data filter */,
		       bool verb    /* verbosity flag */)
/*< initialize >*/
{
    int n2;

    n2 = n*nw;
    
    sf_trianglen_init(ndim, nbox, ndat);
    sf_repeat_init(n,nw,sf_trianglen_lop);

    sf_conjgrad_init(n2, n2, n, n, 1., 1.e-6, verb, false);
    p = sf_floatalloc (n2);
    sf_weight2_init(nw,n,den);

    prec = (bool) (NULL != aa);
    if (prec) sf_helicon_init(aa);
}

void sf_multidivn_close (void)
/*< free allocated storage >*/
{
    sf_trianglen_close();
    sf_conjgrad_close();
    free (p);
    sf_weight2_close();
}

void sf_multidivn (float* num  /* numerator */, 
		   float* rat  /* ratio */, 
		   int niter   /* number of iterations */)
/*< smoothly divide num/rat >*/
{
    sf_conjgrad(prec? sf_helicon_lop: NULL,
		sf_weight2_lop,sf_repeat_lop,p,rat,num,niter);
}

void sf_multidivn_adj (bool adj /* adjoint flag */,
		       float* num  /* numerator */, 
		       float* rat  /* ratio */, 
		       int niter   /* number of iterations */)
/*< smoothly divide num/rat >*/
{
    sf_conjgrad_adj(adj,sf_weight2_lop,sf_repeat_lop,p,rat,num,niter);
}

/* 	$Id: multidivn.c 1136 2005-04-20 20:43:14Z fomels $	 */
