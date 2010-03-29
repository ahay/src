/* Smooth division with several components and nonstationary smoothing */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "ntriangle2.h"

static float *p;

void nmultidivn_init(int nw       /* number of components */, \
		     int ndim     /* number of dimensions */,
		     int n        /* data size */,
		     int *ndat    /* data dimensions [ndim] */, 
		     int *nbox    /* smoothing radius [nw] */,
		     float* den   /* denominator [nw*nd] */,
		     bool verb    /* verbosity flag */)
/*< initialize >*/
{
    int n2;

    n2 = n*nw;
    
    ntrianglen2_init(nw,n,ndat[0],nbox);

    sf_conjgrad_init(n2, n2, n, n, 1., 1.e-6, verb, false);
    p = sf_floatalloc (n2);
    sf_weight2_init(nw,n,den);
}

void nmultidivn_close (void)
/*< free allocated storage >*/
{
    ntriangle2_close();
    
    sf_conjgrad_close();
    free (p);
    sf_weight2_close();
}

void nmultidivn (float* num  /* numerator */, 
		   float* rat  /* ratio */, 
		   int niter   /* number of iterations */)
/*< smoothly divide num/rat >*/
{
    sf_conjgrad(NULL,sf_weight2_lop,ntriangle2_lop,p,rat,num,niter);
}

/* 	$Id: multidivn.c 1136 2005-04-20 20:43:14Z fomels $	 */
