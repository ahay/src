/* N-dimensional smooth division */
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

#include "weight2.h"

static int niter, n;
static float *p;

void divn2_init(int ndim   /* number of dimensions */, 
		  int nd     /* data size */, 
		  int *ndat  /* data dimensions [ndim] */, 
		  int *nbox  /* smoothing radius [ndim] */, 
		  int niter1 /* number of iterations */,
		  bool verb  /* verbosity */) 
/*< initialize >*/
{
    niter = niter1;
    n = nd;

    sf_trianglen_init(ndim, nbox, ndat);
    sf_conjgrad_init(nd, nd, 2*nd, 2*nd, 1., 1.e-6, verb, false);
    p = sf_floatalloc (nd);
}

void divn2_close (void)
/*< free allocated storage >*/
{
    sf_trianglen_close();
    sf_conjgrad_close();
    free (p);
}

void divn2 (float** num, float** den,  float* rat)
/*< smoothly divide rat=num/den >*/
{
    weight2_init(den);
    sf_conjgrad(NULL, weight2_lop,sf_trianglen_lop,p,rat,num[0],niter); 
}

void divn2_enhance (float* rat)
/*< enhance contast >*/
{
    int i;
    float p;

    for (i=0; i < n; i++) {
	p = fabsf(rat[i]);
	p += 1.;
	p *= p;
	p *= p/16.;
	rat[i] = p;	
    }
}
   
