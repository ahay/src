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

#include "divn.h"
#include "trianglen.h"

static int niter;
static float *p;

void divn_init(int ndim   /* number of dimensions */, 
	       int nd     /* data size */, 
	       int *ndat  /* data dimensions [ndim] */, 
	       int *nbox  /* smoothing radius [ndim] */, 
	       int niter1 /* number of iterations */,
	       bool verb  /* verbosity */) 
/*< initialize >*/
{
    niter = niter1;

    trianglen_init(ndim, nbox, ndat);
    sf_conjgrad_init(nd, nd, nd, nd, 1., 1.e-6, verb, false);
    p = sf_floatalloc (nd);
}

void divn_close (void)
/*< free allocated storage >*/
{
    trianglen_close();
    sf_conjgrad_close();
    free (p);
}

void divn (float* num, float* den,  float* rat)
/*< smoothly divide rat=num/den >*/
{
    sf_weight_init(den);
    sf_conjgrad(NULL, sf_weight_lop,trianglen_lop,p,rat,num,niter); 
}

/* 	$Id: divn.c 4675 2009-08-24 16:34:57Z sfomel $	 */
