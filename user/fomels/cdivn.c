/* N-dimensional smooth division of complex numbers */
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

#include "cdivn.h"
#include "cweight.h"

static int niter, n;
static sf_complex *p;

void cdivn_init(int ndim   /* number of dimensions */, 
	       int nd     /* data size */, 
	       int *ndat  /* data dimensions [ndim] */, 
	       int *nbox  /* smoothing radius [ndim] */, 
	       int niter1 /* number of iterations */,
	       bool verb  /* verbosity */) 
/*< initialize >*/
{
    niter = niter1;
    n = nd;

    sf_ctrianglen_init(ndim, nbox, ndat);
    sf_cconjgrad_init(nd, nd, nd, nd, 1., 1.e-6, verb, false);
    p = sf_complexalloc (nd);
}

void cdivn_close (void)
/*< free allocated storage >*/
{
    sf_ctrianglen_close();
    sf_cconjgrad_close();
    free (p);
}

void cdivn (sf_complex* num, sf_complex* den,  sf_complex* rat)
/*< smoothly divide rat=num/den >*/
{
    cweight_init(den);
    sf_cconjgrad(NULL, cweight_lop,sf_ctrianglen_lop,p,rat,num,niter); 
}

/* 	$Id: divn.c 5595 2010-03-21 16:54:14Z sfomel $	 */
