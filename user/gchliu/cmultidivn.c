/* Smooth division with several components and complex numbers */
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

static sf_complex *p;

#include "cweight2.h"
#include "crepeat.h"

void cmultidivn_init(int nw            /* number of components */, 
		       int ndim          /* number of dimensions */, 
		       int n             /* data size */, 
		       int *ndat         /* data dimensions [ndim] */, 
		       int *nbox         /* smoothing radius [ndim] */,
		       sf_complex* den   /* denominator [nw*nd] */,
		       bool verb         /* verbosity flag */)
/*< initialize >*/
{
    int n2;

    n2 = n*nw;
    
    sf_ctrianglen_init(ndim, nbox, ndat);
    crepeat_init(n,nw,sf_ctrianglen_lop);

    sf_cconjgrad_init(n2, n2, n, n, 1., 1.e-6, verb, false);
    p = sf_complexalloc (n2);
    cweight2_init(nw,n,den);
}

void cmultidivn_close (void)
/*< free allocated storage >*/
{
    sf_ctrianglen_close();
    sf_cconjgrad_close();
    cweight2_close();
    free (p);
}

void cmultidivn (sf_complex* num  /* numerator */, 
		 sf_complex* rat  /* ratio */, 
		 int niter        /* number of iterations */)
/*< smoothly divide num/rat >*/
{
    sf_cconjgrad(NULL,cweight2_lop,crepeat_lop,p,rat,num,niter);
}

