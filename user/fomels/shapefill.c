/* 2-D missing data interpolation by shaping regularization */
/*
  Copyright (C) 2014 University of Texas at Austin
  
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
/*^*/

static int n;
static float *p;

void shapefill_init (int ndim   /* number of dimensions */, 
		     int nd     /* data size */, 
		     int *ndat  /* data dimensions [ndim] */, 
		     int *nbox  /* smoothing radius [ndim] */,
		     bool verb  /* verbosity */)
/*< initialize >*/
{
    n = nd;
    
    sf_trianglen_init(ndim, nbox, ndat);
    sf_conjgrad_init(nd, nd, nd, nd, 1., 1.e-6, verb, false);
    p = sf_floatalloc (nd);
}

void shapefill_close (void)
/*< free allocated storage >*/
{
    sf_trianglen_close();
    sf_conjgrad_close();
    free (p);
}

void shapefill(int niter   /* number of iterations */, 
	       float* mm   /* model */, 
	       bool *known /* mask for known data */)
/*< interpolate >*/
{
    sf_mask_init(known);
    sf_conjgrad(NULL,sf_mask_lop,sf_trianglen_lop,p,mm,mm,niter); 
}
