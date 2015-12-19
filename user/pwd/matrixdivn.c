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
#include <rsf.h>

#include "matrixdivn.h"
#include "twobytwo.h"

static int niter;
static float *p;

void matrixdivn_init(int ndim     /* number of dimensions */, 
		     int n        /* data size */, 
		     int *ndat    /* data dimensions [ndim] */, 
		     int *nbox    /* smoothing radius [ndim] */,
		     float** mat  /* matrux [4][nd] */,
		     int niter1   /* number of iterations */,
		     bool verb    /* verbosity flag */)
/*< initialize >*/
{
    int n2;

    niter = niter1;
    n2 = n*2;
    
    sf_trianglen_init(ndim, nbox, ndat);
    sf_repeat_init(n,2,sf_trianglen_lop);

    sf_conjgrad_init(n2, n2, n, n, 1., 1.e-6, verb, false);
    p = sf_floatalloc (n2);
    twobytwo_init(mat);
}

void matrixdivn_close (void)
/*< free allocated storage >*/
{
    sf_trianglen_close();
    sf_conjgrad_close();
    free (p);
}

void matrixdivn (float* dat  /* data */, 
		 float* mod /* model */)
/*< smoothly invert >*/
{
    sf_conjgrad(NULL,twobytwo_lop,sf_repeat_lop,p,mod,dat,niter);
}

