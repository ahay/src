/* Smooth derivative by shaping regularization with non-stationary smoothing */

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

static float *tmp;

void nsmoothder_init(int nd     /* data size */, 
        int ndim  /* number of dimensions */, 
        int *nbox /* smoothing radius [ndim] */, 
        float **rct /* triangle lengths [ndim][nd] */,
        float **sft /* triangle shifts [ndim][nd] */,
        int *ndat /* data size [ndim] */,
        int n1_v1, 
        int n2_v1)
/*< initialize >*/
{
    int n1, n2;
    n1 = ndat[0];
    n2 = nd/n1;
    
    sf_repeat_init(n1,n2,sf_causint_lop);

    sf_ntrianglen_init(ndim,nbox,ndat,rct,sft,1);

    tmp = sf_floatalloc(nd);
    sf_conjgrad_init(nd, nd, nd, nd, 1., 1.e-8, true, false);    
}

void nsmoothder_close(void)
/*< free allocated storage >*/
{
    free(tmp);
    sf_conjgrad_close();
    sf_trianglen_close();
}

void nsmoothder(int niter     /* number of iterations */, 
         float* weight /* data weighting */, 
         float* data   /* input data */, 
         float* der    /* output derivative */) 
/*< find the derivative >*/
{ 
    if (NULL != weight) {
  sf_weight_init(weight);
  sf_conjgrad(sf_weight_lop,sf_repeat_lop,sf_ntrianglen_lop,
        tmp,der,data,niter);
    } else {
  sf_conjgrad(NULL,sf_repeat_lop,sf_ntrianglen_lop,
        tmp,der,data,niter);
    }
}

