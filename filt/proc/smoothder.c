/* Smooth derivative by shaping regularization */
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

#include "smoothder.h"

#include "repeat.h"
#include "trianglen.h"
#include "trisl.h"
#include "weight.h"
#include "impl2.h"

static int n, n2;
static float **tmp, *tmp1, *tmp2;
static bool diff, dip;

int smoothder_init(int ndim   /* number of dimensions */, 
		   int *rect  /* smoothing radius [ndim] */, 
		   int *ndat  /* data size [ndim] */, 
		   bool diff1 /* use anisotropic diffusion */, 
		   bool dip1  /* use slope */) 
/*< initialize >*/
{
    int i, n1;
    
    n=1;
    for (i=0; i <ndim; i++) {
	n *= ndat[i];
    }

    n1 = ndat[0];
    n2 = n/n1;
    
    repeat_init(n1,n2,sf_causint_lop);
    trianglen_init(ndim,rect,ndat);

    tmp = sf_floatalloc2(n1,n2);
    tmp1 = tmp[0];
    tmp2 = sf_floatalloc(n);

    diff = diff1;
    dip = dip1;

    sf_conjgrad_init(n, n, n, n, 1., 1.e-8, true, diff);    

    if (dip) {
	trisl_init(n1,n2,rect[0],rect[1]);
    } else if (diff) {
	impl2_init ((float) rect[0], (float) rect[1], 
		    n1, n2, 1., 50., false);
    }

    return n;
}

void smoothder_close(void)
/*< free allocated storage >*/
{
    free(tmp1);
    free(tmp);
    sf_conjgrad_close();
    trianglen_close();
    if (diff) impl2_close();
    if (dip) trisl_close();
}

void smoothder(int niter     /* number of iterations */, 
	       float* weight /* data weighting */, 
	       float* data   /* input data */, 
	       float* der    /* output derivative */) 
/*< find the derivative >*/
{
 
    if (NULL != weight) {
	weight_init(weight);
	sf_conjgrad(weight_lop,repeat_lop,trianglen_lop,tmp2,der,data,niter);
    } else {
	sf_conjgrad(NULL,repeat_lop,trianglen_lop,tmp2,der,data,niter);
    }
}

void smoothdip(int niter     /* number of iterations */, 
	       float** dip   /* slope */, 
	       float* weight /* data weighting */, 
	       float* data   /* input data */, 
	       float* der    /* output derivative */) 
/*< find the derivative along slope >*/
{
    trisl_set(dip);

    if (NULL != weight) {
	weight_init(weight);
	sf_conjgrad(weight_lop,repeat_lop,trisl_lop,tmp2,der,data,niter);
    } else {
	sf_conjgrad(NULL,repeat_lop,trisl_lop,tmp2,der,data,niter);
    }
}


void smoothdiff(int niter     /* number of iterations */, 
		int ncycle    /* number of cycles */, 
		float* weight /* data weighting */, 
		float* data   /* input data */, 
		float* der    /* output derivative */) 
/*< find the derivative using anisotropic diffusion >*/
{
    int i, iter;

    for (i=0; i < n; i++) {
	tmp2[i] = 0.;
    }

    if (NULL != weight) {
	weight_init(weight);
	sf_conjgrad(weight_lop,repeat_lop,impl2_lop,tmp2,der,data,niter);
    } else {
	sf_conjgrad(NULL,repeat_lop,impl2_lop,tmp2,der,data,niter);
    }
	
    for (iter=0; iter < ncycle; iter++) {
	/*
	for (i=0; i < n; i++) {
	    tmp1[i] = der[i];
	}
	impl2_set(tmp);
	*/
	sf_conjgrad((NULL!=weight)? weight_lop: NULL,
		    repeat_lop,impl2_lop,tmp2,der,data,niter/ncycle);
    }    
}

