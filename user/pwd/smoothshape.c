/* Smooth derivative by plane wave shaping regularization */
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

#include "smoothshape.h"
#include "pwdsl.h"

static float *tmp;

void smoothshape_init(int n1, int n2        /* data size */,
		      int order             /* accuracy order */,
		      int rect1, int rect2  /* smoothing radius */,
		      float lam             /* operator scaling */,
		      float **dip           /* dip field */) 
/*< initialize >*/
{
    int n;

    n = n1*n2;

    sf_repeat_init(n1,n2,sf_causint_lop);
    pwdsl_init(n1,n2,order,rect1,rect2, 0.01);
    pwdsl_set(dip);

    tmp = sf_floatalloc(n);
    sf_conjgrad_init(n, n, n, n, lam, 10*FLT_EPSILON, true, false);    
}

void smoothshape_close(void)
/*< free allocated storage >*/
{
    free(tmp);
    sf_conjgrad_close();
    pwdsl_close();
}

void smoothshape(int niter     /* number of iterations */, 
		 float* weight /* data weighting */, 
		 float* data   /* input data */, 
		 float* der    /* output derivative */) 
/*< find the derivative >*/
{ 
    if (NULL != weight) {
	sf_weight_init(weight);
	sf_conjgrad(sf_weight_lop,sf_repeat_lop,pwdsl_lop,tmp,der,data,niter);
    } else {
	sf_conjgrad(NULL,sf_repeat_lop,pwdsl_lop,tmp,der,data,niter);
    }
}
