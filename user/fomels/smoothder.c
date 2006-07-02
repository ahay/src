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
#include "ntrianglen.h"
#include "weight.h"

static float *tmp;
static bool nonstat;
static sf_operator oper;

void smoothder_init(int n     /* data size */, 
		    int ndim  /* number of dimensions */, 
		    int *rect /* smoothing radius [ndim] */, 
		    int *ndat /* data size [ndim] */,
		    int **len /* radius for nonstat smoothing [ndim][nd] */)
/*< initialize >*/
{
    int n1, n2;
    
    n1 = ndat[0];
    n2 = n/n1;
    
    repeat_init(n1,n2,sf_causint_lop);

    nonstat = (bool) (NULL != len);

    if (nonstat) {
	ntrianglen_init(ndim,rect,ndat,len);
	oper = ntrianglen_lop;
    } else {
	trianglen_init(ndim,rect,ndat);
	oper = trianglen_lop;
    }

    tmp = sf_floatalloc(n);
    sf_conjgrad_init(n, n, n, n, 1., 1.e-8, true, false);    
}

void smoothder_close(void)
/*< free allocated storage >*/
{
    free(tmp);
    sf_conjgrad_close();
    if (nonstat) {
	ntrianglen_close();
    } else {
	trianglen_close();
    }
}

void smoothder(int niter     /* number of iterations */, 
	       float* weight /* data weighting */, 
	       float* data   /* input data */, 
	       float* der    /* output derivative */) 
/*< find the derivative >*/
{ 
    if (NULL != weight) {
	weight_init(weight);
	sf_conjgrad(weight_lop,repeat_lop,oper,tmp,der,data,niter);
    } else {
	sf_conjgrad(NULL,repeat_lop,oper,tmp,der,data,niter);
    }
}
