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

#include "smoothder1.h"

static float *tmp;
static float ndati;

void smoothder_init(int n     /* data size */, 
		    int ndim  /* number of dimensions */, 
		    int *rect /* smoothing radius [ndim] */, 
		    int *ndat /* data size [ndim] */)
/*< initialize >*/
{
    int n1, n2;
    
    n1 = ndat[0];
    n2 = n/n1;
    
    sf_repeat_init(n1,n2,sf_causint_lop);
    sf_trianglen_init(ndim,rect,ndat);

    tmp = sf_floatalloc(n);
    //sf_conjgrad_init(n, n, n, n, 1., 1.e-8, true, false);
    ndati=n;
}

void smoothder_close(void)
/*< free allocated storage >*/
{
    free(tmp);
    sf_conjgrad_close();
    sf_trianglen_close();
}

void smoothder(int niter     /* number of iterations */, 
	       float* weight /* data weighting */, 
	       float* data   /* input data */, 
	       float* der    /* output derivative */) 
/*< find the derivative >*/
{ 
    if (NULL != weight) {
	sf_weight_init(weight);
	/*sf_conjgrad(sf_weight_lop,sf_repeat_lop,sf_trianglen_lop,
		    tmp,der,data,niter);*/

	sf_solver_prec(sf_repeat_lop,sf_cgstep,sf_trianglen_lop,
	    		ndati,ndati,ndati,der,data,niter,0.1,"verb",true,"end");

	sf_cgstep_close();

    } else {
	/*sf_conjgrad(NULL,sf_repeat_lop,sf_trianglen_lop,
		    tmp,der,data,niter);*/


    sf_solver_prec(sf_repeat_lop,sf_cgstep,sf_trianglen_lop,
    		ndati,ndati,ndati,der,data,niter,0.1,"verb",true,"end");

    sf_cgstep_close();





    }
}
