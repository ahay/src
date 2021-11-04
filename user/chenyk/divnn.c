/* N-dimensional smooth division with non-stationary smoothing*/
/*
  Copyright (C) 2020 University of Texas at Austin
  By Yangkang Chen, On Jan, 2020
  
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
#include "divnn.h"

static int niter, n;
static float *p;

void divnn_init(int ndim   /* number of dimensions */, 
		  int nd     /* data size */, 
		  int *nbox /* triangle radius [ndim] */, 
		  int *ndat  /* data dimensions [ndim] */, 
		  float **rct /* triangle lengths [ndim][nd] */,
          int **sft /* triangle shifts [ndim][nd] */,
		  int niter1 /* number of iterations */,
		  bool verb  /* verbosity */) 
/*< initialize >*/
{

    niter = niter1;
    n = nd;
	
	/*initialization for the non-stationary triangle smoothing operator*/
    sf_ntrianglen_init(ndim,nbox,ndat,rct,sft,1);
	
    sf_conjgrad_init(nd, nd, nd, nd, 1., 1.e-6, verb, false);
    p = sf_floatalloc (nd);
}

void divnn_close (void)
/*< free allocated storage >*/
{
    sf_trianglen_close();
    sf_conjgrad_close();
    free (p);
}

void divnn (float* num, float* den,  float* rat)
/*< smoothly divide rat=num/den >*/
{
    sf_weight_init(den);
    sf_conjgrad(NULL, sf_weight_lop,sf_ntrianglen_lop,p,rat,num,niter); 
}

void divnne (float* num, float* den,  float* rat, float eps)
/*< smoothly divide rat=num/den with preconditioning >*/
{
    int i;
    double norm;

    if (eps > 0.0f) {
	for (i=0; i < n; i++) {
	    norm = 1.0/hypot(den[i],eps);

	    num[i] *= norm;
	    den[i] *= norm;
	}
    } 

    norm = cblas_dsdot(n,den,1,den,1);
    if (norm == 0.0) {
	for (i=0; i < n; i++) {
	    rat[i] = 0.0;
	}
	return;
    }
    norm = sqrt(n/norm);

    for (i=0; i < n; i++) {
	num[i] *= norm;
	den[i] *= norm;
    }   

    sf_weight_init(den);
    sf_conjgrad(NULL, sf_weight_lop,sf_ntrianglen_lop,p,rat,num,niter); 
}


void divnn_combine (const float* one, const float* two, float *prod)
/*< compute product of two divisions >*/
{
    int i;
    float p;

    for (i=0; i < n; i++) {
	p = sqrtf(fabsf(one[i]*two[i]));
	if ((one[i] > 0. && two[i] < 0. && -two[i] >= one[i]) ||
	    (one[i] < 0. && two[i] > 0. && two[i] >= -one[i])) 
	    p = -p;
	p += 1.;
	p *= p;
	p *= p/16.;
	prod[i] = p;	
    }
}

void divnn_combine_sign (const float* one, const float* two, float *prod)
/*< compute product of two divisions (with a sign) >*/
{
    int i;
    float p;

    for (i=0; i < n; i++) {
	p = sqrtf(fabsf(one[i]*two[i]));
	if (one[i] < 0. || two[i] < 0.) 
	    p = -p;
	prod[i] = p;
    }
}
    


