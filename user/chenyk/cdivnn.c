/* Smooth division with several components and complex numbers
with frequency-dependent smoothing */
/*
  Copyright (C) 2020 Yangkang Chen
  
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
#include "cdivnn.h"
#include "cdivn.h" /*for the repeat_smooth and weight operations*/
#include "cntrianglen.h"

static int n2, nw;
static sf_complex *p;

void cmultidivnn_init(int nw            /* number of components */, 
		       int ndim          /* number of dimensions */, 
		       int n             /* data size */, 
		       int *ndat         /* data dimensions [ndim] */, 
		       int *nbox         /* smoothing radius [ndim] */,
		       float **rct /* triangle lengths [ndim][nd] */,
               int **sft /* triangle shifts [ndim][nd] */,
		       sf_complex* den   /* denominator [nw*nd] */,
		       bool verb         /* verbosity flag */)
/*< initialize cmultidivn with frequency-dependent smoothing >*/
{
    int n2;

    n2 = n*nw;
    
    cntrianglen_init(ndim, nbox, ndat, rct, sft, 1);
    crepeat_init(n,nw,cntrianglen_lop);

    sf_cconjgrad_init(n2, n2, n, n, 1., 1.e-6, verb, false);
    p = sf_complexalloc (n2);
    cweight2_init(nw,n,den);
}

void cmultidivnn_close (void)
/*< free allocated storage >*/
{
	cntrianglen_close();
    sf_cconjgrad_close();
    cweight2_close();
    free (p);
}

void cmultidivnn (sf_complex* num  /* numerator */, 
		 sf_complex* rat  /* ratio */, 
		 int niter        /* number of iterations */)
/*< smoothly divide num/rat >*/
{
    sf_cconjgrad(NULL,cweight2_lop,crepeat_lop,p,rat,num,niter);
}



/**Following is single division*/

static int niter, n;
// static sf_complex *p; /*defined above*/

void cdivnn_init(int ndim   /* number of dimensions */, 
	       int nd     /* data size */, 
	       int *ndat  /* data dimensions [ndim] */, 
	       int *nbox  /* smoothing radius [ndim] */, 
		   float **rct /* triangle lengths [ndim][nd] */,
           int **sft /* triangle shifts [ndim][nd] */,
	       int niter1 /* number of iterations */,
	       bool verb  /* verbosity */) 
/*< initialize >*/
{
    niter = niter1;
    n = nd;
	cntrianglen_init(ndim, nbox, ndat, rct, sft, 1);
//     sf_ctrianglen_init(ndim, nbox, ndat);
    sf_cconjgrad_init(nd, nd, nd, nd, 1., 1.e-6, verb, false);
    p = sf_complexalloc (nd);
}

void cdivnn_close (void)
/*< free allocated storage >*/
{
    cntrianglen_close();
    sf_cconjgrad_close();
    free (p);
}

void cdivnn (sf_complex* num, sf_complex* den,  sf_complex* rat)
/*< smoothly divide rat=num/den >*/
{
    cweight_init(den);
    sf_cconjgrad(NULL, cweight_lop,cntrianglen_lop,p,rat,num,niter); 
//     sf_cconjgrad(NULL, cweight_lop,sf_ctrianglen_lop,p,rat,num,niter); 

}


void cdivnne (sf_complex* num, sf_complex* den,  sf_complex* rat)
/*< smoothly divide rat=num/den with preconditioning >*/
{
    int id;
    float a, norm=0.; 

    for (id=0; id < n; id++) {
	a = cabsf(den[id]);
	norm += a*a;
    }
    norm = sqrtf(n/norm);

    for (id=0; id < n; id++) {
#ifdef SF_HAS_COMPLEX_H
	num[id] *= norm;
	den[id] *= norm;
#else
	num[id] = sf_crmul(num[id],norm);
	den[id] = sf_crmul(sig[id],norm);
#endif
    }
    
    cweight_init(den);
    sf_cconjgrad(NULL, cweight_lop,cntrianglen_lop,p,rat,num,niter); 
}

