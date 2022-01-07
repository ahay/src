/* Smooth division with several components and complex numbers
e.g., regularized non-stationary autoregression */
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
#include "cdivn_rnar.h"
#include "cdivn.h" /*for the repeat_smooth and weight operations*/


static int *n, s[SF_MAX_DIM], nd, dim;
static sf_ctriangle *tr;
static sf_complex *tmp;


static sf_complex *p;

void cmultidivn_rnar_init(int nw            /* number of components */, 
		       int ndim          /* number of dimensions */, 
		       int n             /* data size */, 
		       int *ndat         /* data dimensions [ndim] */, 
			   int *nbox /* triangle radius [ndim-1] */, 
		       sf_complex* den   /* denominator [nw*nd] */,
		       bool verb         /* verbosity flag */)
/*< initialize cmultidivn with frequency-dependent smoothing >*/
{
    int n2;

    n2 = n*nw;
    smooth_rnar_init(ndim, nbox, ndat);
    crepeat_init(n,nw,smooth_rnar_lop);
    sf_cconjgrad_init(n2, n2, n, n, 1., 1.e-6, verb, false);
    p = sf_complexalloc (n2);
    cweight2_init(nw,n,den);
}

void cmultidivn_rnar_close (void)
/*< free allocated storage >*/
{
	smooth_rnar_close();
    sf_cconjgrad_close();
    cweight2_close();
    free (p);
}

void cmultidivn_rnar (sf_complex* num  /* numerator */, 
		 sf_complex* rat  /* ratio */, 
		 int niter        /* number of iterations */)
/*< smoothly divide num/rat >*/
{
    sf_cconjgrad(NULL,cweight2_lop,crepeat_lop,p,rat,num,niter);
}


void smooth_rnar_init (int ndim  /* number of dimensions */, 
			int *nbox /* triangle radius [ndim-1] */, 
			int *ndat /* data dimensions [ndim] */)
/*< initialize >*/
{
    int i;

    dim = ndim;
    n = sf_intalloc(dim);
    
    tr = (sf_ctriangle*) sf_alloc(dim,sizeof(sf_ctriangle));
    nd = 1;
    for (i=0; i < dim; i++) {
	tr[i] = (nbox[i] > 1)? sf_ctriangle_init (nbox[i],ndat[i],false): NULL;
	s[i] = nd;
	n[i] = ndat[i];
	nd *= ndat[i];
    }

    tmp = sf_complexalloc (nd);
}

void smooth_rnar_close(void)
/*< free allocated storage >*/
{
    int i;

	free(n);
    free (tmp);

    for (i=0; i < dim; i++) {
	if (NULL != tr[i]) sf_ctriangle_close (tr[i]);
    }

    free(tr);
}


void smooth_rnar_lop (bool adj, bool add, int nx, int ny, sf_complex* x, sf_complex* y)
/*< linear operator >*/
{
    int i, j, i0;

    if (nx != ny || nx != nd) 
	sf_error("%s: Wrong data dimensions: nx=%d, ny=%d, nd=%d",
		 __FILE__,nx,ny,nd);

    sf_cadjnull (adj,add,nx,ny,x,y);
  
    if (adj) {
	for (i=0; i < nd; i++) {
	    tmp[i] = y[i];
	}
    } else {
	for (i=0; i < nd; i++) {
	    tmp[i] = x[i];
	}
    }

    /*only for x and y when rect[0]=1 (rect1=1 in command line)*/
    for (i=0; i < dim; i++) {
	if (NULL != tr[i]) {
	    for (j=0; j < nd/n[i]; j++) {
		i0 = sf_first_index (i,j,dim,n,s);
		sf_csmooth (tr[i], i0, s[i], false, tmp);
	    }
	}
    }
	
    if (adj) {
	for (i=0; i < nd; i++) {
#ifdef SF_HAS_COMPLEX_H
	    x[i] += tmp[i];
#else
	    x[i] = sf_cadd(x[i],tmp[i]);
#endif
	}
    } else {
	for (i=0; i < nd; i++) {
#ifdef SF_HAS_COMPLEX_H
	    y[i] += tmp[i];
#else
	    y[i] = sf_cadd(y[i],tmp[i]);
#endif
	}
    }    
}



