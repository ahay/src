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
#include "cdivn_fs.h"
#include "cdivn.h" /*for the repeat_smooth and weight operations*/

static int *n, s[SF_MAX_DIM], nd, dim; 
static sf_ctriangle tr0; /*triangle smoother in frequency*/
static sf_ctriangle *tr;
static sf_complex *tmp;

static sf_complex *p;
static int nf; /*nf is the dimension of frequency*/

void cmultidivn_fs_init(int nw            /* number of components */, 
		       int ndim          /* number of dimensions */, 
		       int n             /* data size */, 
		       int *ndat         /* data dimensions [ndim] */, 
			   int nbox0  /* triangle radius in frequency */, 
			   int **nbox /* triangle radius [ndim-1] */, 
		       sf_complex* den   /* denominator [nw*nd] */,
		       bool verb         /* verbosity flag */)
/*< initialize cmultidivn with frequency-dependent smoothing >*/
{
    int n2;

    nf = ndat[0];
    n2 = n*nw;
    smooth_fs_init(ndim, nbox0, nbox, ndat);
    crepeat_init(n,nw,smooth_fs_lop);
    sf_cconjgrad_init(n2, n2, n, n, 1., 1.e-6, verb, false);
    p = sf_complexalloc (n2);
    cweight2_init(nw,n,den);    
}


void cmultidivn_fs_close (void)
/*< free allocated storage >*/
{
	smooth_fs_close();
    sf_cconjgrad_close();
    cweight2_close();
    free (p);
}

void cmultidivn_fs (sf_complex* num  /* numerator */, 
		 sf_complex* rat  /* ratio */, 
		 int niter        /* number of iterations */)
/*< smoothly divide num/rat >*/
{
    sf_cconjgrad(NULL,cweight2_lop,crepeat_lop,p,rat,num,niter);
}


void smooth_fs_init (int ndim  /* number of dimensions */, 
			int nbox0  /* triangle radius in frequency */, 
			int **nbox /* triangle radius [ndim-1] */, 
			int *ndat /* data dimensions [ndim] */)
/*< initialize >*/
{
    int i, j, k;

    dim = ndim;
    n = sf_intalloc(dim);
    
    tr = (sf_ctriangle*) sf_alloc(dim*ndat[0],sizeof(sf_ctriangle));
    nd = 1;
    for (i=0; i < dim; i++) {
	s[i] = nd;
	n[i] = ndat[i];
	nd *= ndat[i];
    }
	/*for frequency*/
    tr0 = (nbox0 > 1)? sf_ctriangle_init (nbox0,ndat[0],false): NULL; /*smoother for frequency*/
	/*for x and y*/
    for (i=1; i < dim; i++) {
	    for (j=0; j < nd/n[i]; j++) {
	    k=j%nf;
	    tr[i*nf+k] = (nbox[i-1][k] > 1)? sf_ctriangle_init (nbox[i-1][k],ndat[i],false): NULL;
	    }
    }
    tmp = sf_complexalloc (nd);
}


void smooth_fs_lop (bool adj, bool add, int nx, int ny, sf_complex* x, sf_complex* y)
/*< linear operator >*/
{
    int i, j, k, i0;

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

	/*for frequency*/
	if (NULL != tr0) {
	    for (j=0; j < nd/n[0]; j++) {
		i0 = sf_first_index (0,j,dim,n,s);
		sf_csmooth (tr0, i0, s[0], false, tmp);
	    }
	}

    /*for x and y*/
    for (i=1; i < dim; i++) {
	    for (j=0; j < nd/n[i]; j++) {
	    k=j%nf;
	    if (NULL != tr[i*nf+k]) {
		i0 = sf_first_index (i,j,dim,n,s);
		sf_csmooth (tr[i*nf+k], i0, s[i], false, tmp);
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

void smooth_fs_close (void)
/*< free allocated storage >*/
{
/*	int i;
	for(i=0;i<nf*dim;i++)
	{
		if (NULL != tr[i]) 
			sf_ctriangle_close(tr[i]);
	}*/ 	/*what's wrong with it, only for 2D denoising ???*/
	free(tr);
	if (NULL != tr0) sf_ctriangle_close(tr0);
	free(n);
	free(tmp);
}

