/* Smooth division with several components and complex numbers (stationary version) */
/*
  Copyright (C) 2020 Zhejiang University
  
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
#include "cdivns.h"

static sf_complex *p, **w;
static int n1, n2, nw, nf;
static sf_coperator oper;

void cmultidivns_init(int nw            /* number of components */, 
		       int ndim          /* number of dimensions */, 
		       int n             /* data size */, 
		       int *ndat         /* data dimensions [ndim] */, 
		       int *nbox         /* smoothing radius [ndim] */,
		       sf_complex* den   /* denominator [nw*nd] */,
		       bool verb         /* verbosity flag */)
/*< initialize >*/
{
    int n2;

    n2 = n*nw;
    nf=ndat[0]; /*frequency points*/
    sf_ctrianglen_init(ndim, nbox, ndat);
    crepeats_init(nf,nw,sf_ctrianglen_lop);

    sf_cconjgrad_init(nw*nf,nw*nf, n, n, 1., 1.e-6, verb, false);
    p = sf_complexalloc (n2);
    cweight2s_init(nw,n,den);
}

void cmultidivns_close (void)
/*< free allocated storage >*/
{
    sf_ctrianglen_close();
    sf_cconjgrad_close();
    cweight2s_close();
    free (p);
}

void cmultidivns (sf_complex* num  /* numerator */, 
		 sf_complex* rat  /* ratio */, 
		 int niter        /* number of iterations */)
/*< smoothly divide num/rat >*/
{
    sf_cconjgrad(NULL,cweight2s_lop,crepeats_lop,p,rat,num,niter);
}



void crepeats_init(int m1            /* trace length */, 
		  int m2            /* number of traces */, 
		  sf_coperator oper1 /* operator */)
/*< initialize >*/
{
    n1 = m1;
    n2 = m2;
    oper = oper1;
}

void crepeats_lop (bool adj, bool add, int nx, int ny, sf_complex *xx, sf_complex *yy)
/*< combined linear operator >*/
{
    int i2;       
    
    if (nx != ny || nx != n1*n2) 
	sf_error("%s: Wrong size (nx=%d ny=%d n1=%d n2=%d)",
		 __FILE__,nx,ny,n1,n2);

    sf_cadjnull (adj, add, nx, ny, xx, yy);

    for (i2=0; i2 < n2; i2++) {
	oper(adj,true,n1,n1,xx+i2*n1,yy+i2*n1);
    }
}

void cweight2s_init(int nw1        /* number of components */, 
		   int n          /* model size */, 
		   sf_complex *ww /* weight [nw*n] */)
/*< initialize >*/
{
    int iw;

    nw = nw1;
    w = (sf_complex**) sf_alloc(nw,sizeof(sf_complex*));

    for (iw=0; iw < nw; iw++) {
	w[iw] = ww+iw*n;
    }
}

void cweight2s_close(void)
/*< free allocated storage >*/
{
    free(w);
}

void cweight2s_lop (bool adj, bool add, int nx, int ny, sf_complex* xx, sf_complex* yy)
/*< linear operator >*/
{
    int i,iff, iw;
	int n22; /*n22=n2*n3;ny=n2*n3*nf;nx=nf*nw*/
	n22=ny/nf;


    sf_cadjnull (adj, add, nx, ny, xx, yy);
  
    for (iw=0; iw < nw; iw++) {
    for (iff=0;iff<nf;iff++){
	for (i=0; i < n22; i++) {
	    if (adj) {
#ifdef SF_HAS_COMPLEX_H
		xx[iw*nf+iff] += yy[i*nf+iff] * conjf(w[iw][i*nf+iff]);
#else
		xx[iw*nf+iff] = sf_cadd(xx[iw*nf+iff],
				      sf_cmul(yy[i*nf+iff],conjf(w[iw][i*nf+iff])));
#endif
	    } else {
#ifdef SF_HAS_COMPLEX_H
		yy[i*nf+iff] += xx[iw*nf+iff] * w[iw][i*nf+iff];
#else
		yy[i*nf+iff] = sf_cadd(yy[i*nf+iff],sf_cmul(xx[iw*nf+iff],w[iw][i*nf+iff]));
#endif
	    }
	}
	}
    }
}
