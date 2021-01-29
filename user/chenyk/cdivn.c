/* Smooth division with several components and complex numbers */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
#include "cdivn.h"

static sf_complex *p, **w;
static int n1, n2, nw;
static sf_coperator oper;

void cmultidivn_init(int nw            /* number of components */, 
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
    
    sf_ctrianglen_init(ndim, nbox, ndat);
    crepeat_init(n,nw,sf_ctrianglen_lop);

    sf_cconjgrad_init(n2, n2, n, n, 1., 1.e-6, verb, false);
    p = sf_complexalloc (n2);
    cweight2_init(nw,n,den);
}

void cmultidivn_close (void)
/*< free allocated storage >*/
{
    sf_ctrianglen_close();
    sf_cconjgrad_close();
    cweight2_close();
    free (p);
}

void cmultidivn (sf_complex* num  /* numerator */, 
		 sf_complex* rat  /* ratio */, 
		 int niter        /* number of iterations */)
/*< smoothly divide num/rat >*/
{
    sf_cconjgrad(NULL,cweight2_lop,crepeat_lop,p,rat,num,niter);
}

void crepeat_init(int m1            /* trace length */, 
		  int m2            /* number of traces */, 
		  sf_coperator oper1 /* operator */)
/*< initialize >*/
{
    n1 = m1;
    n2 = m2;
    oper = oper1;
}

void crepeat_lop (bool adj, bool add, int nx, int ny, sf_complex *xx, sf_complex *yy)
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

void cweight2_init(int nw1        /* number of components */, 
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

void cweight2_close(void)
/*< free allocated storage >*/
{
    free(w);
}

void cweight2_lop (bool adj, bool add, int nx, int ny, sf_complex* xx, sf_complex* yy)
/*< linear operator >*/
{
    int i, iw;

    if (nw*ny != nx) sf_error("%s: size mismatch: %d*%d != %d",
			      __FILE__,nw,ny,nx);

    sf_cadjnull (adj, add, nx, ny, xx, yy);
  
    for (iw=0; iw < nw; iw++) {
	for (i=0; i < ny; i++) {
	    if (adj) {
#ifdef SF_HAS_COMPLEX_H
		xx[i+iw*ny] += yy[i] * conjf(w[iw][i]);
#else
		xx[i+iw*ny] = sf_cadd(xx[i+iw*ny],
				      sf_cmul(yy[i],conjf(w[iw][i])));
#endif
	    } else {
#ifdef SF_HAS_COMPLEX_H
		yy[i] += xx[i+iw*ny] * w[iw][i];
#else
		yy[i] = sf_cadd(yy[i],sf_cmul(xx[i+iw*ny],w[iw][i]));
#endif
	    }
	}
    }
}

/*direct weight*/
static sf_complex* ww;

void cweight_init(sf_complex *w1)
/*< initialize >*/
{
    ww = w1;
}

void cweight_lop (bool adj, bool add, int nx, int ny, 
		  sf_complex* xx, sf_complex* yy)
/*< linear operator >*/
{
    int i;

    if (ny!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,ny,nx);

    sf_cadjnull (adj, add, nx, ny, xx, yy);
  
    for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H
	if (adj) {
	    xx[i] += yy[i] * conjf(ww[i]);
	} else {
	    yy[i] += xx[i] * ww[i];
	}
#else
	if (adj) {
	    xx[i] = sf_cadd(xx[i],sf_cmul(yy[i],conjf(ww[i])));
	} else {
	    yy[i] = sf_cadd(yy[i],sf_cmul(xx[i],ww[i]));
	}
#endif
    }
}
