/* Operations with complex sinusoids */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "sin.h"

static sf_complex z0;
static int k2, k, nu;
static sf_complex **u;
static float *w, *w1;

void sin_init(sf_complex z1) 
/*< initialize >*/
{
    z0 = conjf(z1);
    
}

void sinsmooth_init(sf_complex z1, 
		    int n2 /* data size */,
		    int k1 /* radius */)
/*< initialize >*/
{
    int ik, i1;
    sf_complex *t1, *t2;

    k = k1;
    k2 = 2*k+1;

    nu = k2*n2;
    u = sf_complexalloc2(k2,n2);
    w = sf_floatalloc(k2);
    w1 = sf_floatalloc(n2);

    for (ik=0; ik < k2; ik++) {
	w[ik]=k+1-SF_ABS(ik-k);
    }

    /* Normalization */
    t1 = sf_complexalloc(n2);
    t2 = sf_complexalloc(n2);

    for (i1=0; i1 < n2; i1++) {
	t1[i1]=sf_cmplx(1.0,0.0);
	w1[i1]=1.0f;
    }

    z0=sf_cmplx(1.0,0.0);
    sin_smooth(false,false,n2,n2,t1,t2);

    z0 = conjf(z1);
    for (i1=0; i1 < n2; i1++) {
	w1[i1]=1.0/cabsf(t2[i1]);
    }

    free(t1);  
    free(t2);
}

void sin_destruct(bool adj, bool add, int nx, int ny, 
		  sf_complex *xx, sf_complex *yy)
/*< destruction operator >*/
{
    int i;

    if (nx != ny) sf_error("%s: wrong dimensions %d != %d",__FILE__,nx,ny);

    sf_ccopy_lop(adj,add,nx,nx,xx,yy);

    for (i = 1; i < nx; i++) {
	if(adj) {
#ifdef SF_HAS_COMPLEX_H	 
	    xx[i-1] -= yy[i] * conjf(z0);
#else
	    xx[i-1] = sf_cadd(xx[i-1],sf_cmul(yy[i],sf_cneg(sf_conjf(z0))));
#endif
	} else {
#ifdef SF_HAS_COMPLEX_H	 
	    yy[i] -= xx[i-1] * z0;
#else
	    yy[i] = sf_cadd(yy[i],sf_cmul(xx[i-1],sf_cneg(z0)));
#endif
	}
    }
}

void sin_construct(bool adj, bool add, int nx, int ny, 
		   sf_complex *xx, sf_complex *yy)
/*< construction operator >*/
{
    int i;
    sf_complex t;

    if (nx != ny) sf_error("%s: wrong dimensions %d != %d",__FILE__,nx,ny);

    sf_cadjnull(adj, add, nx, nx, xx, yy);

    t = sf_cmplx(0.0,0.0);
    if(adj) {	
	for (i = nx-1; i >= 0; i--) {
#ifdef SF_HAS_COMPLEX_H
	    t = yy[i] + t * conjf(z0);
	    xx[i] += t;  
#else
	    t = sf_cadd(yy[i],sf_cmul(t,sf_conjf(z0)));
	    xx[i] = sf_cadd(xx[i],t);
#endif
	}
    } else {
	for (i = 0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H	    
	    t = xx[i] + t * z0;
	    yy[i] += t;
#else
	    t = sf_cadd(xx[i],sf_cmul(t,z0));
	    yy[i] = sf_cadd(yy[i],t);
#endif
	}
    }
}

void sinspray_lop(bool adj, bool add, int n, int nu, sf_complex* u1, sf_complex *u)
/*< spraying operator >*/
{
    int i, ik, ip, j;
    sf_complex t;

    if (nu != n*k2) sf_error("%s: wrong size %d != %d*%d",__FILE__,nu,n,k2);

    sf_cadjnull(adj,add,n,nu,u1,u);

    for (i=0; i < n; i++) { 	
	if (adj) {
	    /* predict forward */
	    t = sf_cmplx(0.0,0.0);
	    for (ik=k-1; ik >= 0; ik--) {
		ip = i+ik+1;
		if (ip >= n) continue;
		j = ip*k2+k+ik+1;
		t += u[j];
		t *= z0;
	    }
	    u1[i] += t;
	    
	    /* predict backward */
	    t = sf_cmplx(0.0,0.0);
	    for (ik=k-1; ik >= 0; ik--) {
		ip = i-ik-1;
		if (ip < 0) continue;
		j = ip*k2+k-ik-1;
		t += u[j];
		t *= conjf(z0);
	    }
	    u1[i] += t;

	    t = u[i*k2+k];
	    u1[i] += t;	  
	} else {
	    t = u1[i];
	    u[i*k2+k] += t;

            /* predict forward */
	    for (ik=0; ik < k; ik++) {
		ip = i-ik-1;
		if (ip < 0) break;
		j = ip*k2+k-ik-1;
		t *= z0;
		u[j] += t;
	    }

	    t = u1[i];
	    
	    /* predict backward */
	    for (ik=0; ik < k; ik++) {
		ip = i+ik+1;
		if (ip >= n) break;
		j = ip*k2+k+ik+1;
		t *= conjf(z0);
		u[j] += t;
	    }
	}
    }
}

void sin_smooth(bool adj, bool add, int n1, int n2, sf_complex* trace, sf_complex *smooth)
/*< smoothing operator >*/
{
    float ws;
    int ik, i1;

    if (n1 != n2) sf_error("%s: wrong size %d != %d",__FILE__,n1,n2);

    sf_cadjnull(adj,add,n1,n2,trace,smooth);

    if (adj) {
	for (i1=0; i1 < n1; i1++) {
	    ws=w1[i1]; 
	    for (ik=0; ik < k2; ik++) {
		u[i1][ik] = smooth[i1]*w[ik]*ws;
	    }
	}

	sinspray_lop(true, add,    n1, nu, trace, u[0]);
    } else {
	sinspray_lop(false, false, n1, nu, trace, u[0]);

	for (i1=0; i1 < n1; i1++) {
	    ws=w1[i1]; 
	    for (ik=0; ik < k2; ik++) {
		smooth[i1] += u[i1][ik]*w[ik]*ws;
	    }
	}
    }
}
