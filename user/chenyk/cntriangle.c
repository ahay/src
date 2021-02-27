/* Non-stationary triangle smoothing (fixed version) 
Note Fomel's ntriangle.c may contain bugs in user/fomel/ntriangle.c
t[i] -> t[o+i*d]; s[i]->s[o+i*d] in triple */
/*
  Copyright (C) 2020 University of Texas at Austin
  Corrected by Yangkang Chen, Jan, 2020
  
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

#include "cntriangle.h"

#ifndef _cntriangle_h

typedef struct CNtriangle *cntriangle;
/* abstract data type */
/*^*/

#endif

struct CNtriangle {
    sf_complex *tmp;
    int np, nb, nx;
};

static void cfold (int o, int d, int nx, int nb, int np, 
		  const sf_complex *x, sf_complex* tmp);
static void cfold2 (int o, int d, int nx, int nb, int np, 
		   sf_complex *x, const sf_complex* tmp);
static void cdoubint (int nx, sf_complex *x, bool der);
static void ctriple (int o, int d, int nx, int nb, 
		    const float* t, const int* s, sf_complex* x, const sf_complex* tmp);
static void ctriple2 (int o, int d, int nx, int nb, 
		     const float* t, const int* s, const sf_complex* x, sf_complex* tmp);

cntriangle cntriangle_init (int nbox /* maximum triangle length */, 
			  int ndat /* data length */)
/*< initialize >*/
{
    cntriangle tr;

    tr = (cntriangle) sf_alloc(1,sizeof(*tr));

    tr->nx = ndat;
    tr->nb = nbox;
    tr->np = ndat + 2*nbox;
    
    tr->tmp = sf_complexalloc(tr->np);

    return tr;
}

static void cfold (int o, int d, int nx, int nb, int np, 
		  const sf_complex *x, sf_complex* tmp)
{
    int i, j;

    /* copy middle */
    for (i=0; i < nx; i++) 
	tmp[i+nb] = x[o+i*d];
    
    /* reflections from the right side */
    for (j=nb+nx; j < np; j += nx) {
	for (i=0; i < nx && i < np-j; i++)
	    tmp[j+i] = x[o+(nx-1-i)*d];
	j += nx;
	for (i=0; i < nx && i < np-j; i++)
	    tmp[j+i] = x[o+i*d];
    }
    
    /* reflections from the left side */
    for (j=nb; j >= 0; j -= nx) {
	for (i=0; i < nx && i < j; i++)
	    tmp[j-1-i] = x[o+i*d];
	j -= nx;
	for (i=0; i < nx && i < j; i++)
	    tmp[j-1-i] = x[o+(nx-1-i)*d];
    }
}




static void cfold2 (int o, int d, int nx, int nb, int np, 
		   sf_complex *x, const sf_complex* tmp)
{
    int i, j;

    /* copy middle */
    for (i=0; i < nx; i++) 
	x[o+i*d] = tmp[i+nb];

    /* reflections from the right side */
    for (j=nb+nx; j < np; j += nx) {
	for (i=0; i < nx && i < np-j; i++) {
#ifdef SF_HAS_COMPLEX_H
	    x[o+(nx-1-i)*d] += tmp[j+i];
#else
	    x[o+(nx-1-i)*d] = sf_cadd(x[o+(nx-1-i)*d],tmp[j+i]);
#endif
	}
	j += nx;
	for (i=0; i < nx && i < np-j; i++) {
#ifdef SF_HAS_COMPLEX_H
	    x[o+i*d] += tmp[j+i];
#else
	    x[o+i*d] = sf_cadd(x[o+i*d],tmp[j+i]);
#endif
	}
    }
    
    /* reflections from the left side */
    for (j=nb; j >= 0; j -= nx) {
	for (i=0; i < nx && i < j; i++) {
#ifdef SF_HAS_COMPLEX_H
	    x[o+i*d] += tmp[j-1-i];
#else
	    x[o+i*d] = sf_cadd(x[o+i*d],tmp[j-1-i]);
#endif
	}
	j -= nx;
	for (i=0; i < nx && i < j; i++) {
#ifdef SF_HAS_COMPLEX_H
	    x[o+(nx-1-i)*d] += tmp[j-1-i];
#else
	    x[o+(nx-1-i)*d] = sf_cadd(x[o+(nx-1-i)*d],tmp[j-1-i]);
#endif
	}
    }
    
}
    
static void cdoubint (int nx, sf_complex *xx, bool der)
{
    int i;
    sf_complex t;

    /* integrate backward */
    t = sf_cmplx(0.,0.);
    for (i=nx-1; i >= 0; i--) {
#ifdef SF_HAS_COMPLEX_H
	t += xx[i];
#else
	t = sf_cadd(t,xx[i]);
#endif
	xx[i] = t;
    }

    if (der) return;

    /* integrate forward */
    t=sf_cmplx(0.,0.);
    for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H
	t += xx[i];
#else
	t = sf_cadd(t,xx[i]);
#endif
	xx[i] = t;
    }
}

static void ctriple (int o, int d, int nx, int nb, 
		    const float* t,
		    const int* s,
		    sf_complex* x, const sf_complex* tmp)
{
    int i, nt, nt1, ns;
    float tt, wt, wt1;

    for (i=0; i < nx; i++) {
	tt = t[o+i*d];
	nt = floorf(tt);
	nt1 = nt+1;
	ns = nb + s[o+i*d];
	wt  = (nt1*nt1-tt*tt)/(nt*nt*(nt+nt1));
	wt1 = (tt*tt-nt*nt)/(nt1*nt1*(nt+nt1));
	
// 	sf_warning("ix/nx=%d/%d,tt=%g,nt=%d,ns=%d,wt=%g,wt1=%g",i,nx,tt,nt,ns,wt,wt1);
	#ifdef SF_HAS_COMPLEX_H
	x[o+i*d] = 2*(wt+wt1)*tmp[i+ns] - 
	    (tmp[i+ns-nt1] + tmp[i+ns+nt1])*wt1 - 
	    (tmp[i+ns-nt]  + tmp[i+ns+nt ])*wt;
	#else
	x[o+i*d] = 	sf_cadd(sf_cadd(sf_crmul(tmp[i+ns],2*(wt+wt1)), sf_cneg(sf_crmul(sf_cadd(tmp[i+ns-nt1],tmp[i+ns+nt1]),wt1))), sf_cneg(sf_crmul(sf_cadd(tmp[i+ns-nt],tmp[i+ns+nt]),wt)));
	#endif
    }
}

static void ctriple2 (int o, int d, int nx, int nb, 
		     const float* t,
		     const int* s,
		     const sf_complex* x, sf_complex* tmp)
{
    int i, nt, nt1, ns;
	float tt, wt, wt1;

    for (i=0; i < nx + 2*nb; i++) {
	tmp[i] = sf_cmplx(0.,0.);
    }

    for (i=0; i < nx; i++) {
	tt = t[i];
	nt = floorf(tt);
	nt1 = nt+1;
	ns = nb + s[i];
	#ifdef SF_HAS_COMPLEX_H
	wt  = x[o+i*d]*(nt1*nt1-tt*tt)/(nt*nt*(nt+nt1));
	wt1 = x[o+i*d]*(tt*tt-nt*nt)/(nt1*nt1*(nt+nt1));
	tmp[i+ns-nt1] -= wt1; 
	tmp[i+ns-nt]  -= wt; 
	tmp[i+ns]     += 2*(wt+wt1);
	tmp[i+ns+nt]  -= wt;
	tmp[i+ns+nt1] -= wt1;
	#else
	wt=sf_crmul(x[o+i*d],(nt1*nt1-tt*tt)/(nt*nt*(nt+nt1)));
	wt1=sf_crmul(x[o+i*d],(tt*tt-nt*nt)/(nt1*nt1*(nt+nt1));
	tmp[i+ns-nt1] = sf_cadd(tmp[i+ns-nt1],sf_cneg(wt1));
	tmp[i+ns-nt] = sf_cadd(tmp[i+ns-nt],sf_cneg(wt));
	tmp[i+ns] = sf_cadd(tmp[i+ns],sf_crmul(sf_cadd(wt,wt1),2)));
	tmp[i+ns+nt] = sf_cadd(tmp[i+ns+nt],sf_cneg(wt));
	tmp[i+ns+nt1] = sf_cadd(tmp[i+ns+nt1],sf_cneg(wt1));
	#endif
    }
}

void cnsmooth (cntriangle tr /* smoothing object */, 
	      int o, int d /* sampling. o: starting index, d: stride in samples for 1/2/3rd dimension; all refer to a correct index in a 1D vector  */, 
	      bool der     /* derivative flag */, 
	      const float *t /* triangle lengths */, 
	      const int *s /* triangle shifts */,
	      sf_complex *x     /* data (smoothed in place) */)
/*< smooth >*/
{
    cfold (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp); 
    cdoubint (tr->np,tr->tmp,der);
    ctriple (o,d,tr->nx,tr->nb,t,s,x,tr->tmp);
}

void cnsmooth2 (cntriangle tr /* smoothing object */, 
	       int o, int d /* sampling */, 
	       bool der     /* derivative flag */, 
	       const float *t /* triangle lengths */,
	       const int *s /* triangle shifts */,
	       sf_complex *x     /* data (smoothed in place) */)
/*< alternative smooth >*/
{
    ctriple2 (o,d,tr->nx,tr->nb,t,s,x,tr->tmp);
    cdoubint (tr->np,tr->tmp,der);
    cfold2 (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp);
}

void  cntriangle_close(cntriangle tr)
/*< free allocated storage >*/
{
    free (tr->tmp);
    free (tr);
}

