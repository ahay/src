/* Non-stationary triangle smoothing */
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

#include "ntriangle.h"

#ifndef _ntriangle_h

typedef struct NTriangle *ntriangle;
/* abstract data type */
/*^*/

#endif

struct NTriangle {
    float *tmp;
    int np, nb, nx;
};

static void fold (int o, int d, int nx, int nb, int np, 
		  const float *x, float* tmp);
static void fold2 (int o, int d, int nx, int nb, int np, 
		   float *x, const float* tmp);
static void doubint (int nx, float *x, bool der);
static void triple (int o, int d, int nx, int nb, 
		    const int* t, const int* s, float* x, const float* tmp);
static void triple2 (int o, int d, int nx, int nb, 
		     const int* t, const int* s, const float* x, float* tmp);

ntriangle ntriangle_init (int nbox /* maximum triangle length */, 
			  int ndat /* data length */)
/*< initialize >*/
{
    ntriangle tr;

    tr = (ntriangle) sf_alloc(1,sizeof(*tr));

    tr->nx = ndat;
    tr->nb = nbox;
    tr->np = ndat + 2*nbox;
    
    tr->tmp = sf_floatalloc(tr->np);

    return tr;
}

static void fold (int o, int d, int nx, int nb, int np, 
		  const float *x, float* tmp)
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

static void fold2 (int o, int d, int nx, int nb, int np, 
		   float *x, const float* tmp)
{
    int i, j;

    /* copy middle */
    for (i=0; i < nx; i++) 
	x[o+i*d] = tmp[i+nb];

    /* reflections from the right side */
    for (j=nb+nx; j < np; j += nx) {
	for (i=0; i < nx && i < np-j; i++)
	    x[o+(nx-1-i)*d] += tmp[j+i];
	j += nx;
	for (i=0; i < nx && i < np-j; i++)
	    x[o+i*d] += tmp[j+i];
    }
    
    /* reflections from the left side */
    for (j=nb; j >= 0; j -= nx) {
	for (i=0; i < nx && i < j; i++)
	    x[o+i*d] += tmp[j-1-i];
	j -= nx;
	for (i=0; i < nx && i < j; i++)
	    x[o+(nx-1-i)*d] += tmp[j-1-i];
    }
}
    
static void doubint (int nx, float *xx, bool der)
{
    int i;
    float t;

    /* integrate backward */
    t = 0.;
    for (i=nx-1; i >= 0; i--) {
	t += xx[i];
	xx[i] = t;
    }

    if (der) return;

    /* integrate forward */
    t=0.;
    for (i=0; i < nx; i++) {
	t += xx[i];
	xx[i] = t;
    }
}

static void triple (int o, int d, int nx, int nb, 
		    const int* t,
		    const int* s,
		    float* x, const float* tmp)
{
    int i, nt, ns;
    float wt;

    for (i=0; i < nx; i++) {
	nt = t[i];
	ns = nb + s[i];
	wt = 1./(nt*nt);
	x[o+i*d] = (2.*tmp[i+ns] - tmp[i+ns-nt] - tmp[i+ns+nt])*wt;
    }
}

static void triple2 (int o, int d, int nx, int nb, 
		     const int* t,
		     const int* s,
		     const float* x, float* tmp)
{
    int i, nt, ns;
    float wt;

    for (i=0; i < nx + 2*nb; i++) {
	tmp[i] = 0;
    }

    for (i=0; i < nx; i++) {
	nt = t[i];
	ns = nb + s[i];
	wt = x[o+i*d]/(nt*nt);
	tmp[i+ns-nt] -= wt; 
	tmp[i+ns]     += 2.*wt;
	tmp[i+ns+nt] -= wt;
    }
}

void nsmooth (ntriangle tr /* smoothing object */, 
	      int o, int d /* sampling */, 
	      bool der     /* derivative flag */, 
	      const int *t /* triangle lengths */, 
	      const int *s /* triangle shifts */,
	      float *x     /* data (smoothed in place) */)
/*< smooth >*/
{
    fold (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp); 
    doubint (tr->np,tr->tmp,der);
    triple (o,d,tr->nx,tr->nb,t,s,x,tr->tmp);
}

void nsmooth2 (ntriangle tr /* smoothing object */, 
	       int o, int d /* sampling */, 
	       bool der     /* derivative flag */, 
	       const int *t /* triangle lengths */,
	       const int *s /* triangle shifts */,
	       float *x     /* data (smoothed in place) */)
/*< alternative smooth >*/
{
    triple2 (o,d,tr->nx,tr->nb,t,s,x,tr->tmp);
    doubint (tr->np,tr->tmp,der);
    fold2 (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp);
}

void  ntriangle_close(ntriangle tr)
/*< free allocated storage >*/
{
    free (tr->tmp);
    free (tr);
}

/* 	$Id: ntriangle.c 691 2004-07-04 19:28:08Z fomels $	 */
