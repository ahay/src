/* Triangle smoothing */
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

#include "triangle.h"

#include "_bool.h"
/*^*/

#include "alloc.h"

#ifndef _sf_triangle_h

typedef struct sf_Triangle *sf_triangle;
/* abstract data type */
/*^*/

#endif

struct sf_Triangle {
    float *tmp;
    int np, nb, nx;
};

static void fold (int o, int d, int nx, int nb, int np, 
		  const float *x, float* tmp);
static void fold2 (int o, int d, int nx, int nb, int np, 
		   float *x, const float* tmp);
static void doubint (int nx, float *x, bool der);
static void triple (int o, int d, int nx, int nb, 
		    float* x, const float* tmp);
static void triple2 (int o, int d, int nx, int nb, const float* x, float* tmp);

sf_triangle sf_triangle_init (int nbox /* triangle length */, 
			      int ndat /* data length */)
/*< initialize >*/
{
    sf_triangle tr;

    tr = (sf_triangle) sf_alloc(1,sizeof(*tr));

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

static void triple (int o, int d, int nx, int nb, float* x, const float* tmp)
{
    int i;
    const float *tmp1, *tmp2;
    float wt;

    tmp1 = tmp + nb;
    tmp2 = tmp + 2*nb;
    wt = 1./(nb*nb);
    
    for (i=0; i < nx; i++) {
	x[o+i*d] = (2.*tmp1[i] - tmp[i] - tmp2[i])*wt;
    }
}

static void triple2 (int o, int d, int nx, int nb, const float* x, float* tmp)
{
    int i;
    float *tmp1, *tmp2;
    float wt;

    tmp1 = tmp + nb;
    tmp2 = tmp + 2*nb;
    wt = 1./(nb*nb);
    
    for (i=0; i < nx + 2*nb; i++) {
	tmp[i] = 0;
    }

    for (i=0; i < nx; i++) {
	tmp1[i] += 2.*x[o+i*d]*wt;
	tmp[i]  -= x[o+i*d]*wt; 
	tmp2[i] -= x[o+i*d]*wt;
    }
}

void sf_smooth (sf_triangle tr  /* smoothing object */, 
		int o, int d    /* trace sampling */, 
		bool der        /* if derivative */, 
		float *x        /* data (smoothed in place) */)
/*< apply triangle smoothing >*/
{
    fold (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp);
    doubint (tr->np,tr->tmp,der);
    triple (o,d,tr->nx,tr->nb,x,tr->tmp);
}

void sf_smooth2 (sf_triangle tr  /* smoothing object */, 
		 int o, int d    /* trace sampling */, 
		 bool der        /* if derivative */, 
		 float *x        /* data (smoothed in place) */)
/*< apply adjoint triangle smoothing >*/
{
    triple2 (o,d,tr->nx,tr->nb,x,tr->tmp);
    doubint (tr->np,tr->tmp,der);
    fold2 (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp);
}

void  sf_triangle_close(sf_triangle tr)
/*< free allocated storage >*/
{
    free (tr->tmp);
    free (tr);
}

/* 	$Id$	 */
