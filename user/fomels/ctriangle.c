/* Triangle smoothing for complex numbers */
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

#include "ctriangle.h"

#ifndef _ctriangle_h

typedef struct cTriangle *ctriangle;
/* abstract data type */
/*^*/

#endif

struct cTriangle {
    sf_complex *tmp;
    int np, nb, nx;
};

static void fold (int o, int d, int nx, int nb, int np, 
		 sf_complex *x, const sf_complex* tmp);
static void doubint (int nx, sf_complex *x, bool der);
static void triple (int o, int d, int nx, int nb, const sf_complex* x, sf_complex* tmp, bool box);

ctriangle ctriangle_init (int nbox /* triangle length */, 
			  int ndat /* data length */)
/*< initialize >*/
{
    ctriangle tr;

    tr = (ctriangle) sf_alloc(1,sizeof(*tr));

    tr->nx = ndat;
    tr->nb = nbox;
    tr->np = ndat + 2*nbox;
    
    tr->tmp = sf_complexalloc(tr->np);

    return tr;
}

static void fold (int o, int d, int nx, int nb, int np, 
		   sf_complex *x, const sf_complex* tmp)
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
    
static void doubint (int nx, sf_complex *xx, bool der)
{
    int i;
    sf_complex t;


    /* integrate forward */
    t=0.;
    for (i=0; i < nx; i++) {
	t += xx[i];
	xx[i] = t;
    }

    if (der) return;

    /* integrate backward */
    t = 0.;
    for (i=nx-1; i >= 0; i--) {
	t += xx[i];
	xx[i] = t;
    }
}

static void triple (int o, int d, int nx, int nb, const sf_complex* x, sf_complex* tmp, bool box)
{
    int i;
    float wt;
    sf_complex xi;

    for (i=0; i < nx + 2*nb; i++) {
	tmp[i] = 0;
    }

    if (box) {
	wt = 1.0/(2*nb-1);
   
	for (i=0; i < nx; i++) {
	    xi = wt*x[o+i*d];

	    tmp[i+1]    += xi;
	    tmp[i+2*nb] -= xi;
	}
    } else {
	wt = 1.0/(nb*nb);
    
	for (i=0; i < nx; i++) {
	    xi = wt*x[o+i*d];

	    tmp[i]      -=   xi;
	    tmp[i+nb]   += 2*xi;
	    tmp[i+2*nb] -=   xi;
	}
    }
}

void csmooth (ctriangle tr    /* smoothing object */, 
	      int o, int d    /* trace sampling */, 
	      bool der        /* if derivative */,
	      bool box        /* if box filter */,
	      sf_complex *x   /* data (smoothed in place) */)
/*< apply adjoint triangle smoothing >*/
{
    triple (o,d,tr->nx,tr->nb,x,tr->tmp,box);
    doubint (tr->np,tr->tmp,(bool) (box || der));
    fold (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp);
}

void  ctriangle_close(ctriangle tr)
/*< free allocated storage >*/
{
    free (tr->tmp);
    free (tr);
}

/* 	$Id: triangle.c 5070 2009-12-04 00:45:45Z sfomel $	 */
