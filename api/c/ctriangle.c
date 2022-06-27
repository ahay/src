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

#include "ctriangle.h"
#include "alloc.h"

#include "_bool.h"
#include "komplex.h"
/*^*/

#ifndef _sf_ctriangle_h

typedef struct sf_Ctriangle *sf_ctriangle;
/* abstract data type */
/*^*/

#endif

struct sf_Ctriangle {
    sf_complex *tmp;
    float wt;
    int np, nb, nx;
    bool box;
};

static void fold (int o, int d, int nx, int nb, int np, 
		 sf_complex *x, const sf_complex* tmp);
static void doubint (int nx, sf_complex *x, bool der);
static void triple (int o, int d, int nx, int nb, const sf_complex* x, sf_complex* tmp, bool box, float wt);

sf_ctriangle sf_ctriangle_init (int nbox /* triangle length */, 
				int ndat /* data length */,
				bool box /* if box instead of triangle */)
/*< initialize >*/
{
    sf_ctriangle tr;

    tr = (sf_ctriangle) sf_alloc(1,sizeof(*tr));

    tr->nx = ndat;
    tr->nb = nbox;
    tr->box = box;
    tr->np = ndat + 2*nbox;

    if (box) {
	tr->wt = 1.0/(2*nbox-1);
    } else {
	tr->wt = 1.0/(nbox*nbox);
    }
    
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
    
static void doubint (int nx, sf_complex *xx, bool der)
{
    int i;
    sf_complex t;


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

    if (der) return;

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
}

static void triple (int o, int d, int nx, int nb, const sf_complex* x, sf_complex* tmp, bool box, float wt)
{
    int i;
    sf_complex xi;

    for (i=0; i < nx + 2*nb; i++) {
	tmp[i] = sf_cmplx(0.,0.);
    }

    if (box) {
	for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H
	    xi = wt*x[o+i*d];

	    tmp[i+1]    += xi;
	    tmp[i+2*nb] -= xi;
#else
	    xi = sf_crmul(x[o+i*d],wt);

	    tmp[i+1]    = sf_cadd(tmp[i+1],xi);
	    tmp[i+2*nb] = sf_cadd(tmp[i+2*nb],sf_cneg(xi));
#endif
	}
    } else {
	for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H
	    xi = wt*x[o+i*d];

	    tmp[i]      -=   xi;
	    tmp[i+nb]   += 2*xi;
	    tmp[i+2*nb] -=   xi;
#else
	    xi = sf_crmul(x[o+i*d],wt);

	    tmp[i]      = sf_cadd(tmp[i],sf_cneg(xi));
	    tmp[i+nb]   = sf_cadd(tmp[i+nb],sf_crmul(xi,2.));
	    tmp[i+2*nb] = sf_cadd(tmp[i+2*nb],sf_cneg(xi));
#endif
	}
    }
}

void sf_csmooth (sf_ctriangle tr    /* smoothing object */, 
		 int o, int d    /* trace sampling */, 
		 bool der        /* if derivative */,
		 sf_complex *x   /* data (smoothed in place) */)
/*< apply adjoint triangle smoothing >*/
{
    triple (o,d,tr->nx,tr->nb,x,tr->tmp,tr->box,tr->wt);
    doubint (tr->np,tr->tmp,(bool) (tr->box || der));
    fold (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp);
}

void  sf_ctriangle_close(sf_ctriangle tr)
/*< free allocated storage >*/
{
    free (tr->tmp);
    free (tr);
}

