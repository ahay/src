/* interpolation functions */

/*
  Copyright (C) 2011 KAUST
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "interp.h"

#ifndef _interp_h

typedef struct lint1d *lint1dp;
/*^*/

typedef struct lint2d *lint2dp;
/*^*/

struct lint1d{
    float *x; /* known coord */
    float *y; /* target coord */
    int nx,ny;
    int *i;   /* x_i < y < x_{i+1} */
    float *w; /* weighting */
};
/*^*/

struct lint2d{
    int n;
    float *x1,*x2; /* target coord */
    int   *i1,*i2; /* low-low neighbour */
    float *w1,*w2; /* weighting */
};
/*^*/

#endif

lint1dp lint1d_init(float *x, int nx, /* source coord */
		    float *y, int ny  /* target coord */)
/*< initialize 1-D linear interpolation >*/
{
    int iy,j;
    lint1dp L;

    L = (lint1dp)sf_alloc(1,sizeof(*L));
    L->x = sf_floatalloc(nx);
    L->y = sf_floatalloc(ny);
    L->nx = nx;
    L->ny = ny;
    L->i = sf_intalloc(ny);
    L->w = sf_floatalloc(ny);

    /* assuming x0 < x1 < .. and y0 < y1 < .. */
    /* otherwise need sorting */

    /* if y <= x[..], f(y) = f(x[0]) */
    for (iy=0; iy < ny; iy++) {
	if (y[iy] > x[0])
		break;
	else {
	    L->i[iy] = 0;
	    L->w[iy] = 1.0f;
	}
    }

    j = 0;
    while (iy < ny && y[iy] <= x[nx-1]) {
	while (j < nx-1) {
	    if (x[j] < y[iy] && y[iy] <= x[j+1]) break;
	    j++;
	}
	L->i[iy] = j;
	L->w[iy] = (x[j+1] - y[iy])/(x[j+1] - x[j]);
	iy++;
    }

    /* if y > x[..], f(y) = f(x[nx-1]) */
    while (iy < ny) {
	L->i[iy] = nx-2;
	L->w[iy] = 0.0f;
	iy++;
    }

    return L;
}


void lint1d_extract(float *fx, /* source */
		    float *fy, /* target */
		    lint1dp L)
/*< 1D >*/
{
    int j,iy;
    float w;

    for (iy=0; iy < L->ny; iy++) {
	j = L->i[iy];
	w = L->w[iy];
	fy[iy] = w*fx[j] + (1.0f-w)*fx[j+1];
    }
}

void lint2d_info(lint2dp L,int i)
/*< debug >*/
{
    sf_warning("n=%d",L->n);
    sf_warning("x1[%d]=%g,x2[%d]=%g",i,L->x1[i],i,L->x2[i]);
    sf_warning("i1[%d]=%d,i2[%d]=%d",i,L->i1[i],i,L->i2[i]);
    sf_warning("w1[%d]=%g,w2[%d]=%g",i,L->w1[i],i,L->w2[i]);
}

lint2dp lint2d_init(float *x1, /* [n] */
		    float *x2, /* [n] */
		    int n,
		    sf_axis a1 /* axis 1 */,
		    sf_axis a2 /* axis 2 */)
/*< initialize 2D >*/
{
    int n1,n2,k,i;
    float c1,c2,l1,l2,h1,h2,t,o1,o2,d1,d2;

    lint2dp L;
    L = (lint2dp)sf_alloc(1,sizeof(*L));

    L->n  = n;
    L->x1 = sf_floatalloc(n);
    L->x2 = sf_floatalloc(n);
    L->i1 = sf_intalloc(n);
    L->i2 = sf_intalloc(n);
    L->w1 = sf_floatalloc(n);
    L->w2 = sf_floatalloc(n);

    n1 = sf_n(a1); o1 = sf_o(a1); d1 = sf_d(a1);;
    n2 = sf_n(a2); o2 = sf_o(a2); d2 = sf_d(a2);;

    t = o1 + (n1-1)*d1;
    l1 = (d1 >= 0) ? o1 : t;
    h1 = (d1 >= 0) ? t  : o1;

    t = o2 + (n2-1)*d2;
    l2 = (d2 >= 0) ? o2 : t;
    h2 = (d2 >= 0) ? t  : o2;

    for (k=0; k<n; k++) {
	L->x1[k] = c1 = x1[k];
	L->x2[k] = c2 = x2[k];
	
	if (c1 < l1) {
	    L->i1[k] = 0;
	    L->w1[k] = 1;
	} else if (c1 >= h1) {
	    L->i1[k] = n1-2;
	    L->w1[k] = 0;
	} else {
	    L->i1[k] = i = (int)((c1-o1)/d1);
	    L->w1[k] = (c1-(o1+i*d1))/d1;
	}
	
	if (c2 < l2) {
	    L->i2[k] = 0;
	    L->w2[k] = 1;
	} else if (c2 >= h2) {
	    L->i2[k] = n2-2;
	    L->w2[k] = 0;
	} else {
	    L->i2[k] = i = (int)((c2-o2)/d2);
	    L->w2[k] = (c2-(o2+i*d2))/d2;
	}
    }
    return L;
}

void lint2d_extract(float **a,/* [n1*n2] = in */
		    float *b, /* [n] = out */
		    lint2dp L)
/*< extract from neighbours >*/
{
    int k,i1,i2;
    float w1,w2;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) \
    private(k,i1,i2,w1,w2) shared(L,a,b)
#endif
    for (k=0; k < L->n; k++) {
	i1 = L->i1[k];
	i2 = L->i2[k];
	w1 = L->w1[k];
	w2 = L->w2[k];
	b[k]= a[i2  ][i1  ]*(1.0-w1)*(1.0-w2)
	    + a[i2  ][i1+1]*w1*(1.0-w2)
	    + a[i2+1][i1  ]*(1.0-w1)*w2
	    + a[i2+1][i1+1]*w1*w2;
    }
}

void lint2d_inject(float **a,/* [n1*n2] = out */
		   float *b, /* [n] = in */
		   lint2dp L)
/*< inject to neighbours >*/
{
    int k,i1,i2;
    float w1,w2;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) \
    private(k,i1,i2,w1,w2) shared(L,a,b)
#endif
    for (k=0; k < L->n; k++) {
	i1 = L->i1[k];
	i2 = L->i2[k];
	w1 = L->w1[k];
	w2 = L->w2[k];
	a[i2  ][i1  ] += b[k]*(1.0-w1)*(1.0-w2);
	a[i2  ][i1+1] += b[k]*w1*(1.0-w2);
	a[i2+1][i1  ] += b[k]*(1.0-w1)*w2;
	a[i2+1][i1+1] += b[k]*w1*w2;
    }
}

void lint2d_set(float **a,/* [n1*n2] = out */
		float *b, /* [n] = in */
		lint2dp L)
/*< set values to neighbours >*/
{
    int k,i1,i2;
    float w1,w2;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) \
    private(k,i1,i2,w1,w2) shared(L,a,b)
#endif
    for (k=0; k < L->n; k++) {
	i1 = L->i1[k];
	i2 = L->i2[k];
	w1 = L->w1[k];
	w2 = L->w2[k];
	a[i2  ][i1  ] = b[k]*(1.0-w1)*(1.0-w2);
	a[i2  ][i1+1] = b[k]*w1*(1.0-w2);
	a[i2+1][i1  ] = b[k]*(1.0-w1)*w2;
	a[i2+1][i1+1] = b[k]*w1*w2;
    }
}

void lint1_tau(float *x, int nx, float *fx,/* given points x0<x1<..*/
	       float *y, int ny, float *fy /* target points y0<y1<..*/)
/*< 1D linear intepolation >*/
/* used in Mpseudodepth */
{
    float w;
    int iy,ne=0;

    for (iy=0; iy<ny; iy++) { /* y[..]<x[0] */
	if (y[iy]>=x[0]) break;
	else fy[iy]=fx[0];
    }
    while (iy<ny) {
	while (ne < nx-1) {
	    if (x[ne] <= y[iy] && x[ne+1] > y[iy]) break;
	    ne++;
	}
	if (ne==nx-1) {      /* y[..]>x[nx-1] */
	    fy[iy] = fx[nx-1];
	} else {
	    w=(y[iy]-x[ne])/(x[ne+1]-x[ne]);
	    fy[iy] =  (1.0f-w)*fx[ne] + w*fx[ne+1];
	}
	    iy++;
    }

}
