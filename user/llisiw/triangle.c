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

#include <rsf.h>
/*^*/

#include "triangle.h"

static void fold (int o, int d, int nx, int nb, 
		  float *x, const float* tmp);
static void doubint (int nx, float *xx, bool der);
static void triple (int o, int d, int nx, int nb, 
		    const float* x, float* tmp, bool box);

static void fold (int o, int d, int nx, int nb, 
		  float *x, const float* tmp)
{
    int i, j;

    /* copy middle */
    for (i=0; i < nx; i++) 
	x[o+i*d] = tmp[i+nb];

    /* reflections from the right side */
    for (j=nb+nx; j < nx+2*nb; j += nx) {
	for (i=0; i < nx && i < nx+2*nb-j; i++)
	    x[o+(nx-1-i)*d] += tmp[j+i];
	j += nx;
	for (i=0; i < nx && i < nx+2*nb-j; i++)
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

static void triple (int o, int d, int nx, int nb, 
		    const float* x, float* tmp, bool box)
{
    int i;
    float wt;

    for (i=0; i < nx + 2*nb; i++) {
	tmp[i] = 0;
    }

    if (box) {
	wt = 1.0/(2*nb-1);

	cblas_saxpy(nx,  +wt,x+o,d,tmp+1   ,1);
	cblas_saxpy(nx,  -wt,x+o,d,tmp+2*nb,1);
    } else {
	wt = 1.0/(nb*nb);
    
	cblas_saxpy(nx,  -wt,x+o,d,tmp     ,1);
	cblas_saxpy(nx,2.*wt,x+o,d,tmp+nb  ,1);
	cblas_saxpy(nx,  -wt,x+o,d,tmp+2*nb,1);
    }
}

void smooth (float *tr       /* smoothing object */, 
	     int o, int d    /* trace sampling */,
	     int nx, int nb  /* data and box length */,
	     bool der        /* if derivative */,
	     bool box        /* if box filter */,
	     float *x        /* data (smoothed in place) */)
/*< apply triangle smoothing >*/
{
    triple (o,d,nx,nb,x,tr,box);
    doubint (nx+2*nb,tr,(bool) (box || der));
    fold (o,d,nx,nb,x,tr);
}

/* 	$Id: triangle.c 7107 2011-04-10 02:04:14Z ivlad $	 */
