/* Box smoothing. */
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

#include "box.h"

#include "_bool.h"
/*^*/

#include "alloc.h"
#include "adjnull.h"

static  int np, nb, nx;
static float *tmp;

static void causint (int nx, float *x);
static void causint2 (int nx, float *x);
static void duble (int o, int d, int nx, int nb, 
		    float* x, const float* tmp);
static void duble2 (int o, int d, int nx, int nb, const float* x, float* tmp);

void sf_box_init (int nbox /* box length */, 
		  int ndat /* data length */, 
		  bool lin /* true for linear operator */)
/*< initialize >*/
{
    nx = ndat;
    nb = nbox;
    np = ndat + nbox;
    
    if (lin) tmp = sf_floatalloc(np);
}

void sf_box_close (void)
/*< free allocated storage >*/
{
    free (tmp);
}

static void causint (int nx, float *xx)
/* anti-causal integration */
{
    int i;
    float t;

    /* integrate backward */
    t = 0.;
    for (i=nx-1; i >= 0; i--) {
	t += xx[i];
	xx[i] = t;
    }
}

static void causint2 (int nx, float *xx)
/* causal integration */
{
    int i;
    float t;

    /* integrate forward */
    t=0.;
    for (i=0; i < nx; i++) {
	t += xx[i];
	xx[i] = t;
    }
}

static void duble (int o, int d, int nx, int nb, float* x, const float* tmp)
{
    int i;
    const float *tmp1;
    float wt;

    tmp1 = tmp + nb;
    wt = 1./nb;
    
    for (i=0; i < nx; i++) {
	x[o+i*d] = (tmp[i] - tmp1[i])*wt;
    }
}

static void duble2 (int o, int d, int nx, int nb, const float* x, float* tmp)
{
    int i;
    float *tmp1;
    float wt;

    tmp1 = tmp + nb;
    wt = 1./nb;
    
    for (i=0; i < nx + nb; i++) {
	tmp[i] = 0;
    }

    for (i=0; i < nx; i++) {
	tmp1[i] -= x[o+i*d]*wt;
	tmp[i]  += x[o+i*d]*wt; 
    }
}

void sf_boxsmooth (int o    /* start */, 
		   int d    /* increment */, 
		   float *x /* output */, 
		   float *y /* input */)
/*< adjoint smoothing >*/
{
    causint (np,y);
    duble (o,d,nx,nb,x,y);
}

void sf_boxsmooth2 (int o    /* start */, 
		    int d    /* increment */, 
		    float *x /* input */, 
		    float *y /* output */)
/*< smoothing >*/
{
    duble2 (o,d,nx,nb,x,y);
    causint2 (np,y);
}

void sf_box_lop(bool adj, bool add, int nx, int ny, float* x, float* y) 
/*< smoothing as linear operator >*/
{
    int i;

    sf_adjnull(adj,add,nx,ny,x,y);

    if (adj) {
	sf_boxsmooth(0,1,tmp,y);
	for (i=0; i < nx; i++) {
	    x[i] += tmp[i];
	}
    } else {
	sf_boxsmooth2(0,1,x,tmp);
	for (i=0; i < ny; i++) {
	    y[i] += tmp[i];
	}
    }
}

/* 	$Id: box.c 11806 2014-02-19 00:22:55Z sfomel $	 */

