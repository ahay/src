/* 2-D interpolation */
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

#include <math.h>

#include <rsf.h>
/*^*/

#include "interp.h"
/*^*/

#include "int2.h"

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif

static int nd, nf, m1, m2, **nxy;
static bool *mask, allocated=false;
static float **w1, **w2;

void  int2_init (float** coord       /* coordinates [nd][2] */, 
		 float o1, float o2, 
		 float d1, float d2,
		 int n1, int n2      /* axes */, 
		 interpolator interp /* interpolation function */, 
		 int nf_in           /* interpolator length */, 
		 int nd_in           /* number of data points */)
/*< initialize >*/
{
    int id, i1, i2; 
    float x1, x2, rx;

    nf = nf_in;
    nd = nd_in;
    m1 = n1;
    m2 = n2;

    if (!allocated) {
	nxy = sf_intalloc2(2,nd);
	mask = sf_boolalloc(nd);
	w1 = sf_floatalloc2(nf,nd);
	w2 = sf_floatalloc2(nf,nd);
    }

    for (id = 0; id < nd; id++) {
	rx = (coord[id][0] - o1)/d1;
	i1 = (int) floorf(rx + 1. - 0.5*nf);
	x1 = rx - floorf(rx);
	
	rx = (coord[id][1] - o2)/d2;
	i2 = (int) floorf(rx + 1. - 0.5*nf);
	x2 = rx - floorf(rx);
   
	if (i1 > - nf && i1 < n1 &&
	    i2 > - nf && i2 < n2) {
	    mask[id] = false; 
	    interp (x1, nf, w1[id]);
	    interp (x2, nf, w2[id]);
	    nxy[id][0] = i1;
	    nxy[id][1] = i2;
	} else {
	    mask[id] = true;
	}
    }
}

void  int2_lop (bool adj, bool add, int nm, int ny, float* x, float* ord)
/*< linear operator >*/
{ 
    int id, i0, j0, i, j, im;
    float w;

    if (ny != nd) sf_error("%s: wrong dimensions: %d != %d",__FILE__,ny,nd);

    sf_adjnull (adj,add,nm,nd,x,ord);

    for (id=0; id < nd; id++) {
	if (mask[id]) continue;
	i0 = nxy[id][0]; 
	j0 = nxy[id][1]; 
	for (j = MAX(0,-j0); j < MIN(nf,m2-j0); j++) {
	    w = w2[id][j];
	    for (i = MAX(0,-i0); i < MIN(nf,m1-i0); i++) { 
		im = (i+i0) + (j+j0)*m1;
		if( adj) { 
		    x[im] += ord[id] * w * w1[id][i];
		} else {
		    ord[id] += x[im] * w * w1[id][i];
		}
	    }
	}
    }
}

void int2_close (void)
/*< free allocated storage >*/
{
    if (allocated) {
	allocated = false;
	free (*nxy);
	free (nxy);
	free (mask);
	free (*w1);
	free (*w2);
	free (w1);
	free (w2);
    }
}

/* 	$Id$	 */
