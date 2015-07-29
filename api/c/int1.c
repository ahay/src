/* 1-D Interpolation */
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

#include "_bool.h"
#include "_defs.h"
#include "alloc.h"
#include "error.h"
#include "adjnull.h"
#include "komplex.h"

#include "interp.h"
/*^*/

#include "int1.h"

static int nd, nf, m1, *nx;
static bool *mask, allocated=false;
static float **w1;

void  sf_int1_init (float* coord               /* cooordinates [nd] */, 
		    float o1, float d1, int n1 /* axis */, 
		    sf_interpolator interp     /* interpolation function */, 
		    int nf_in                  /* interpolator length */, 
		    int nd_in                  /* number of data points */,
		    float tau                  /* interpolation shift */)
/*< initialize >*/
{
    int id, i1; 
    float rx, sx;

    nf = nf_in;
    nd = nd_in;
    m1 = n1;
    sx = 1.0 - 0.5*nf;

    if (!allocated) {
	nx = sf_intalloc(nd);
	mask = sf_boolalloc(nd);
	w1 = sf_floatalloc2(nf,nd);
	allocated = true;
    }

    for (id = 0; id < nd; id++) {
	rx = sx + (coord[id] - o1)/d1 - tau;
	i1 = floorf(rx);
	rx -= i1;
	
	if (i1 > - nf && i1 < n1) {
	    mask[id] = false; 
	    interp (rx, nf, w1[id]);
	    nx[id] = i1;
	} else {
	    mask[id] = true;
	}
    }
}

void  sf_int1_lop (bool adj, bool add, int nm, int ny, float* x, float* ord)
/*< linear operator >*/
{ 
    int id, i0, i, im;
    
    if (ny != nd) sf_error("%s: wrong data size: %d != %d",__FILE__,ny,nd);

    sf_adjnull (adj,add,nm,nd,x,ord);

    for (id=0; id < nd; id++) {
	if (mask[id]) continue;
	i0 = nx[id];

	for (i = SF_MAX(0,-i0); i < SF_MIN(nf,m1-i0); i++) { 
	    im = i+i0;
	    if( adj) { 
		x[im] += ord[id] * w1[id][i];
	    } else {
		ord[id] += x[im] * w1[id][i];
	    }
	}
    }
}

void  sf_cint1_lop (bool adj, bool add, int nm, int ny, sf_complex* x, sf_complex* ord)
/*< linear operator for complex numbers >*/
{ 
    int id, i0, i, im;
    
    if (ny != nd) sf_error("%s: wrong data size: %d != %d",__FILE__,ny,nd);

    sf_cadjnull (adj,add,nm,nd,x,ord);

    for (id=0; id < nd; id++) {
	if (mask[id]) continue;
	i0 = nx[id];

	for (i = SF_MAX(0,-i0); i < SF_MIN(nf,m1-i0); i++) { 
	    im = i+i0;
	    if( adj) { 
#ifdef SF_HAS_COMPLEX_H
		x[im] += ord[id] * w1[id][i];
#else
		x[im] = sf_cadd(x[im],sf_crmul(ord[id],w1[id][i]));
#endif
	    } else {
#ifdef SF_HAS_COMPLEX_H
		ord[id] += x[im] * w1[id][i];
#else
		ord[id] = sf_cadd(ord[id],sf_crmul(x[im],w1[id][i]));
#endif
	    }
	}
    }
}


void sf_int1_close (void)
/*< free allocated storage >*/
{
    if (allocated) {
	free (nx);
	free (mask);
	free (*w1);
	free (w1);
    }
    allocated = false;
}

/* 	$Id: int1.c 13985 2015-03-26 13:56:59Z sfomel $	 */
