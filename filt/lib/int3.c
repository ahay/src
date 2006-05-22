/* 3-D interpolation */
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

#include "interp.h"
/*^*/

#include "int3.h"
#include "alloc.h"
#include "adjnull.h"
#include "error.h"
#include "_defs.h"

#include "_bool.h"
/*^*/

static int nd, nf, **nxyz;
static bool *mask, allocated=false;
static   int   m1,   m2,   m3;
static float **w1, **w2, **w3;

void  sf_int3_init (float** coord          /* coordinates [nd][3] */, 
		    float o1, float o2, float o3,
		    float d1, float d2, float d3,
		    int   n1, int   n2,   int n3 /* axes */, 
		    sf_interpolator interp /* interpolation function */, 
		    int nf_in              /* interpolator length */, 
		    int nd_in              /* number of data points */)
/*< initialize >*/
{

    int   id;
    int   i1, i2, i3; 
    float x1, x2, x3, ss;

    nf = nf_in;
    nd = nd_in;
    m1 = n1;
    m2 = n2;
    m3 = n3;
    ss = 1. - 0.5*nf;

    if (!allocated) {
	nxyz = sf_intalloc2  ( 3,nd);
	mask = sf_boolalloc  (   nd);
	w1   = sf_floatalloc2(nf,nd);
	w2   = sf_floatalloc2(nf,nd);
	w3   = sf_floatalloc2(nf,nd);
    }

    for (id = 0; id < nd; id++) {

	x1 = ss + (coord[id][0] - o1)/d1;
	i1 = floorf(x1);
	x1 -= i1;
	
	x2 = ss + (coord[id][1] - o2)/d2;
	i2 = floorf(x2);
	x2 -= i2;

	x3 = ss + (coord[id][2] - o3)/d3;
	i3 = floorf(x3);
	x3 -= i3;
   
	if (i1 > - nf && i1 < n1 &&
	    i2 > - nf && i2 < n2 &&
	    i3 > - nf && i3 < n3) {
	    mask[id] = false; 
	    interp (x1, nf, w1[id]);
	    interp (x2, nf, w2[id]);
	    interp (x3, nf, w3[id]);
	    nxyz[id][0] = i1;
	    nxyz[id][1] = i2;
	    nxyz[id][2] = i3;
	} else {
	    mask[id] = true;
	}
    }
}

void  sf_int3_lop (bool adj, bool add, int nm, int ny, float* mm, float* dd)
/*< linear operator >*/
{ 
    int id, im;
    int i0, j0, k0;
    int i,  j,  k;
    
    if (ny != nd) sf_error("%s: wrong dimensions: %d != %d",__FILE__,ny,nd);
    
    sf_adjnull(adj,add,nm,nd,mm,dd);
    
    for (id=0; id < nd; id++) {

	if (mask[id]) continue;

	i0 = nxyz[id][0]; 
	j0 = nxyz[id][1];
	k0 = nxyz[id][2];

	for         (k = SF_MAX(0,-k0); k < SF_MIN(nf,m3-k0); k++) {
	    for     (j = SF_MAX(0,-j0); j < SF_MIN(nf,m2-j0); j++) {
		for (i = SF_MAX(0,-i0); i < SF_MIN(nf,m1-i0); i++) { 

		    im =(i+i0) + 
			(j+j0)*m1 + 
			(k+k0)*m1*m2;

		    if( adj) { 
			mm[im] += dd[id] * w3[id][k] * w2[id][j] * w1[id][i];
		    } else {
			dd[id] += mm[im] * w3[id][k] * w2[id][j] * w1[id][i];
		    }

		}
	    }
	}
    } /* end id */
}

void int3_close (void)
/*< free allocated storage >*/
{
    if (allocated) {
	allocated = false;
	free (mask);
	free (*nxyz); free (nxyz);
	free (*w1);   free (w1); 
	free (*w2);   free (w2);	
    }
}
