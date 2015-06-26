/* Resampling with triangle weights */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

#include "tristack2.h"

static bool gauss;
static int n1, n2, nc1, rect1, rect2;
static float *t1, *t2;
static sf_triangle tr1, tr2;

void tristack2_init (bool gauss1 /* pseudo-gaussian */,
		     int ndat1, int ndat2 /* dense data length */,
		     int nbox1, int nbox2 /* triangle length */)
/*< initialize >*/
{
    n1 = ndat1;
    n2 = ndat2;

    rect1 = nbox1;
    rect2 = nbox2;

    gauss = gauss1;

    nc1 = (n1-1)/rect1+1;

    t1 = sf_floatalloc(n1*n2);
    t2 = sf_floatalloc(nc1*n2);

    tr1 = sf_triangle_init (rect1,n1,false);
    tr2 = sf_triangle_init (rect2,n2,false);
}

void  tristack2_close(void)
/*< free allocated storage >*/
{
    free(t1);
    free(t2);
    sf_triangle_close(tr1);
    sf_triangle_close(tr2);
}

void tristack2 (bool adj, bool add, int nc, int nd, float *c, float *d) 
/*< linear operator >*/
{
    int ic, id, i1, i2;

    if (nd != n1*n2) sf_error("%s: wrong size",__FILE__);

    sf_adjnull(adj,add,nc,nd,c,d);

    if (adj) {
	for (id=0; id < nd; id++) {
	    t1[id] = d[id];
	}
	for (i2=0; i2 < n2; i2++) {
	    if (gauss) {
		sf_smooth2 (tr1, i2*n1, 1, false, t1);
		sf_smooth2 (tr1, i2*n1, 1, false, t1);
	    } else {
		sf_smooth2 (tr1, i2*n1, 1, false, t1);
	    }

	    for (i1=0; i1 < n1; i1+= rect1) {
		t2[i2*nc1+i1/rect1] = t1[i2*n1+i1];
	    } 
	}
	for (i1=0; i1 < nc1; i1++) {
	    if (gauss) {
		sf_smooth2 (tr2, i1, nc1, false, t2);
		sf_smooth2 (tr2, i1, nc1, false, t2);
	    } else {
		sf_smooth2 (tr2, i1, nc1, false, t2);
	    }
	    
	    for (i2=0; i2 < n2; i2+= rect2) {
		c[i2*nc1/rect2+i1] += t2[i2*nc1+i1];
	    } 
	}
    } else {
	for (i1=0; i1 < nc1; i1++) {
	    for (ic=i2=0; i2 < n2; i2++) {
		if (0==i2%rect2) {
		    t2[i2*nc1+i1] = c[ic*nc1+i1];
		    ic++;
		} else {
		    t2[i2*nc1+i1] = 0.0f;
		}
	    } 

	    if (gauss) {
		sf_smooth2 (tr2, i1, nc1, false, t2);
		sf_smooth2 (tr2, i1, nc1, false, t2);
	    } else {
		sf_smooth2 (tr2, i1, nc1, false, t2);
	    }
	}
	for (i2=0; i2 < n2; i2++) {
	    for (ic=i1=0; i1 < n1; i1++) {
		if (0==i1%rect1) {
		    t1[i2*n1+i1] = t2[i2*nc1+ic];
		    ic++;
		} else {
		    t1[i2*n1+i1] = 0.0f;
		}
	    }
	    if (gauss) {
		sf_smooth2 (tr1, i2*n1, 1, false, t1);
		sf_smooth2 (tr1, i2*n1, 1, false, t1);
	    } else {
		sf_smooth2 (tr1, i2*n1, 1, false, t1);
	    }
	}
	for (id=0; id < nd; id++) {
	    d[id] += t1[id];
	}
    }
}
