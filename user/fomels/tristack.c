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

#include "tristack.h"

static bool gauss;
static int nx, rect;
static float *t;
static sf_triangle tr;

void tristack_init (bool gauss1 /* pseudo-gaussian */,
		    int ndat /* dense data length */,
		    int nbox /* triangle length */)
/*< initialize >*/
{
    nx = ndat;
    rect = nbox;
    gauss = gauss1;

    t = sf_floatalloc(nx);
    tr = sf_triangle_init (nbox,ndat,false);
}

void  tristack_close(void)
/*< free allocated storage >*/
{
    free(t);
    sf_triangle_close(tr);
}

void tristack (bool adj, bool add, int nc, int nd, float *c, float *d) 
/*< linear operator >*/
{
    int ic, id;

    if (nd != nx) sf_error("%s: wrong size",__FILE__);

    sf_adjnull(adj,add,nc,nd,c,d);

    if (adj) {
	for (id=0; id < nd; id++) {
	    t[id] = d[id];
	}
	if (gauss) {
	    sf_smooth2 (tr, 0, 1, false, t);
	    sf_smooth2 (tr, 0, 1, false, t);
	} else {
	    sf_smooth2 (tr, 0, 1, false, t);
	}
	for (ic=id=0; id < nd; id++) {
	    if (0==id%rect) {
		c[ic] += t[id];
		ic++;
	    } 
	}
    } else {
	for (ic=id=0; id < nd; id++) {
	    if (0==id%rect) {
		t[id] = c[ic];
		ic++;
	    } else {
		t[id] = 0.;
	    }
	}
	if (gauss) {
	    sf_smooth2 (tr, 0, 1, false, t);
	    sf_smooth2 (tr, 0, 1, false, t);
	} else {
	    sf_smooth2 (tr, 0, 1, false, t);
	}
	for (id=0; id < nd; id++) {
	    d[id] += t[id];
	}
    }
}
