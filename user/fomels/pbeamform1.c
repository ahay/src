/* 1-D beam forming. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "pbeamform1.h"

static bool gauss;
static int n1, rect;
static float *t, *t2;

void pbeamform1_init(bool gauss1    /* pseudo-gaussian */,
		     int m1         /* data dimension */, 
		     int rect1      /* triangle radius */)
/*< initialize >*/
{
    n1=m1;
    rect=rect1;
    gauss=gauss1;

    t = sf_floatalloc(n1);
    t2 = gauss? sf_floatalloc(n1): NULL;

    sf_triangle1_init(rect,n1);
}

void pbeamform1_close(void)
/*< free allocated storage >*/
{
    free(t);
    if (gauss) free(t2);
    sf_triangle1_close();
}

void pbeamform1_lop(bool adj, bool add, int nc, int nd, float* c, float* d)
/*< linear operator >*/
{
    int ic, id;

    if (nd != n1) sf_error("%s: wrong size",__FILE__);

    sf_adjnull(adj,add,nc,nd,c,d);

    if (adj) {
	if (gauss) {
	    sf_triangle1_lop(false,false,nd,nd,d,t2);
	    sf_triangle1_lop(true,false,nd,nd,t,t2);
	} else {
	    sf_triangle1_lop(true,false,nd,nd,t,d);
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
	    sf_triangle1_lop(false,false,nd,nd,t,t2);
	    sf_triangle1_lop(true,true,nd,nd,d,t2);
	} else {
	    sf_triangle1_lop(false,true,nd,nd,t,d);
	}
    }
}


