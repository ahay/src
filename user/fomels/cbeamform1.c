/* Beam forming with triangle smoothing: complex data. */
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

#include "ctriangle2.h"

static int n1,n2, rect;
static sf_complex *t;

void cbeamform1_init(int nrep    /* repeat smoothing */,
		    int m1,
		    int m2      /* data dimensions */,
		    int rect1   /* triangle radius */)
/*< initialize >*/  
{
    n1=m1;
    n2=m2;
    rect=rect1;

    t = sf_complexalloc(n1*n2);
    ctriangle2_init(1,rect,n1,n2,nrep);
}

void cbeamform1_close(void)
/*< free allocated storage >*/
{
    free(t);
    ctriangle2_close();
}

void cbeamform1_lop(bool adj, bool add, int nc, int nd, sf_complex* c, sf_complex* d)
/*< linear operator >*/
{
    int i1, ic, id;

    if (nd != n1*n2) sf_error("%s: wrong size",__FILE__);

    sf_cadjnull(adj,add,nc,nd,c,d);

    if (adj) {
	ctriangle2_lop(true,false,nd,nd,t,d);
	for (ic=id=0; id < n2; id++) {
	    if (0==id%rect) {
		for (i1=0; i1 < n1; i1++) {
		    c[ic*n1+i1] += t[id*n1+i1];
		}
		ic++;
	    } 
	}
    } else {
	for (ic=id=0; id < n2; id++) {
	    if (0==id%rect) {
		for (i1=0; i1 < n1; i1++) {
		    t[id*n1+i1] = c[ic*n1+i1];
		}
		ic++;
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    t[id*n1+i1] = sf_cmplx(0.,0.);
		}
	    }
	}
	ctriangle2_lop(false,true,nd,nd,t,d);
    }
}


