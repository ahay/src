/* 2-D triangle smoothing as a linear operator */
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

#include "triangle.h"
#include "triangle2.h"

#include <rsf.h>

static int n1, n2, nd;
static triangle tr1, tr2;
static float *tmp;

void triangle2_init (int nbox1, int nbox2 /* triangle size */, 
		     int ndat1, int ndat2 /* data size */)
/*< initialize >*/
{
    n1 = ndat1;
    n2 = ndat2;
    nd = n1*n2;
    tr1 = (nbox1>1)? triangle_init (nbox1,ndat1): NULL;
    tr2 = (nbox2>1)? triangle_init (nbox2,ndat2): NULL;
    tmp = sf_floatalloc (nd);
}

void triangle2_lop (bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{
    int i, i1, i2;

    if (nx != ny || nx != nd) 
	sf_error("%s: Wrong data dimensions: nx=%d, ny=%d, nd=%d",
		 __FILE__,nx,ny,nd);

    sf_adjnull (adj,add,nx,ny,x,y);
  
    if (adj) {
	for (i=0; i < nd; i++) {
	    tmp[i] = y[i];
	}
  
	if (NULL != tr1) {
	    for (i2=0; i2 < n2; i2++) {    
		smooth (tr1, i2*n1, 1, false, tmp);
	    }
	}
  
	if (NULL != tr2) {
	    for (i1=0; i1 < n1; i1++) { 
		smooth (tr2, i1, n1, false, tmp);
	    }
	}

	for (i=0; i < nd; i++) {
	    x[i] += tmp[i];
	}
    } else {
	for (i=0; i < nd; i++) {
	    tmp[i] = x[i];
	}
  
	if (NULL != tr2) {
	    for (i1=0; i1 < n1; i1++) { 
		smooth (tr2, i1, n1, false, tmp);
	    }
	}

	if (NULL != tr1) {
	    for (i2=0; i2 < n2; i2++) {    
		smooth (tr1, i2*n1, 1, false, tmp);
	    }
	}

	for (i=0; i < nd; i++) {
	    y[i] += tmp[i];
	}
    }    
}

void triangle2_close(void)
/*< free allocated storage >*/
{
    free (tmp);
    if (NULL != tr1) triangle_close (tr1);
    if (NULL != tr2) triangle_close (tr2);
}

/* 	$Id$	 */
