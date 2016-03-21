/* Nonstationary triangle smoothing as a linear operator, applied in 2-D */
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

#include "ntriangle1.h"
#include "ntriangle.h"

#include <rsf.h>
/*^*/

static int n1, n2, nd, **ns;
static ntriangle tr;
static float *tmp, **nr;

void ntriangle1_init (int nbox1            /* maximum triangle size */, 
		      int ndat1, int ndat2 /* data size */,
		      float **rect           /* triangle sizes */,
                      int **shift          /* triangle shifts */)
/*< initialize >*/
{
    n1 = ndat1;
    n2 = ndat2;
    nd = n1*n2;
    nr = rect;
    ns = shift;

    tr = ntriangle_init (nbox1,ndat1);
    tmp = sf_floatalloc (n1);
}

void ntriangle1_lop (bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{
    int i1, i2;

    if (nx != ny || nx != nd) 
	sf_error("%s: Wrong data dimensions: nx=%d, ny=%d, nd=%d",
		 __FILE__,nx,ny,nd);

    sf_adjnull (adj,add,nx,ny,x,y);
  
    for (i2=0; i2 < n2; i2++) {
	if (adj) {
	    for (i1=0; i1 < n1; i1++) {
		tmp[i1] = y[i1+i2*n1];
	    }
  
	    nsmooth2 (tr, 0, 1, false, nr[i2], ns[i2], tmp);

	    for (i1=0; i1 < n1; i1++) {
		x[i1+i2*n1] += tmp[i1];
	    }
	} else {
	    for (i1=0; i1 < n1; i1++) {
		tmp[i1] = x[i1+i2*n1];
	    }
  
	    nsmooth (tr, 0, 1, false, nr[i2], ns[i2], tmp);

	    for (i1=0; i1 < n1; i1++) {
		y[i1+i2*n1] += tmp[i1];
	    }
	}
    }    
}

void ntriangle1_close(void)
/*< free allocated storage >*/
{
    free (tmp);
    ntriangle_close (tr);
}

/* 	$Id: ntriangle1.c 839 2004-10-25 11:54:43Z fomels $	 */
