/* 2-D complex triangle smoothing as a linear operator */
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
#include <rsf.h>

static int n1, n2, nd, nr;
static sf_ctriangle tr1, tr2;
static sf_complex *tmp;

void ctriangle2_init (int nbox1, int nbox2 /* triangle size */, 
		      int ndat1, int ndat2 /* data size */,
		      int nrep /* repeat smoothing */)
/*< initialize >*/
{
    n1 = ndat1;
    n2 = ndat2;
    nd = n1*n2;
    nr = nrep;

    tr1 = (nbox1>1)? sf_ctriangle_init (nbox1,ndat1,false): NULL;
    tr2 = (nbox2>1)? sf_ctriangle_init (nbox2,ndat2,false): NULL;
    tmp = sf_complexalloc (nd);
}

void ctriangle2_lop (bool adj, bool add, int nx, int ny, sf_complex* x, sf_complex* y)
/*< linear operator >*/
{
    int i, i1, i2, ir;

    if (nx != ny || nx != nd) 
	sf_error("%s: Wrong data dimensions: nx=%d, ny=%d, nd=%d",
		 __FILE__,nx,ny,nd);

    sf_cadjnull (adj,add,nx,ny,x,y);
  
    if (adj) {
	for (i=0; i < nd; i++) {
	    tmp[i] = y[i];
	}
  
	if (NULL != tr1) {
	    for (i2=0; i2 < n2; i2++) {
		for (ir=0; ir < nr; ir++) {
		    sf_csmooth (tr1, i2*n1, 1, false, tmp);
		}
	    }
	}
  
	if (NULL != tr2) {
	    for (i1=0; i1 < n1; i1++) {
		for (ir=0; ir < nr; ir++) {
		    sf_csmooth (tr2, i1, n1, false, tmp);
		}
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
		for (ir=0; ir < nr; ir++) {
		    sf_csmooth (tr2, i1, n1, false, tmp);
		}
	    }
	}

	if (NULL != tr1) {
	    for (i2=0; i2 < n2; i2++) { 
		for (ir=0; ir < nr; ir++) {
		    sf_csmooth (tr1, i2*n1, 1, false, tmp);
		}
	    }
	}

	for (i=0; i < nd; i++) {
	    y[i] += tmp[i];
	}
    }    
}

void ctriangle2_close(void)
/*< free allocated storage >*/
{
    free (tmp);
    if (NULL != tr1) sf_ctriangle_close (tr1);
    if (NULL != tr2) sf_ctriangle_close (tr2);
}

/* 	$Id: triangle2.c 4092 2009-01-29 21:16:20Z sfomel $	 */
