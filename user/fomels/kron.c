/* Kroneker product of square matrices */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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

static int n;
static float **a, **b;

void kron_init(int n1, float **a1, float**b1)
/*< initialize >*/
{
    n = n1;
    a = a1;
    b = b1;
}

void kron_lop(bool adj, bool add, int nf, int ng, float *f, float *g)
/*< linear operator: G = B * F * A' >*/
{
    int i, j, ii, k, l, kk;

    if (nf != ng || nf != n*n) sf_error("%s: wrong size",__FILE__);

    sf_adjnull(adj,add,nf,ng,f,g);

    for (kk=k=0; k < n; k++) {
	for (l=0; l < n; l++, kk++) {
	    for (ii=i=0; i < n; i++) {
		for (j=0; j < n; j++, ii++) {
		    if (adj) {
			f[kk] += b[i][k]*a[j][l]*g[ii];
		    } else {
			g[ii] += b[i][k]*a[j][l]*f[kk];
		    }
		}
	    }
	}
    }
}
