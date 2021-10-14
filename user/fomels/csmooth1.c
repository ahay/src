/* Smoothing two components */
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

static sf_ctriangle tr;
static int n, nc;
static sf_complex *tmp;

void csmooth1_init(int n1   /* data size */,
		   int nc1  /* number of components */,
		   int nb   /* smoothing radius */)
/*< initialize >*/
{
    n = n1;
    nc = nc1;
    tr = sf_ctriangle_init (nb,n,false);
    tmp = sf_complexalloc(n);
}

void csmooth1_close (void)
/*< free allocated storage >*/
{
    free (tr);
    free (tmp);
}

void csmooth1_lop (bool adj, bool add, int nx, int ny, sf_complex* x, sf_complex* y)
/*< smooth each component >*/
{
    int i, ic;

    sf_cadjnull (adj,add,nx,ny,x,y);

    for (ic=0; ic < nc; ic++) {
	for (i=0; i < n; i++) {
	    if (adj) {
		tmp[i] = y[i+ic*n];
	    } else {
		tmp[i] = x[i+ic*n];
	    }
	}

	sf_csmooth (tr, 0, 1, false, tmp);

	for (i=0; i < n; i++) {
	    if (adj) {
		x[i+ic*n] += tmp[i];
	    } else {
		y[i+ic*n] += tmp[i];
	    }
	}
    }

    for (i=nc*n; i < nx; i++) {
	if (adj) {
	    x[i] += y[i];
	} else {
	    y[i] += x[i];
	}
    }
}
    
