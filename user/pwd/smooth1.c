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

static int n, nc;
static float *tmp;

void smooth1_init(int n1, int n2     /* data size */,
		  int nc1            /* number of components */,
		  int nb1, int nb2   /* smoothing radius */)
/*< initialize >*/
{
    n = n1*n2;
    nc = nc1;
    sf_triangle2_init (nb1,nb2,n1,n2,1);
    tmp = sf_floatalloc(n);
}

void smooth1_close (void)
/*< free allocated storage >*/
{
    sf_triangle2_close();
    free (tmp);
}

void smooth1_lop (bool adj, bool add, int nx, int ny, float* x, float* y)
/*< smooth each component >*/
{
    int i, ic;

    sf_adjnull (adj,add,nx,ny,x,y);

    for (ic=0; ic < nc; ic++) {
	sf_triangle2_lop (adj,true,n,n,x+ic*n,y+ic*n);
    }

    for (i=nc*n; i < nx; i++) {
	if (adj) {
	    x[i] += y[i];
	} else {
	    y[i] += x[i];
	}
    }
}
    
