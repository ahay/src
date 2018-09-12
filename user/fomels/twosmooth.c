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

static sf_triangle *tr;
static int n, ns;
static float *tmp;

void twosmooth_init(int n1, int n2   /* data size */,
		    int nb1, int nb2 /* smoothing radius */,
		    int skip /* skip dimensions */)
/*< initialize >*/
{
    tr = (sf_triangle*) sf_alloc(2,sizeof(sf_triangle));
    tr[0] = (nb1 > 1)? sf_triangle_init (nb1,n1,false): NULL;
    tr[1] = (nb2 > 1)? sf_triangle_init (nb2,n2,false): NULL;
    n = n1;
    ns = skip;
    tmp = sf_floatalloc(ns+n1+n2);
}

void twosmooth_close (void)
/*< free allocated storage >*/
{
    if (NULL != tr[0]) sf_triangle_close (tr[0]);
    if (NULL != tr[1]) sf_triangle_close (tr[1]);
    free (tr);
    free (tmp);
}

void twosmooth_lop (bool adj, bool add, int nx, int ny, float* x, float* y)
/*< smooth both components >*/
{
    int i;

    sf_adjnull (adj,add,nx,ny,x,y);

    if (adj) {
	for (i=0; i < ny; i++) {
	    tmp[i] = y[i];
	}
    } else {
	for (i=0; i < nx; i++) {
	    tmp[i] = x[i];
	}
    }

    if (NULL != tr[0]) sf_smooth2 (tr[0], 0, 1, false, tmp+ns);
    if (NULL != tr[1]) sf_smooth2 (tr[1], 0, 1, false, tmp+ns+n);
	
    if (adj) {
	for (i=0; i < nx; i++) {
	    x[i] += tmp[i];
	}
    } else {
	for (i=0; i < ny; i++) {
	    y[i] += tmp[i];
	}
    }      
}
