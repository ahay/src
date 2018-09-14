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

static sf_triangle *tr1, *tr2;
static int n, ns, ny1, ny2, nx1, nx2;
static float *tmp;

void twosmooth2_init(int n1, int n2   /* data size */,
		     int nt1, int nt2 /* trace length */,
		     int nb1, int nb2 /* smoothing radius */,
		     int nc1, int nc2 /* smoothing radius */,
		     int skip /* skip dimensions */)
/*< initialize >*/
{
    tr1 = (sf_triangle*) sf_alloc(2,sizeof(sf_triangle));
    tr2 = (sf_triangle*) sf_alloc(2,sizeof(sf_triangle));

    ny1 = nt1; nx1=n1/nt1;
    ny2 = nt2; nx2=n2/nt2;
    
    tr1[0] = (nb1 > 1)? sf_triangle_init (nb1,ny1,false): NULL;
    tr1[1] = (nb2 > 1)? sf_triangle_init (nb2,nx1,false): NULL;
    tr2[0] = (nc1 > 1)? sf_triangle_init (nc1,ny2,false): NULL;
    tr2[1] = (nc2 > 1)? sf_triangle_init (nc2,nx2,false): NULL;
    
    n = n1;
    ns = skip;
    tmp = sf_floatalloc(ns+n1+n2);
}

void twosmooth2_close (void)
/*< free allocated storage >*/
{
    if (NULL != tr1[0]) sf_triangle_close (tr1[0]);
    if (NULL != tr1[1]) sf_triangle_close (tr1[1]);
    if (NULL != tr2[0]) sf_triangle_close (tr2[0]);
    if (NULL != tr2[1]) sf_triangle_close (tr2[1]);
    free (tr1);
    free (tr2);
    free (tmp);
}

void twosmooth2_lop (bool adj, bool add, int nx, int ny, float* x, float* y)
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

    if (NULL != tr1[0]) {
	for (i=0; i < nx1; i++) {
	    sf_smooth2 (tr1[0], i*ny1, 1, false, tmp+ns);
	}
    }
    if (NULL != tr1[1]) {
	for (i=0; i < ny1; i++) {
	    sf_smooth2 (tr1[1], i, ny1, false, tmp+ns);
	}
    }

    if (NULL != tr2[0]) {
	for (i=0; i < nx2; i++) {
	    sf_smooth2 (tr2[0], i*ny2, 1, false, tmp+ns+n);
	}
    }
    if (NULL != tr2[1]) {
	for (i=0; i < ny2; i++) {
	    sf_smooth2 (tr2[1], i, ny2, false, tmp+ns+n);
	}
    }
	
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
