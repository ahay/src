/* Nonstationary smoothing */
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

#include "ntriangle2.h"

static int n0, n1, n2;
static sf_triangle *tr;
static float *tmp;

void ntrianglen2_init(int nw      /* number of components */, 
		      int n       /* data size */,
		      int nf      /* first dimension */,
		      int *nbox   /* smoothing radius [nw] */)
/*< initialization >*/
{
    int i2;

    n0 = nf;
    n1 = n;
    n2 = nw;

    tmp = sf_floatalloc(n1*n2);

    tr = (sf_triangle *) sf_alloc(nw,sizeof(*tr));

    for (i2=0; i2 < n2; i2++) {
	tr[i2] = sf_triangle_init (nbox[i2],n0,false);
    }
}

void ntriangle2_close(void) 
/*< free allocated storage >*/
{
    int i2;

    for (i2=0; i2 < n2; i2++) {
		sf_triangle_close(tr[i2]);
    }

    free(tr);
    free(tmp);
}

void ntriangle2_lop (bool adj, bool add, int nx, int ny, float *x, float *y)
/*< combined linear operator >*/
{
    int i, j, i2;       
    
    if (nx != ny || nx != n1*n2) 
	sf_error("%s: Wrong size (nx=%d ny=%d n1=%d n2=%d)",
		 __FILE__,nx,ny,n1,n2);

    sf_adjnull (adj, add, nx, ny, x, y);

    if (adj) {
		for (i=0; i < nx; i++) {
			tmp[i] = y[i];
		}
    } else {
		for (i=0; i < nx; i++) {
			tmp[i] = x[i];
		}
    }

    for (i2=0; i2 < n2; i2++) {
		for (j=0; j < n1/n0; j++) {
			sf_smooth2 (tr[i2],j*n0,1,false,tmp+i2*n1);
		}
    }

    if (adj) {
    	for (i=0; i < nx; i++) {
			x[i] += tmp[i];
		}
    } else {
		for (i=0; i < nx; i++) {
			y[i] += tmp[i];
		}
    }
}
