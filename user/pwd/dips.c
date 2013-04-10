/* Estimating constant dips */
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

#include "dips.h"

#include "callp2.h"

static int nd, n1, n2, n12, skip;
static callpass2 *ap;
static float *x, **tmp1, **tmp2;

void dips_init(int nd1        /* number of dips */, 
	       int nw         /* accuracy order */, 
	       int nj         /* filter stretch for aliasing */, 
	       int nx, int ny /* data dimensions */, 
	       float** x1     /* data [ny][nx] */)
/*< initialize >*/
{
    int id;

    x = x1[0];
    nd = nd1;
    ap = (callpass2 *) sf_alloc(nd,sizeof(callpass2));
    for (id=0; id < nd; id++) {
	ap[id] = callpass2_init(nw, nj, nx, ny);
    }
    tmp1 = sf_floatalloc2(nx,ny);
    tmp2 = sf_floatalloc2(nx,ny);

    n1 = nx;
    n2 = ny;
    n12 = nx*ny;

    skip = nd*nw*nj;
}

void dips_close(void)
/*< free allocated storage >*/
{
    int id;

    for (id=0; id < nd; id++) {
	free (ap[id]);
    }
    free (ap);
    free (tmp1[0]);
    free (tmp1);
    free (tmp2[0]);
    free (tmp2);
}

void dips(const float *d /* initial dip [nd] */, 
	  float *b       /* right-hand side [nx*ny] */, 
	  float **aa     /* matrix to invert [nx*ny][nd] */)
/*< estimate dips >*/
{
    int id, jd, i1, i2, i;
    float **tmp;

    for (i=0; i < n12; i++) {
	tmp2[0][i] = x[i];
    }

    for (id=0; id < nd; id++) {
	tmp = tmp1; tmp1 = tmp2; tmp2 = tmp;

	callpass21_set (ap[id], d[id]);
	callpass21 (false, ap[id], tmp1, tmp2);
    }

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    i = i2*n1+i1;
	    if (i2 < n2-nd && i1 >= skip && i1 < n1-skip) {
		b[i] = tmp2[i2][i1];
	    } else {
		b[i] = 0;
	    }
	}
    }

    for (id=0; id < nd; id++) {
	for (i=0; i < n12; i++) {
	    tmp2[0][i] = x[i];
	}

	for (jd=0; jd < nd; jd++) {
	    tmp = tmp1; tmp1 = tmp2; tmp2 = tmp;
	    callpass21 ((bool) (jd == id), ap[jd], tmp1, tmp2);
	}
	
	
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		i = i2*n1+i1;
		if (i2 < n2-nd && i1 >= skip && i1 < n1-skip) {
		    aa[i][id] = tmp2[i2][i1];
		} else {
		    aa[i][id] = 0.;
		}
	    }
	}
    }
}


