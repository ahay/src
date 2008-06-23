/* ENO interpolation in second dimension + adaptive grid in first dimension */
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

#include "enogrid.h"

#include "grid1.h"
/*^*/

#ifndef _enogrid_h

typedef struct Enogrid *enogrid;
/* abstract data type */
/*^*/

#endif

struct Enogrid {
    int order, ng, n2, nf;
    sf_eno jnt;
    grid1 *grid;
    float **f, *f0;
};
/* concrete data type */

enogrid enogrid_init (int order  /* interpolation order */, 
		      int n2     /* second data dimensions */,
		      int nf     /* number of values */,
		      grid1 *grid /* first dimension grids */)
/*< Initialize interpolation object >*/
{
    enogrid pnt;
    
    pnt = (enogrid) sf_alloc(1,sizeof(*pnt));
    pnt->order = order; 
    pnt->n2 = n2;
    pnt->ng = 2*order-2;
    if (pnt->ng > pnt->n2) sf_error("%s: ng=%d is too big",__FILE__,pnt->ng);
    pnt->jnt = sf_eno_init (order, pnt->ng);
    pnt->nf = nf;
    pnt->f  = sf_floatalloc2(nf,pnt->ng);
    pnt->f0 = sf_floatalloc(pnt->ng);
    pnt->grid = grid;
    
    return pnt;
}

void enogrid_close (enogrid pnt)
/*< Free internal storage >*/
{
    sf_eno_close (pnt->jnt);
    free (pnt->f[0]);
    free (pnt->f);
    free (pnt->f0);
    free (pnt);
}

void enogrid_apply (enogrid pnt, 
		    int j     /* grid location */, 
		    float y /* offset from grid */, 
		    float x /* coordinate on the 1st axis */,
		    float* f  /* output data values */)
/*< Apply interpolation. >*/
{
    int k, b2, i;
    float *f1;
    
    if (j-pnt->order+2 < 0) {
	b2 = 0;
    } else if (j+pnt->order-1 > pnt->n2-1) {
	b2 = pnt->n2 - pnt->ng;
    } else {
	b2 = j-pnt->order+2;
    }
    
    j -= b2;
    
    for (k = 0; k < pnt->ng; k++) {
	grid1_interpolate(pnt->grid[b2+k],x,pnt->nf,pnt->f[k]);
    }
    
    for (i = 0; i < pnt->nf; i++) {
	for (k = 0; k < pnt->ng; k++) {
	    pnt->f0[k] = pnt->f[k][i];
	}

	sf_eno_set (pnt->jnt,pnt->f0);
	sf_eno_apply (pnt->jnt,j,y,f+i,f1,FUNC);
    }
}

/* 	$Id: eno2.c 2086 2006-07-27 23:13:17Z sfomel $	 */
