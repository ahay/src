/* 2-D adaptive grid */
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

#include "agrid2.h"
#include "agrid.h"

#ifndef _agrid2_h

typedef struct AGrid2 *agrid2;
/* abstract data type */
/*^*/

#endif

/* concrete data type */
struct AGrid2 {
    int order, ng, n1, n2, nd;
    sf_eno jnt;
    agrid *ent;
    float **f, *f1;
};

agrid2 agrid2_init (int order      /* interpolation order */, 
		    int n1, int n2 /* data dimensions */, 
		    int nd, int max)
/*< Initialize interpolation object >*/
{
    agrid2 pnt;
    int i2;

    pnt = (agrid2) sf_alloc(1,sizeof(*pnt));
    pnt->order = order; 
    pnt->n1 = n1; 
    pnt->n2 = n2;
    pnt->nd = nd;
    pnt->ng = 2*order-2;
    if (pnt->ng > pnt->n2) sf_error("ng is too big in agrid2");
    pnt->jnt = sf_eno_init (order, pnt->ng);
    pnt->f  = sf_floatalloc2(pnt->nd,pnt->ng);
    pnt->f1 = sf_floatalloc(pnt->ng);
    pnt->ent = (agrid*) sf_alloc(n2,sizeof(agrid));
    for (i2 = 0; i2 < n2; i2++) {
	pnt->ent[i2] = agrid_init (n1, nd, max);
    }

    return pnt;
}

void agrid2_set (agrid2 pnt, float*** c /* [n2][n1] */)
/*< Set the interpolation table >*/
{
    int i2;

    for (i2 = 0; i2 < pnt->n2; i2++) {
	agrid_set (pnt->ent[i2], c[i2]);
    }
}

void agrid2_set1 (agrid2 pnt, float** c /* [n2*n1] */)
/*< Set the interpolation table >*/
{
    int i2;

    for (i2 = 0; i2 < pnt->n2; i2++) {
	agrid_set (pnt->ent[i2], c+i2*(pnt->n1));
    }
}

void agrid2_close (agrid2 pnt)
/*< free internal storage >*/
{
    int i2;

    sf_eno_close (pnt->jnt);
    for (i2 = 0; i2 < pnt->n2; i2++) {
	agrid_close (pnt->ent[i2]);
    }
    free (pnt->f[0]);
    free (pnt->f);
    free (pnt->f1);
    free (pnt->ent);
    free (pnt);
}

void agrid2_interp (agrid2 pnt, int i, int j, float x, float y, float* f)
/*< Interpolation >*/
{
    int id, k, b2;
    float g;
    
    if (j-pnt->order+2 < 0) {
	b2 = 0;
    } else if (j+pnt->order-1 > pnt->n2-1) {
	b2 = pnt->n2 - pnt->ng;
    } else {
	b2 = j-pnt->order+2;
    }
    
    j -= b2;
    
    for (k = 0; k < pnt->ng; k++) {
	agrid_interp (pnt->ent[b2+k],i,x,pnt->f[k]);
    }

    for (id=0; id < pnt->nd; id++) {
	for (k = 0; k < pnt->ng; k++) {
	    pnt->f1[k] = pnt->f[k][id];
	}
	sf_eno_set (pnt->jnt,pnt->f1);
	sf_eno_apply (pnt->jnt,j,y,f+id,&g,FUNC);
    }
}

void agrid2_fill (agrid2 pnt, float min1, float max1, float min2, float max2,
		  void (*fill)(float x, void* dat, float* f))
/*< Populate the grid. x is relative to the start of the grid >*/
{
    int i2;

    for (i2 = 0; i2 < pnt->n2; i2++) {
	fill_grid (pnt->ent[i2],min1,max1,min2,max2,&i2,fill);
    }
}

void agrid2_size (agrid2 pnt, float* size)
/*< determine grid size >*/
{
    int i2;

    for (i2 = 0; i2 < pnt->n2; i2++) {
	size[i2] = (float) grid_size (pnt->ent[i2]);
    }
}

void agrid2_write (agrid2 pnt, float*** dat) 
/*< dump the grid in an array >*/
{
    int i2;

    for (i2 = 0; i2 < pnt->n2; i2++) {
	dat[i2] = write_grid (pnt->ent[i2]);
    }
}
