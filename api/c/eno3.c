/* ENO interpolation in 3-D */
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

#include "eno2.h"
/*^*/

#include "eno3.h"

#include "alloc.h"
#include "error.h"

#ifndef _sf_eno3_h

typedef struct Eno3 *sf_eno3;
/* abstract data type */
/*^*/

#endif

struct Eno3 {
    int order, ng, n1, n2, n3;
    sf_eno **ent;
    sf_eno2 jnt;
    float **f, **f1;
};
/* concrete data type */

sf_eno3 sf_eno3_init (int order              /* interpolation order */, 
		      int n1, int n2, int n3 /* data dimensions */)
/*< Initialize interpolation object >*/
{
    sf_eno3 pnt;
    int i2, i3;
    
    pnt = (sf_eno3) sf_alloc(1,sizeof(*pnt));
    pnt->order = order; 
    pnt->n1 = n1; 
    pnt->n2 = n2;
    pnt->n3 = n3;
    pnt->ng = 2*order-2;
    if (pnt->ng > n2 || pnt->ng > n3) 
	sf_error("%s: ng=%d is too big",__FILE__,pnt->ng);
    pnt->jnt = sf_eno2_init (order, pnt->ng, pnt->ng);
    pnt->f  = sf_floatalloc2(pnt->ng,pnt->ng);
    pnt->f1 = sf_floatalloc2(pnt->ng,pnt->ng);
    pnt->ent = (sf_eno**) sf_alloc(n3,sizeof(sf_eno*));
    for (i3 = 0; i3 < n3; i3++) {
	pnt->ent[i3] = (sf_eno*) sf_alloc(n2,sizeof(sf_eno));
	for (i2 = 0; i2 < n2; i2++) {
	    pnt->ent[i3][i2] = sf_eno_init (order, n1);
	}
    }

    return pnt;
}

void sf_eno3_set (sf_eno3 pnt, float*** c /* data [n3][n2][n1] */)
/*< Set the interpolation table. c can be changed or freed afterwords. >*/
{
    int i2, i3;
    
    for (i3 = 0; i3 < pnt->n3; i3++) {
	for (i2 = 0; i2 < pnt->n2; i2++) {
	    sf_eno_set (pnt->ent[i3][i2], c[i3][i2]);
	}
    }
}

void sf_eno3_set1 (sf_eno3 pnt, float* c /* data [n3*n2*n1] */)
/*< Set the interpolation table. c can be changed or freed afterwords. >*/
{
    int i2, i3;
    
    for (i3 = 0; i3 < pnt->n3; i3++) {
	for (i2 = 0; i2 < pnt->n2; i2++) {
	    sf_eno_set (pnt->ent[i3][i2], c+(pnt->n1)*(i2+(pnt->n2)*i3));
	}
    }
}

void sf_eno3_close (sf_eno3 pnt)
/*< Free internal storage. >*/
{
    int i2, i3;
    
    sf_eno2_close (pnt->jnt);
    for (i3 = 0; i3 < pnt->n3; i3++) {
	for (i2 = 0; i2 < pnt->n2; i2++) {
	    sf_eno_close (pnt->ent[i3][i2]);
	}
	free (pnt->ent[i3]);
    }
    free (pnt->ent);
    free (pnt->f[0]);
    free (pnt->f);
    free (pnt->f1[0]);
    free (pnt->f1);
    free (pnt);
}

void sf_eno3_apply (sf_eno3 pnt, 
		    int i, int j, int k       /* grid location */, 
		    float x, float y, float z /* offsets from grid */,
		    float* f                  /* output data */, 
		    float* f1                 /* output derivative [3] */, 
		    der what                  /* to compute [FUNC|DER|BOTH] */)
/*< Apply interpolation. >*/
{
    int i2, i3, b2, b3;
    float g;
    
    if (j-pnt->order+2 < 0) {
	b2 = 0;
    } else if (j+pnt->order-1 > pnt->n2-1) {
	b2 = pnt->n2 - pnt->ng;
    } else {
	b2 = j-pnt->order+2;
    }
    
    j -= b2;
    
    
    if (k-pnt->order+2 < 0) {
	b3 = 0;
    } else if (k+pnt->order-1 > pnt->n3-1) {
	b3 = pnt->n3 - pnt->ng;
    } else {
	b3 = k-pnt->order+2;
    }
    
    k -= b3;
    
    for (i3 = 0; i3 < pnt->ng; i3++) {
	for (i2 = 0; i2 < pnt->ng; i2++) {
	    sf_eno_apply (pnt->ent[b3+i3][b2+i2],i,x,
		       &(pnt->f[i3][i2]),
		       &(pnt->f1[i3][i2]),
		       (what==FUNC? FUNC: BOTH));
	}
    }
    
    sf_eno2_set (pnt->jnt,pnt->f);
    sf_eno2_apply (pnt->jnt,j,k,y,z,f,f1+1,what);
    
    if (what != FUNC) {
	sf_eno2_set (pnt->jnt,pnt->f1);
	sf_eno2_apply(pnt->jnt,j,k,y,z,f1,&g,FUNC);
    }
}

/* 	$Id: eno3.c 4148 2009-02-09 03:55:32Z sfomel $	 */
