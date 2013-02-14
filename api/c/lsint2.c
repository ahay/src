/* Least-squares interpolation in 2-D */
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
#include <stdlib.h>

#include "eno.h"
/*^*/

#include "lsint2.h"
#include "alloc.h"

#ifndef _sf_lsint2_h

typedef struct LSint2 *sf_lsint2;
/* abstract data type */
/*^*/

#endif

struct LSint2 {
    int n1, n2;
    float **data2;
};
/* concrete data type */

sf_lsint2 sf_lsint2_init (int n1, int n2 /* data dimensions */)
/*< Initialize interpolation object >*/
{
    sf_lsint2 pnt;

    pnt = (sf_lsint2) sf_alloc(1,sizeof(*pnt));
    pnt->n1 = n1; 
    pnt->n2 = n2;
    pnt->data2 = sf_floatalloc2(n1,n2);
    return pnt;
}

void sf_lsint2_set (sf_lsint2 pnt, float** c /* data [n2][n1] */)
/*< Set the interpolation table. c can be changed or freed afterwords. >*/
{
    int i1, i2;

    for (i2=0; i2 < pnt->n2; i2++) {
	for (i1=0; i1 < pnt->n1; i1++) {
	    pnt->data2[i2][i1] = c[i2][i1];
	}
    }
}

void sf_lsint2_close (sf_lsint2 pnt)
/*< Free internal storage >*/
{
    free(*pnt->data2);
    free(pnt->data2);
    free(pnt);
}

void sf_lsint2_set1 (sf_lsint2 pnt, float* c /* data [n2*n1] */)
/*< Set the interpolation table. c can be changed or freed afterwords. >*/
{
    int i1, i2;

    for (i2=0; i2 < pnt->n2; i2++) {
	for (i1=0; i1 < pnt->n1; i1++) {
	    pnt->data2[i2][i1] = c[i2*pnt->n1+i1];
	}
    }
}

void sf_lsint2_apply (sf_lsint2 pnt, 
		      int i, int j     /* grid location */, 
		      float x, float y /* offset from grid */, 
		      float* f         /* output data value */, 
		      float* f1        /* output derivative [2] */,
		      der what         /* what to compute [FUNC,DER,BOTH] */)
/*< Apply interpolation. >*/
{
    float a, b, c, **d;

    d = pnt->data2;

    if (0==i) {
	if (0==j) {
	    a = (d[j+1][i+1]+d[j+1][i]+
		 d[j  ][i+1]+d[j  ][i])/4.0f;
	    b = (d[j+1][i+1]-d[j+1][i]+
		 d[j  ][i+1]-d[j  ][i])/2.0f;
	    c = (d[j+1][i+1]+d[j+1][i]-
		 d[j  ][i+1]-d[j  ][i])/2.0f;
	} else if (pnt->n2-1==j) {
	    a = (d[j  ][i+1]+d[j  ][i]+
		 d[j-1][i+1]+d[j-1][i])/4.0f;
	    b = (d[j  ][i+1]-d[j  ][i]+
		 d[j-1][i+1]-d[j-1][i])/2.0f;
	    c = (d[j  ][i+1]+d[j  ][i]-
		 d[j-1][i+1]-d[j-1][i])/2.0f;
	} else {
	    a = (d[j+1][i+1]+d[j+1][i]+
		 d[j  ][i+1]+d[j  ][i]+
		 d[j-1][i+1]+d[j-1][i])/6.0f;
	    b = (d[j+1][i+1]-d[j+1][i]+
		 d[j  ][i+1]-d[j  ][i]+
		 d[j-1][i+1]-d[j-1][i])/3.0f;
	    c = (d[j+1][i+1]+d[j+1][i]-
		 d[j-1][i+1]-d[j-1][i])/4.0f;
	}
    } else if (pnt->n1-1==i) {
	if (0==j) {
	    a = (d[j+1][i]+d[j+1][i-1]+
		 d[j  ][i]+d[j  ][i-1])/4.0f;
	    b = (d[j+1][i]-d[j+1][i-1]+
		 d[j  ][i]-d[j  ][i-1])/2.0f;
	    c = (d[j+1][i]+d[j+1][i-1]-
		 d[j  ][i]-d[j  ][i-1])/2.0f;
	} else if (pnt->n2-1==j) {
	    a = (d[j  ][i]+d[j  ][i-1]+
		 d[j-1][i]+d[j-1][i-1])/4.0f;
	    b = (d[j  ][i]-d[j  ][i-1]+
		 d[j-1][i]-d[j-1][i-1])/2.0f;
	    c = (d[j  ][i]+d[j  ][i-1]-
		 d[j-1][i]-d[j-1][i-1])/2.0f;
	} else {
	    a = (d[j+1][i]+d[j+1][i-1]+
		 d[j  ][i]+d[j  ][i-1]+
		 d[j-1][i]+d[j-1][i-1])/6.0f;
	    b = (d[j+1][i]-d[j+1][i-1]+
		 d[j  ][i]-d[j  ][i-1]+
		 d[j-1][i]-d[j-1][i-1])/3.0f;
	    c = (d[j+1][i]+d[j+1][i-1]-
		 d[j-1][i]-d[j-1][i-1])/4.0f;
	} 
    } else {
	if (0==j) {
	    a = (d[j+1][i+1]+d[j+1][i]+d[j+1][i-1]+
		 d[j  ][i+1]+d[j  ][i]+d[j  ][i-1])/6.0f;
	    b = (d[j+1][i+1]-d[j+1][i-1]+
		 d[j  ][i+1]-d[j  ][i-1])/4.0f;
	    c = (d[j+1][i+1]+d[j+1][i]+d[j+1][i-1]-
		 d[j  ][i+1]-d[j  ][i]-d[j  ][i-1])/3.0f;
	} else if (pnt->n2-1==j) {
	    a = (d[j  ][i+1]+d[j  ][i]+d[j  ][i-1]+
		 d[j-1][i+1]+d[j-1][i]+d[j-1][i-1])/6.0f;
	    b = (d[j  ][i+1]-d[j  ][i-1]+
		 d[j-1][i+1]-d[j-1][i-1])/4.0f;
	    c = (d[j  ][i+1]+d[j  ][i]+d[j  ][i-1]-
		 d[j-1][i+1]-d[j-1][i]-d[j-1][i-1])/3.0f;
	} else {
	    a = (d[j+1][i+1]+d[j+1][i]+d[j+1][i-1]+
		 d[j  ][i+1]+d[j  ][i]+d[j  ][i-1]+
		 d[j-1][i+1]+d[j-1][i]+d[j-1][i-1])/9.0f;
	    b = (d[j+1][i+1]-d[j+1][i-1]+
		 d[j  ][i+1]-d[j  ][i-1]+
		 d[j-1][i+1]-d[j-1][i-1])/6.0f;
	    c = (d[j+1][i+1]+d[j+1][i]+d[j+1][i-1]-
		 d[j-1][i+1]-d[j-1][i]-d[j-1][i-1])/6.0f;
	}
    }

    if (DER != what) *f = a+b*x+c*y;
    if (FUNC != what) {
	f1[0] = b;
	f1[1] = c;
    }
}

