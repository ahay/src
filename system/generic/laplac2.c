/* 2-D Laplacian operator */
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
/*^*/

#include "laplac2.h"

static int n1, n2;

void laplac2_init(int m1, int m2)
/*< initialize with data size >*/
{
    n1 = m1;
    n2 = m2;
}

void laplac2_lop(bool adj, bool add, 
		 int np, int nr, float *p, float *r)
/*< linear operator >*/
{
    int i1,i2,j;

    sf_adjnull(adj,add,np,nr,p,r);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    j = i1+i2*n1;
	    if (i1 > 0) {
		if (adj) {
		    p[j-1] -= r[j];
		    p[j]   += r[j];
		} else {
		    r[j] += p[j] - p[j-1];
		}
	    }
	    if (i1 < n1-1) {
		if (adj) {
		    p[j+1] -= r[j];
		    p[j]   += r[j];
		} else {
		    r[j] += p[j] - p[j+1];
		}
	    }
	    if (i2 > 0) {
		if (adj) {
		    p[j-n1] -= r[j];
		    p[j]    += r[j];
		} else {
		    r[j] += p[j] - p[j-n1];
		}
	    }
	    if (i2 < n2-1) {
		if (adj) {
		    p[j+n1] -= r[j];
		    p[j]    += r[j];
		} else {
		    r[j] += p[j] - p[j+n1];
		}
	    }
	}
    }
}

/* 	$Id: laplac2.c 7107 2011-04-10 02:04:14Z ivlad $	 */
