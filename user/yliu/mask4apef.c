/* Boundary masks for adaptive PEFs */
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
/*^*/

#include "mask4apef.h"

void mask4apef (int *a                 /* filter size */, 
		int jump               /* dealiasing stretch */, 
		int *nd                /* data size */, 
		float *yy              /* data [nz][ny][nx] */, 
		bool *m                /* first dip mask [nz][ny][nx] */)
/*< data masks for 3-D adaptive PEFs >*/
{
    int i1, i2, i3, j1, j2, j3, n;
    bool *xx;
    
    n = nd[2]*nd[1]*nd[0];
    xx = sf_boolalloc(n);
    
    for (i3=0; i3 < n; i3++) {
	xx[i3] = (bool) (yy[i3] != 0.);
	m[i3] = true;
    }
    
    for (i3=0; i3 < a[2]; i3++) {
	for (i2=0; i2 < a[1]; i2++) {
	    for (i1=-a[0]/2; i1 < (a[0]+1)/2; i1++) {
		for (j3=0; j3 < nd[2]; j3++) {
		    for (j2=0; j2 < nd[1]; j2++) {
			for (j1=0; j1 < nd[0]; j1++) {
			    if ((j1+i1*jump)<0 || (j1+i1*jump)>=nd[0] || 
				(j2+i2*jump)<0 || (j2+i2*jump)>=nd[1] ||
				(j3+i3*jump)<0 || (j3+i3*jump)>=nd[2]) {
				m[j3*nd[1]*nd[0]+j2*nd[0]+j1] = false;
			    } else {
				m[j3*nd[1]*nd[0]+j2*nd[0]+j1] = 
				    (bool) (m[j3*nd[1]*nd[0]+j2*nd[0]+j1] && 
					    xx[(j3+i3*jump)*nd[1]*nd[0]+
					       (j2+i2*jump)*nd[0]+
					       (j1+i1*jump)]);
			    }  
			}
		    }
		}
	    }
	}
    }
    
    free(xx); 
}

/* 	$Id: mask4apef.c 5518 2010-03-15 21:01:37Z yang_liu $	 */

