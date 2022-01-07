/* Boundary masks for adaptive PEFs */
/*
  Copyright (C) 2020 University of Texas at Austin
  
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

#include "npef_mask.h"

void mask4apef (int *a                 /* filter size */, 
		int *jumps               /* dealiasing stretch */, 
		int *nd                /* data size */, 
		float *yy              /* data [nz][ny][nx][nhy][nhx] */, 
		bool *m                /* first dip mask [nz][ny][nx][nhy][nhx] */)
/*< data masks for 5-D adaptive PEFs >*/
{
    int i1, i2, i3, i4, i5, j1, j2, j3, j4, j5, n;
    bool *xx;
    
    n = nd[4]*nd[3]*nd[2]*nd[1]*nd[0];
    xx = sf_boolalloc(n);
    
    for (i3=0; i3 < n; i3++) {
	xx[i3] = (bool) (yy[i3] != 0.);
	m[i3] = true;
    }

    for (i5=0; i5 < a[4]; i5++) {
    for (i4=0; i4 < a[3]; i4++) {
    for (i3=0; i3 < a[2]; i3++) {
	for (i2=0; i2 < a[1]; i2++) {
	    for (i1=-a[0]/2; i1 < (a[0]+1)/2; i1++) {
		for (j5=0; j5 < nd[4]; j5++) {
		for (j4=0; j4 < nd[3]; j4++) {	
		for (j3=0; j3 < nd[2]; j3++) {
		    for (j2=0; j2 < nd[1]; j2++) {
			for (j1=0; j1 < nd[0]; j1++) {
			    if ((j1+i1*jumps[0])<0 || (j1+i1*jumps[0])>=nd[0] || 
				(j2+i2*jumps[1])<0 || (j2+i2*jumps[1])>=nd[1] ||
				(j3+i3*jumps[2])<0 || (j3+i3*jumps[2])>=nd[2] ||
				(j4+i4*jumps[3])<0 || (j4+i4*jumps[3])>=nd[3] || 
				(j5+i5*jumps[4])<0 || (j5+i5*jumps[4])>=nd[4]) {
				m[j5*nd[3]*nd[2]*nd[1]*nd[0]+ j4*nd[2]*nd[1]*nd[0]+j3*nd[1]*nd[0]+j2*nd[0]+j1] = false;
			    } else {
				m[j5*nd[3]*nd[2]*nd[1]*nd[0]+ j4*nd[2]*nd[1]*nd[0]+j3*nd[1]*nd[0]+j2*nd[0]+j1] = 
				    (bool) (m[j5*nd[3]*nd[2]*nd[1]*nd[0]+ j4*nd[2]*nd[1]*nd[0]+j3*nd[1]*nd[0]+j2*nd[0]+j1] && 
					    xx[	(j5+i5*jumps[4])*nd[3]*nd[2]*nd[1]*nd[0]+
						(j4+i4*jumps[3])*nd[2]*nd[1]*nd[0]+
						(j3+i3*jumps[2])*nd[1]*nd[0]+
					        (j2+i2*jumps[1])*nd[0]+
					        (j1+i1*jumps[0])]);
			    }  
			}
		    }
		}
	    }
	}
    }
    }
    }
    }
    }
    free(xx); 
}


