/* 1-D triangle smoothing as a linear operator */
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

#include "triangleh.h"

static int n1,n2,n3;
static sf_triangle tr;
static float *tmp;

void sf_triangleh_init(int nbox /* triangle size */, 
		       int nt /* time number */,
                       int nm /* cmp number */,
                       int nh /* offset number */)
/*< initialize >*/
{
    n3 = nh;
    n2 = nm;
    n1 = nt;
    tr = sf_triangle_init (nbox,n3,false);
    tmp = sf_floatalloc (n3);
}

void sf_triangleh_lop(bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{
    int i1,i2,i3;

    if (nx != ny) sf_error("%s: Wrong data dimensions",__FILE__);

    sf_adjnull (adj,add,nx,ny,x,y);
    if (adj) {
	for (i2=0; i2 < n2; i2++) {
            for (i1=0; i1 < n1; i1++) {
                for (i3=0; i3 < n3; i3++) {
                    tmp[i3] = y[i3*n2*n1+i2*n1+i1];
                }
	        sf_smooth2 (tr, 0, 1, false, tmp);
                for (i3=0; i3 < n3; i3++) {
                    x[i3*n2*n1+i2*n1+i1] += tmp[i3];
                }
             }
	}
    } else {
	for (i2=0; i2 < n2; i2++) {
            for (i1=0; i1 < n1; i1++) {
                for (i3=0; i3 < n3; i3++) {
	            tmp[i3] = x[i3*n2*n1+i2*n1+i1];
                }
	        sf_smooth2 (tr, 0, 1, false, tmp);
                for (i3=0; i3 < n3; i3++) {
                    y[i3*n2*n1+i2*n1+i1] += tmp[i3];
                }
             }
        }
    }
}

void sf_triangleh_close(void)
/*< free allocated storage >*/
{
    free (tmp);
    sf_triangle_close (tr);
}

/* 	$Id: triangleh.c 6676 2012-4-11 08:51:11Z sfomel $	 */
