/* time-squared warping as a linear operator */
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

static int n1, n2, n3;
static sf_map4 mo;

void t2warp_init(int m1 /* input trace length */, 
		 int m2 /* output trace length */,
		 int m3 /* number of traces */,
		 float o1, float d1 /* input sampling */,
		 float o2, float d2 /* output sampling */,
		 float eps /* regularization */)
/*< initialize >*/
{
    int i2;
    float t, *t2;

    n1 = m1;
    n2 = m2;
    n3 = m3;

    mo = sf_stretch4_init (n1, o1, d1, n2, eps);

    t2 = sf_floatalloc(n2);

    for (i2=0; i2 < n2; i2++) {
	t = o2+i2*d2;
	t2[i2] = sqrtf(t);
    } 

    sf_stretch4_define (mo,t2,false);

    free(t2);
}

void t2warp_close(void)
/*< free allocated storage >*/
{
    sf_stretch4_close(mo);
}

void t2warp(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator >*/
{
    int i3;

    if (nx != n3*n1 || ny != n3*n2) sf_error("%s: fwd adj=%d wrong dimensions",__FILE__,adj);

    for (i3=0; i3 < n3; i3++) {
	if (adj) {
	    sf_stretch4_invert_adj (add,mo,y+i3*n2,x+i3*n1);
	} else {
	    sf_stretch4_invert (add,mo,y+i3*n2,x+i3*n1);
	}
    }

}

void t2warp_inv(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator for the inverse warping >*/
{
    int i3;

    if (nx != n3*n2 || ny != n3*n1) sf_error("%s: inv adj=%d wrong dimensions",__FILE__,adj);

    for (i3=0; i3 < n3; i3++) {
	if (adj) {
	    sf_stretch4_apply_adj (add,mo,x+i3*n2,y+i3*n1);
	} else {
	    sf_stretch4_apply (add,mo,x+i3*n2,y+i3*n1);
	}
    } 
}


