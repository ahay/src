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


#include "triangle1.h"

#include "_bool.h"
/*^*/

#include "triangle.h"
#include "alloc.h"
#include "adjnull.h"
#include "error.h"

static int nd;
static sf_triangle tr;
static float *tmp;

void sf_triangle1_init (int nbox /* triangle size */, 
						int ndat /* data size */)
/*< initialize >*/
{
    nd = ndat;
    tr = sf_triangle_init (nbox,ndat);
    tmp = sf_floatalloc (ndat);
}

void sf_triangle1_lop (bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{
    int i;

    if (nx != ny || nx != nd) sf_error("%s: Wrong data dimensions",__FILE__);

    sf_adjnull (adj,add,nx,ny,x,y);
    if (adj) {
		for (i=0; i < nd; i++) {
			tmp[i] = y[i];
		}
		sf_smooth2 (tr, 0, 1, false, false, tmp);
		for (i=0; i < nd; i++) {
			x[i] += tmp[i];
		}
    } else {
		for (i=0; i < nd; i++) {
			tmp[i] = x[i];
		}
		sf_smooth2 (tr, 0, 1, false, false, tmp);
		for (i=0; i < nd; i++) {
			y[i] += tmp[i];
		}
    }
}

void sf_triangle1_close(void)
/*< free allocated storage >*/
{
    free (tmp);
    sf_triangle_close (tr);
}

/* 	$Id: triangle1.c 8858 2012-07-23 16:33:06Z saragiotis $	 */
