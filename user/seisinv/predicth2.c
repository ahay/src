/* Predict two slopes. */
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

#include "predicth.h"

static int n12, n1, n2, n3;
static float **p, **q, *tmp;

void predicth2_init(int m1, int m2, int m3   /* data dimensions */, 
		   float eps                /* regularization parameter */,
		   int order                /* accuracy order */,
		   float** pp, float **qq   /* slopes [m1][m2] */)
/*< initialize >*/
{
    p=pp;
    q=qq;
    n12=m1*m2;
    n1=m1;
    n2=m2;
    n3=m3;

    predicth_init(m1,m2,m3,eps,order,1,false);
    tmp = sf_floatalloc(n12);
}

void predict2_close(void)
/*< free allocated storage >*/
{
    predict_close();
    free(tmp);
}

void predicth2_lop(bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{
    if (nx != n12*n3 || ny != nx) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj,add,nx,ny,x,y);

    if (adj) {
	predict_set(q);
	predicth_lop (true, false, nx, nx, tmp, y);
	predict_set(p);
	predicth_lop (true, add, nx, nx, x, tmp);
    } else {
	predict_set(p);
	predicth_lop (false, false, nx, nx, x, tmp);
	predict_set(q);
	predicth_lop (false, add, nx, nx, tmp, y);
    }
}
