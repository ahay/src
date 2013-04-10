/* Apply seislet transform on each component. */
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

#include "seislet.h"

static int n, n12;
static float ***p;

void seisletk_init(int nk      /* number of components */,
		   int n1      /* trace length */, 
		   int n2      /* number of traces */, 
		   bool inv    /* inversion flag */, 
		   float eps   /* regularization parameter */,
		   int order   /* accuracy order */,
		   char type   /* transform type */,
		   float*** pk /* slopes [nk][m1][m2] */)
/*< initialize >*/
{
    n=nk;
    p=pk;
    n12=n1*n2;

    seislet_init(n1,n2,inv,false,eps,order,type);
}

void seisletk_close(void)
/*< free allocated storage >*/
{
    seislet_close();
}

void seisletk_lop(bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{
    int i;

    if (nx != n*n12 || ny != nx) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj,add,nx,ny,x,y);

    for (i=0; i < n; i++) { /** loop over components **/
	seislet_set(p[i]);
	seislet_lop(adj,true,n12,n12,x+i*n12,y+i*n12);
    }
}



