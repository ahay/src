/* Diplet transform */
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

#include "diplet.h"
#include "seislet.h"

static int np;
static float ***p;

void diplet_init(int n1      /* trace length */, 
		 int n2      /* number of traces */,
		 int np1     /* number of slopes */,
		 float ***p1  /* slopes */,
		 bool inv    /* inversion flag */, 
		 float eps   /* regularization parameter */,
		 char type   /* transform type */) 
/*< allocate space >*/
{
    seislet_init(n1, n2, inv, true, eps, type);
    np = np1;
    p = p1;
}

void diplet_close(void) 
/*< deallocate space >*/
{
    seislet_close();
}

void diplet_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator >*/
{
    int ip;

    if (nx != ny*np) sf_error("%s: wrong size",__FILE__);

    sf_adjnull (adj,add,nx,ny,x,y);
    for (ip=0; ip < np; ip++) {
	seislet_set(p[ip]);
	seislet_lop(adj, true, ny, ny, x+ip*ny, y);
    }
}
