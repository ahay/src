/* Linear operator for linearized complex eikonal equation  */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
#include "upgradcpx.h"

static upgrad upgreal, upgimag;
static float *temp1, *temp2, *temp;

void cpxeiko_init(int dim,
		  int *n,
		  int nm, /* model size */
		  float *d,
		  float *tr  /* realtime */,
		  float *ti  /* imagtime */)
/*< initialize operator >*/
{
    temp1 = sf_floatalloc(nm);
    temp2 = sf_floatalloc(nm);
    temp  = sf_floatalloc(nm);

    upgreal = upgrad_init(dim,n,d);
    upgimag = upgrad_init(dim,n,d);

    upgrad_set(upgreal,tr);
    upgrad_set(upgimag,ti);
}

void cpxeiko_operator(bool adj, bool add, int nm, int nd, float *m, float *d)
/*< linear operator >*/
{
    int it;

    if (nm != nd)
	sf_error("dim(model) != dim(data)");

    sf_adjnull(adj,add,nm,nd,m,d);

    if (adj) {
	upgrad_adj(upgreal,temp,d);
	upgrad_inverse(upgimag,temp1,temp,NULL);

	upgrad_adj(upgimag,temp,d);
	upgrad_inverse(upgreal,temp2,temp,NULL);

	for (it=0; it < nm; it++) {
	    m[it] += temp1[it] + temp2[it]; 
	}
	/* adjoint operator: F'd=m */
    } else {
	upgrad_solve(upgimag,m,temp,NULL);
	upgrad_forw(upgreal,temp,temp1);

	upgrad_solve(upgreal,m,temp,NULL);
	upgrad_forw(upgimag,temp,temp2);

	for (it=0; it < nd; it++) {
	    d[it] += temp1[it] + temp2[it]; 
	}
	/* forward operator: Fm=d */
    }

}

void cpxeiko_forw(bool flag   /* real=true / imag=false */,
		  float *in   /* data */,
		  float *rhs  /* right-hand side */)
/*< calculate right-hand side >*/
{
    if (flag) {
	upgrad_forw(upgreal,in,rhs);
    } else {
	upgrad_forw(upgimag,in,rhs);
    }
}

void cpxeiko_close()
/*< free allocated memory >*/
{
    free(upgreal);
    free(upgimag);
    free(temp1);
    free(temp2);
    free(temp);
}
