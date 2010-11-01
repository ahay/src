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

#include "upgrad.h"

static upgrad upgreal, upgimag;
static int nt;
static float *temp;

void cpxeikonal_init(int dim   /* number of dimension */,
		     int *n    /* length */,
		     float *d  /* sampling */,
		     int nt1   /* model size */)
/*< initialize operator >*/
{
    upgreal = upgrad_init(dim,n,d);
    upgimag = upgrad_init(dim,n,d);

    nt = nt1;
    temp = sf_floatalloc(nt);
}

void cpxeikonal_set(bool flag  /* real=true / imag=false */,
		    float *t   /* imagtime */)
/*< supply stencil >*/
{
    if (flag)
	upgrad_set(upgreal,t);
    else
	upgrad_set(upgimag,t);
}

void cpxeikonal_pseudo(float *t  /* imag time */,
		       float *w  /* pseudo slowness */)
/*< calculate pseudo slowness from imag time stencil >*/
{
    upgrad_forw(upgimag,t,w);
}

void cpxeikonal_rhs(float *tr   /* real time */,
		    float *ti   /* imag time */,
		    float *rhs  /* right-hand side */)
/*< calculate right-hand side >*/
{
    int it;

    upgrad_forw(upgimag,tr,temp);
    upgrad_forw(upgreal,ti,rhs);

    for (it=0; it < nt; it++) {
	rhs[it] += temp[it];
    }
}

void cpxeikonal_loop(bool adj, bool add, int nm, int nd, float *m, float *d)
/*< linear operator >*/
{
    int it;

    sf_adjnull(adj,add,nm,nd,m,d);

    if (adj) {
	/* adjoint operator: F'd=m */

	upgrad_adj(upgimag,temp,d);
	upgrad_inverse(upgreal,m,temp,NULL);
	upgrad_adj(upgimag,temp,m);
	
	upgrad_adj(upgreal,m,d);
	
	for (it=0; it < nt; it++) {
	    m[it] += temp[it];
	}
    } else {
	/* forward operator: Fm=d */

	upgrad_forw(upgimag,m,temp);
	upgrad_solve(upgreal,temp,d,NULL);
	upgrad_forw(upgimag,d,temp);
	
	upgrad_forw(upgreal,m,d);
	
	for (it=0; it < nt; it++) {
	    d[it] += temp[it];
	}
    }    
}

void cpxeikonal_close()
/*< free memory >*/
{
    upgrad_close(upgreal);
    upgrad_close(upgimag);
    free(temp);
}
