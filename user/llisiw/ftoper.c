/* First-arrival Traveltime Tomography (linear operator) */
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

#include "ftoper.h"

static int *n, ns, **m;
static float *d, *tempr, *tempx;
static sf_upgrad *upglist;

void ftoper_init(int *nn /* dimension */,
		 float *dd /* sampling */,
		 int nshot /* number of shots */,
		 int **mm /* mask */)
/*< initialize >*/
{
    int is;

    n = nn;
    d = dd;
    ns = nshot;
    m = mm;

    upglist = (sf_upgrad *) sf_alloc(ns,sizeof(sf_upgrad));

    for (is=0; is < ns; is++) {
	upglist[is] = sf_upgrad_init(2,n,d);
    }

    tempr = sf_floatalloc(n[0]*n[1]*n[2]);
    tempx = sf_floatalloc(n[0]*n[1]*n[2]);
}

void ftoper_set(float **t /* time */)
/*< set-up >*/
{
    int is;

    for (is=0; is < ns; is++) {
	sf_upgrad_set(upglist[is],t[is]);
    }
}

void ftoper_close(void)
/*< close >*/
{
    int is;

    for (is=0; is < ns; is++) {
	sf_upgrad_close(upglist[is]);
    }

    free(upglist);
}

void ftoper_oper(bool adj, bool add, int nx, int nr, float *x, float *r)
/*< linear operator >*/
{
    int i, k, is;

    sf_adjnull(adj,add,nx,nr,x,r);

    if (adj) {
	for (is=0; is < ns; is++) {
	    for (k=0; k < n[1]*n[2]; k++) {
		tempr[k*n[0]] = (m[is][k]==1)? r[is*n[0]*n[1]*n[2]+k*n[0]]: 0.;

		for (i=1; i < n[0]; i++) {
		    tempr[k*n[0]+i] = 0.;
		}
	    }

	    sf_upgrad_inverse(upglist[is],tempx,tempr,NULL);

	    for (i=0; i < n[0]*n[1]*n[2]; i++) {
		x[i] += tempx[i];
	    }
	}
    } else {
	for (is=0; is < ns; is++) {
	    sf_upgrad_solve(upglist[is],x,tempx,NULL);

	    for (k=0; k < n[1]*n[2]; k++) {
		tempr[k*n[0]] = (m[is][k]==1)? tempx[k*n[0]]: 0.;

		for (i=1; i < n[0]; i++) {
		    tempr[k*n[0]+i] = 0.;
		}
	    }

	    for (i=0; i < n[0]*n[1]*n[2]; i++) {
		r[is*n[0]*n[1]*n[2]+i] = tempr[i];
	    }
	}
    }
}
