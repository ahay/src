/* DSR tomography */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

#include "upgraddsr0.h"
#include "dsrtomo0.h"

static int *nn, *dp, *mp;
static long nt;
static upgrad upg;
static float *temp, *pstk;

void dsrtomo_init(int dim  /* model dimension */,
		  int *n   /* model size */,
		  float *d /* model sampling */)
/*< initialize >*/
{
    int i;

    nn = n;

    nt = 1;
    for (i=0; i < dim; i++) {
	nt = nt*n[i];
    }

    upg = upgrad_init(dim,n,d);

    temp = sf_floatalloc(nt);
    pstk = sf_floatalloc(nt);
}

long* dsrtomo_order(void)
/*< upwind order >*/
{
    return upgrad_order(upg);
}

void dsrtomo_set(float *t /* stencil time */,
		 float *w /* stencil slowness-squared */,
		 int *f   /* stencil flag */,
		 int *dp0 /* receiver mask */,
		 int *mp0 /* model mask */)
/*< set operator >*/
{
    upgrad_set(upg,t,w,f);

    dp = dp0;
    mp = mp0;
}

void dsrtomo_close(void)
/*< free allocated space >*/
{
    upgrad_close(upg);
    free(temp);
}

void dsrtomo_oper(bool adj, bool add, int nx, int nr, float *x, float *r)
/*< linear operator >*/
{
    int i,j;

    sf_adjnull(adj,add,nx,nr,x,r);

    if (adj) {
	/* given dt solve dw */
	
	/* data precon */
	for (i=0; i < nn[1]*nn[2]; i++) {
	    if (dp != NULL && dp[i] != 1)
		pstk[(long) i*nn[0]] = 0.;
	    else
		pstk[(long) i*nn[0]] = r[i];
	    
	    for (j=1; j < nn[0]; j++)
		pstk[(long) i*nn[0]+j] = 0.;
	}
	
	/* linear operator */
	upgrad_inverse(upg,temp,pstk,NULL);
	upgrad_spread(upg,x,temp);

	/* model precon */
	if (mp != NULL) {
	    for (i=0; i < nn[0]*nn[1]; i++) {
		if (mp[i] != 1) x[i] = 0.;
	    }
	}
    } else {
	/* given dw solve dt */
	
	/* model precon */
	if (mp != NULL) {
	    for (i=0; i < nn[0]*nn[1]; i++) {
		if (mp[i] != 1) x[i] = 0.;
	    }
	}

	/* linear operator */
	upgrad_collect(upg,x,temp);
	upgrad_solve(upg,temp,pstk,NULL);

	/* data precon */
	for (i=0; i < nn[1]*nn[2]; i++) {
	    if (dp != NULL && dp[i] != 1)
		r[i] = 0.;
	    else
		r[i] = pstk[(long) i*nn[0]];
	}
    }
}
