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

#include "upgraddsr.h"
#include "dsrtomo.h"

static int nt, *m;
static upgrad upg;
static float *temp;

void dsrtomo_init(int dim  /* model dimension */,
		  int *n   /* model size */,
		  float *d /* model sampling */)
/*< initialize >*/
{
    int i;

    nt = 1;
    for (i=0; i < dim; i++) {
	nt = nt*n[i];
    }

    upg = upgrad_init(dim,n,d);

    temp = sf_floatalloc(nt);
}

void dsrtomo_set(float *t /* stencil time */,
		 float *w /* stencil slowness-squared */,
		 int *f   /* stencil flag */,
		 int *m0  /* receiver mask */)
/*< set operator >*/
{
    upgrad_set(upg,t,w,f);

    m = m0;
}

void dsrtomo_close(void)
/*< free allocated space >*/
{
    upgrad_close(upg);
}

void dsrtomo_oper(bool adj, bool add, int nx, int nr, float *x, float *r)
/*< linear operator >*/
{
    int i;

    sf_adjnull(adj,add,nx,nr,x,r);

    if (adj) {
	/* given dt solve dw */
	
	if (m != NULL) {
	    for (i=0; i < nt; i++) {
		if (m[i] != 1)
		    r[i] = 0.;
	    }
	}
	
	upgrad_paste(r);
	upgrad_inverse(upg,temp,r,NULL);
	upgrad_spread(upg,x,temp);
    } else {
	/* given dw solve dt */
	
	upgrad_collect(upg,x,temp);
	upgrad_solve(upg,temp,r,NULL);
	upgrad_copy(r);
	
	if (m != NULL) {
	    for (i=0; i < nt; i++) {
		if (m[i] != 1)
		    r[i] = 0.;
	    }
	}
    }
}
