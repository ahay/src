/* First-arrival tomography. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
#include "fatomo.h"

static int nt, *mask;
static float *tempt;
static upgrad upg;

void fatomo_init(int dim    /* model dimension */,
		 int *n     /* model size */,
		 float *d   /* model sampling */)
/*< initialize >*/
{
    int i;

    nt = 1;
    for (i=0; i < dim; i++) {
	nt = nt*n[i];
    }

    tempt = sf_floatalloc(nt);

    upg = upgrad_init(dim,n,d);
}

void fatomo_set(float *t   /* stencil time */,
		int *m     /* mask */)
/*< set fatomo operator and right-hand side >*/
{
    /* set stencil */
    upgrad_set(upg,t);

    /* set mask */
    mask = m;
}

void fatomo_close(void)
/*< free allocated storage >*/
{
    upgrad_close(upg);
}

void fatomo_lop(bool adj, bool add, int nx, int nr, float *x, float *r)
/*< linear operator >*/
{
    int it, count;

    sf_adjnull(adj,add,nx,nr,x,r);

    if (adj) {
	/* given dt solve ds */

	count = 0;
	for (it=0; it < nt; it++) {
	    if (mask[it] == 1) {
		tempt[it] = r[count];
		count++;
	    } else {
		tempt[it] = 0.;
	    }
	}
	
	upgrad_inverse(upg,x,tempt,NULL);
    } else {
	/* given ds solve dt */

	upgrad_solve(upg,x,tempt,NULL);

	count = 0;
	for (it=0; it < nt; it++) {
	    if (mask[it] == 1) {
		r[count] = tempt[it];
		count++;
	    }
	}
    }
}
