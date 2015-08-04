/* First-arrival traveltime tomography. */
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

#include "fatomo.h"

static int nt, **mask, ns;
static float *tempt, *tempx;
static sf_upgrad *upglist;

void fatomo_init(int dim    /* model dimension */,
		 int *n     /* model size */,
		 float *d   /* model sampling */,
		 int nshot  /* number of shots */)
/*< initialize >*/
{
    int i, is;

    nt = 1;
    for (i=0; i < dim; i++) {
		nt = nt*n[i];
    }

    ns = nshot;

    tempt = sf_floatalloc(nt);
    tempx = sf_floatalloc(nt);

    upglist = (sf_upgrad *) sf_alloc(ns,sizeof(sf_upgrad));
    
    for (is=0; is < ns; is++) {
		upglist[is] = sf_upgrad_init(dim,n,d);
    }
}

void fatomo_set(float **t  /* stencil time */,
		int **m     /* mask */)
/*< set fatomo operator and right-hand side >*/
{
    int is;

    /* set stencil */
    for (is=0; is < ns; is++) {
		sf_upgrad_set(upglist[is],t[is]);
    }
    
    /* set mask */
    mask = m;    
}

void fatomo_close(void)
/*< free allocated storage >*/
{
    int is;

    for (is=0; is < ns; is++) {
		sf_upgrad_close(upglist[is]);
    }
}

void fatomo_lop(bool adj, bool add, int nx, int nr, float *x, float *r)
/*< linear operator >*/
{
    int it, is, count;

    sf_adjnull(adj,add,nx,nr,x,r);

    if (adj) {
		/* given dt solve ds */

		count = 0;	
		for (is=0; is < ns; is++) {
			for (it=0; it < nt; it++) {
				if (mask[is][it] == 1) {
					tempt[it] = r[count];
					count++;
				} else {
					tempt[it] = 0.;
				}
			}
	    
			sf_upgrad_inverse(upglist[is],tempx,tempt,NULL);
			for (it=0; it < nt; it++) x[it]+=tempx[it];
		}
    } else {
		/* given ds solve dt */

		count = 0;
		for (is=0; is < ns; is++) {
			sf_upgrad_solve(upglist[is],x,tempt,NULL);

			for (it=0; it < nt; it++) {
				if (mask[is][it] == 1) {
					r[count] = tempt[it];
					count++;
				}
			}
		}
    }
}
