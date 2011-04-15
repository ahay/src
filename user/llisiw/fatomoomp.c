/* First-arrival traveltime tomography (OMP) */
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
#include <omp.h>

#include "upgrad.h"
#include "fatomoomp.h"

static int nt, **mask, ns, *list;
static float *tempt, *tempx;
static upgrad *upglist;

void fatomo_init(int dim      /* model dimension */,
		 int *n       /* model size */,
		 float *d     /* model sampling */,
		 int nshot    /* number of shots */,
		 int *rhslist /* rhs list */)
/*< initialize >*/
{
    int i, is;

    nt = 1;
    for (i=0; i < dim; i++) {
	nt = nt*n[i];
    }

    ns = nshot;

    list = sf_intalloc(ns);
    tempx = sf_floatalloc(nt);

    upglist = (upgrad *)malloc(ns*sizeof(upgrad));
    
    for (is=0; is < ns; is++) {
	upglist[is] = upgrad_init(dim,n,d);
    }
    
    list = rhslist;
}

void fatomo_set(float **t  /* stencil time */,
		int **m     /* mask */)
/*< set fatomo operator and right-hand side >*/
{
    int is;

    /* set stencil */
    for (is=0; is < ns; is++) {
	upgrad_set(upglist[is],t[is]);
    }
    
    /* set mask */
    mask = m;    
}

void fatomo_close(void)
/*< free allocated storage >*/
{
    int is;
#pragma omp parallel for
    for (is=0; is < ns; is++) {
	upgrad_close(upglist[is]);
    }
}

void fatomo_lop(bool adj, bool add, int nx, int nr, float *x, float *r)
/*< linear operator >*/
{
    int it, is, i;

    sf_adjnull(adj,add,nx,nr,x,r);

    if (adj) {
	/* given dt solve ds */

#pragma omp parallel private(tempt,i,it)
	{
	    tempt = sf_floatalloc(nt);

#pragma omp for
	    for (is=0; is < ns; is++) {
		i = list[is];
		for (it=nt-1; it >= 0; it--) {
		    if (mask[is][it] == 1) {
			tempt[it] = r[i-1];
			i--;
		    } else {
			tempt[it] = 0.;
		    }
		}
		
		upgrad_inverse(upglist[is],tempx,tempt,NULL);

#pragma omp critical
		{
		    for (it=0; it < nt; it++) x[it]+=tempx[it];
		}
	    }

	    free(tempt);
	}
    } else {
	/* given ds solve dt */
	
#pragma omp parallel private(tempt,i,it)
	{
	    tempt = sf_floatalloc(nt);

#pragma omp for
	    for (is=0; is < ns; is++) {
		upgrad_solve(upglist[is],x,tempt,NULL);
		
		i = list[is];
		for (it=nt-1; it >= 0; it--) {
		    if (mask[is][it] == 1) {
			r[i-1] = tempt[it];
			i--;
		    }
		}
	    }
	    
	    free(tempt);
	}
    }
}
