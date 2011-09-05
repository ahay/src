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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "upgradomp.h"
#include "fastmarchomp.h"
#include "fatomoomp.h"

static int nt, **mask, ns, *list, maxrecv;
static float **tempt, **tempx, **psum, **data;
static upgrad *upglist;

void fatomo_init(int dim      /* model dimension */,
		 int *n       /* model size */,
		 float *o     /* model origin */,
		 float *d     /* model sampling */,
		 int order    /* fast march order */,
		 int nshot    /* number of shots */,
		 int *rhslist /* rhs list */,
		 int **recv   /* receiver list */,
		 int nrecv    /* max recv count */,
		 float **reco /* record list */)
/*< initialize >*/
{
    int i, is, mts;

    nt = 1;
    for (i=0; i < dim; i++) {
	nt = nt*n[i];
    }

    ns = nshot;

    upgrad_init(dim,n,d);
    upglist = (upgrad *)malloc(ns*sizeof(upgrad));

#ifdef _OPENMP    
#pragma omp parallel for
#endif
    for (is=0; is < ns; is++) {
	upglist[is] = upgrad_alloc();
    }
    
#ifdef _OPENMP
    mts = omp_get_max_threads();
#else
    mts = 1;
#endif

    tempt = sf_floatalloc2(nt,mts);
    tempx = sf_floatalloc2(nt,mts);
    psum  = sf_floatalloc2(nt,mts);

    list    = rhslist;
    mask    = recv;
    maxrecv = nrecv;
    data    = reco;

    fastmarch_init(n,o,d,order);
}

void fatomo_set(float **t  /* stencil time */)
/*< set fatomo operator and right-hand side >*/
{
    int is;

    /* set stencil */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (is=0; is < ns; is++) {
	upgrad_set(upglist[is],t[is]);
    }
}

void fatomo_close(void)
/*< free allocated storage >*/
{
    int is;
    for (is=0; is < ns; is++) {
	upgrad_close(upglist[is]);
    }

    free(upglist);
    free(tempt);
    free(tempx);
    free(psum);

    fastmarch_close();
}

void fatomo_fastmarch(float *slow    /* slowness squared */,
		      float **time   /* time */,
		      float **source /* source */,
		      float *rhs     /* rhs */)
/*< fast marching >*/
{
    int is, i, it;

    fastmarch_set(slow);

#ifdef _OPENMP
#pragma omp parallel for private(i,it)
#endif
    for (is=0; is < ns; is++) {
	fastmarch(time[is],source[is]);
	
	i = list[is];
	for (it=maxrecv-1; it >= 0; it--) {
	    if (mask[is][it] >= 0) {
		rhs[i-1] = data[is][it]-time[is][mask[is][it]];
		i--;
	    }
	}
    }
}

void fatomo_lop(bool adj, bool add, int nx, int nr, float *x, float *r)
/*< linear operator >*/
{
    int it, is, i, its;

    sf_adjnull(adj,add,nx,nr,x,r);

    if (adj) {
	/* given dt solve ds */

#ifdef _OPENMP
#pragma omp parallel private(its,i,it)
#endif
	{
#ifdef _OPENMP
	    its = omp_get_thread_num();
#else
	    its = 0;
#endif

	    for (it=0; it < nt; it++) {
		psum[its][it] = 0.;
	    }

#ifdef _OPENMP
#pragma omp for
#endif
	    for (is=0; is < ns; is++) {
		for (it=0; it <= nt; it++)
		    tempt[its][it] = 0.;
		
		i = list[is];
		for (it=maxrecv-1; it >= 0; it--) {
		    if (mask[is][it] >= 0) {
			tempt[its][mask[is][it]] = r[i-1];
			i--;
		    }
		}
		
		upgrad_inverse(upglist[is],tempx[its],tempt[its],NULL);
		
		for (it=0; it < nt; it++) psum[its][it]+=tempx[its][it];
	    }
#ifdef _OPENMP
#pragma omp critical
#endif
	    {
		for (it=0; it < nt; it++) x[it]+=psum[its][it];
	    }
	}
	
    } else {
	/* given ds solve dt */

#ifdef _OPENMP	
#pragma omp parallel private(its,i,it)
#endif
	{
#ifdef _OPENMP
	    its = omp_get_thread_num();
#else
	    its = 0;
#endif

#ifdef _OPENMP
#pragma omp for
#endif
	    for (is=0; is < ns; is++) {
		upgrad_solve(upglist[is],x,tempt[its],NULL);
		
		i = list[is];
		for (it=maxrecv-1; it >= 0; it--) {
		    if (mask[is][it] >= 0) {
			r[i-1] = tempt[its][mask[is][it]];
			i--;
		    }
		}
	    }
	}
    }
}
