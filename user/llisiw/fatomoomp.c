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

static int nt, **mask, ns, **list, *upgnum;
static float **tempt, **tempx, **psum, **data, *wght;
static upgrad *upglist;

void fatomo_init(int dim      /* model dimension */,
		 int *n       /* model size */,
		 float *o     /* model origin */,
		 float *d     /* model sampling */,
		 int order    /* fast march order */,
		 int nshot    /* number of shots */,
		 int **rhslist /* rhs list */,
		 int **recv   /* receiver list */,
		 float **reco /* record list */,
		 float *weight /* data weighting */)
/*< initialize >*/
{
    int i, is, mts;
    
    /* read in parameters */
    nt = 1;
    for (i=0; i < dim; i++) {
	nt = nt*n[i];
    }

    ns   = nshot;
    list = rhslist;
    mask = recv;
    data = reco;
    wght = weight;

    /* initialize upwind stencil and fast marching */
    upgrad_init(dim,n,d);
    fastmarch_init(n,o,d,order);

    /* allocate shared memory */
    upglist = (upgrad *)malloc(ns*sizeof(upgrad));

#ifdef _OPENMP    
#pragma omp parallel for
#endif
    for (is=0; is < ns; is++) {
	upglist[is] = upgrad_alloc();
    }
    
    upgnum = sf_intalloc(ns);

#ifdef _OPENMP
    mts = omp_get_max_threads();
#else
    mts = 1;
#endif

    tempt = sf_floatalloc2(nt,mts);
    tempx = sf_floatalloc2(nt,mts);
    psum  = sf_floatalloc2(nt,mts);
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
    int is;

    /* set background slowness squared */
    fastmarch_set(slow);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (is=0; is < ns; is++) {
	upgnum[is] = fastmarch(time[is],source[is],
			       list[is],mask[is],data[is],
			       rhs,upglist[is]);
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
		
		i = list[is][0];
		for (it=0; it < list[is][1]; it++) {
		    tempt[its][mask[is][it]] = (wght!=NULL)? wght[i]*r[i]: r[i];
		    i++;
		}
		
		upgrad_inverse(upgnum[is],upglist[is],tempx[its],tempt[its],NULL);
		
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
		upgrad_solve(upgnum[is],upglist[is],x,tempt[its],NULL);
		
		i = list[is][0];
		for (it=0; it < list[is][1]; it++) {
		    r[i] = tempt[its][mask[is][it]];
		    if (wght != NULL) r[i] *= wght[i];
		    i++;
		}
	    }
	}
    }
}

void fatomo_ray(float **ray)
/*< extract ray density >*/
{
    int it, is;

#ifdef _OPENMP
#pragma omp for private(it)
#endif
    for (is=0; is < ns; is++) {
	for (it=0; it <= nt; it++)
	    ray[is][it] = 0.;
	
	for (it=0; it < list[is][1]; it++) {
	    ray[is][mask[is][it]] = 1.;
	}
	
	upgrad_ray(upgnum[is],upglist[is],ray[is]);
    }
}
