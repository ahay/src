/* Triangle smoothing */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

/*^*/
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "pblas.h"
/*^*/



struct sf_Triangle {
    float *tmp, wt;
    int np, nb, nx;
    bool box;
};

static void foldp (int o, int d, int nx, int nb, int np, 
		  const float *x, float* tmp);
static void fold2p (int o, int d, int nx, int nb, int np, 
		   float *x, const float* tmp);
static void doubintp (int nx, float *x, bool der);
static void triplep (int o, int d, int nx, int nb, 
		    float* x, const float* tmp, bool box, float wt);
static void triple2p (int o, int d, int nx, int nb, 
		     const float* x, float* tmp, bool box, float wt);

sf_triangle sf_trianglep_init (int  nbox /* triangle length */, 
			      int  ndat /* data length */,
                              bool box  /* if box instead of triangle */)
/*< initialize >*/
{
    sf_triangle tr;

    tr = (sf_triangle) sf_alloc(1,sizeof(*tr));

    tr->nx = ndat;
    tr->nb = nbox;
    tr->box = box;
    tr->np = ndat + 2*nbox;
    
    if (box) {
	tr->wt = 1.0/(2*nbox-1);
    } else {
	tr->wt = 1.0/(nbox*nbox);
    }
    
    tr->tmp = sf_floatalloc(tr->np);

    return tr;
}

static void foldp (int o, int d, int nx, int nb, int np, 
		  const float *x, float* tmp)
{
    int i, j;

    /* copy middle */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i=0; i < nx; i++) 
	tmp[i+nb] = x[o+i*d];
    
    /* reflections from the right side */
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for (j=nb+nx; j < np; j += nx) {
	for (i=0; i < nx && i < np-j; i++)
	    tmp[j+i] = x[o+(nx-1-i)*d];
	j += nx;
	for (i=0; i < nx && i < np-j; i++)
	    tmp[j+i] = x[o+i*d];
    }
    
    /* reflections from the left side */
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for (j=nb; j >= 0; j -= nx) {
	for (i=0; i < nx && i < j; i++)
	    tmp[j-1-i] = x[o+i*d];
	j -= nx;
	for (i=0; i < nx && i < j; i++)
	    tmp[j-1-i] = x[o+(nx-1-i)*d];
    }
}

static void fold2p (int o, int d, int nx, int nb, int np, 
		   float *x, const float* tmp)
{
    int i, j;

    /* copy middle */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i=0; i < nx; i++) 
	x[o+i*d] = tmp[i+nb];

    /* reflections from the right side */
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for (j=nb+nx; j < np; j += nx) {
	for (i=0; i < nx && i < np-j; i++)
	    x[o+(nx-1-i)*d] += tmp[j+i];
	j += nx;
	for (i=0; i < nx && i < np-j; i++)
	    x[o+i*d] += tmp[j+i];
    }
    
    /* reflections from the left side */
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for (j=nb; j >= 0; j -= nx) {
	for (i=0; i < nx && i < j; i++)
	    x[o+i*d] += tmp[j-1-i];
	j -= nx;
	for (i=0; i < nx && i < j; i++)
	    x[o+(nx-1-i)*d] += tmp[j-1-i];
    }
}
    
static void doubintp (int nx, float *xx, bool der)
{
    int i;
    float t;

    /* integrate backward */
    t = 0.;
    for (i=nx-1; i >= 0; i--) {
	t += xx[i];
	xx[i] = t;
    }

    if (der) return;

    /* integrate forward */
    t=0.;
    for (i=0; i < nx; i++) {
	t += xx[i];
	xx[i] = t;
    }
}

static void doubint2p (int nx, float *xx, bool der)
{
    int i;
    float t;
    int myid, threads, intl, e;
    int *assignments, *myrem;
    float split, rem;
    float *checkpoints;
    int tot = 0;
/*
 * #ifdef _OPENMP
#pragma omp parallel private(myid,t)
#endif
{

#ifdef _OPENMP
    myid = omp_get_thread_num();
#pragma omp single
{
    threads = omp_get_max_threads();
    assignments = sf_intalloc(threads+1);
    myrem = sf_intalloc(threads);
    checkpoints = sf_floatalloc(threads+1);
    split = (float)nx/threads;
    intl = floor(split);
    rem = split-(float)intl;

    for (i=0; i<threads; i++){
      if (rem*threads>(float)i){
         myrem[i] = 1;
      }else{
         myrem[i] = 0;
      }
    }
    for (i=0;i<threads;i++){
        tot = tot+myrem[i];
    }
    e = tot + intl*threads - nx;
    if (e == 0 ){}else{ sf_error("BAD OMP SPLIT IN SMOOTHING CAUSAL INTEGRATION, LOOK AT TRIANGLEP.C");}
    assignments[0] = 0;
    checkpoints[0] = 0.;
    for (i=1; i<threads+1; i++){
       assignments[i] = assignments[i-1] + intl + myrem[i-1];
    }
}

    for (i=assignments[myid]; i<assignments[myid+1]; i++){
       checkpoints[myid+1] += xx[i];
    }
//sf_warning("CHECK %g %i\n",checkpoints[myid+1],myid);
#pragma omp single
{
    for (i=1; i<threads+1; i++){
       checkpoints[i] += checkpoints[i-1];
    }
}
    t = checkpoints[myid];
    for (i=assignments[myid]; i < assignments[myid+1] ; i++){
       t += xx[i];
       xx[i] = t;
//sf_warning("%i\n",i);
    }
#endif
}

sf_warning("A");
*/
   t = 0.;
   for (i=0; i<nx; i++){
       t += xx[i];
       xx[i]  = t;
   }
    if (der) return;

    /* integrate backward */
    t = 0.;
    for (i=nx-1; i >= 0; i--) {
	t += xx[i];
	xx[i] = t;
    }
//sf_warning("B");
}

static void triplep (int o, int d, int nx, int nb, float* x, const float* tmp, bool box, float wt)
{
    int i;
    const float *tmp1, *tmp2;
    
    if (box) {
	tmp2 = tmp + 2*nb;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i=0; i < nx; i++) {
	    x[o+i*d] = (tmp[i+1] - tmp2[i])*wt;
	}
    } else {
	tmp1 = tmp + nb;
	tmp2 = tmp + 2*nb;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i=0; i < nx; i++) {
	    x[o+i*d] = (2.*tmp1[i] - tmp[i] - tmp2[i])*wt;
	}
    }
}

static void dtriplep (int o, int d, int nx, int nb, float* x, const float* tmp, float wt)
{
    int i;
    const float *tmp2;

    tmp2 = tmp + 2*nb;
#ifdef _OPENMP
#pragma omp parallel for
#endif    
    for (i=0; i < nx; i++) {
	x[o+i*d] = (tmp[i] - tmp2[i])*wt;
    }
}

static void triple2p (int o, int d, int nx, int nb, const float* x, float* tmp, bool box, float wt)
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i=0; i < nx + 2*nb; i++) {
	tmp[i] = 0;
    }

    if (box) {
	pblas_saxpy(nx,  +wt,x+o,d,tmp+1   ,1);
	pblas_saxpy(nx,  -wt,x+o,d,tmp+2*nb,1);
    } else {
	pblas_saxpy(nx,  -wt,x+o,d,tmp     ,1);
	pblas_saxpy(nx,2.*wt,x+o,d,tmp+nb  ,1);
	pblas_saxpy(nx,  -wt,x+o,d,tmp+2*nb,1);
    }
}

static void dtriple2p (int o, int d, int nx, int nb, const float* x, float* tmp, float wt)
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i=0; i < nx + 2*nb; i++) {
	tmp[i] = 0;
    }

    pblas_saxpy(nx,  wt,x+o,d,tmp     ,1);
    pblas_saxpy(nx, -wt,x+o,d,tmp+2*nb,1);
}

void sf_smoothp (sf_triangle tr  /* smoothing object */, 
		int o, int d    /* trace sampling */, 
		bool der        /* if derivative */, 
		float *x        /* data (smoothed in place) */)
/*< apply triangle smoothing >*/
{
    foldp (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp);
    doubintp (tr->np,tr->tmp,(bool) (tr->box || der));
    triplep (o,d,tr->nx,tr->nb,x,tr->tmp, tr->box, tr->wt);
}

void sf_dsmoothp (sf_triangle tr  /* smoothing object */, 
		int o, int d    /* trace sampling */, 
		bool der        /* if derivative */, 
		float *x        /* data (smoothed in place) */)
/*< apply triangle smoothing >*/
{
    foldp (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp);
    doubintp (tr->np,tr->tmp,(bool) (tr->box || der));
    dtriplep (o,d,tr->nx,tr->nb,x,tr->tmp, tr->wt);
}

void sf_smooth2p (sf_triangle tr  /* smoothing object */, 
		 int o, int d    /* trace sampling */, 
		 bool der        /* if derivative */,
		 float *x        /* data (smoothed in place) */)
/*< apply adjoint triangle smoothing >*/
{
    triple2p (o,d,tr->nx,tr->nb,x,tr->tmp, tr->box, tr->wt);
    doubint2p (tr->np,tr->tmp,(bool) (tr->box || der));
    fold2p (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp);
}

void sf_dsmooth2p (sf_triangle tr  /* smoothing object */, 
		 int o, int d    /* trace sampling */, 
		 bool der        /* if derivative */,
		 float *x        /* data (smoothed in place) */)
/*< apply adjoint triangle smoothing >*/
{
    dtriple2p (o,d,tr->nx,tr->nb,x,tr->tmp, tr->wt);
    doubint2p (tr->np,tr->tmp,(bool) (tr->box || der));
    fold2p (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp);
}

void  sf_trianglep_close(sf_triangle tr)
/*< free allocated storage >*/
{
    free (tr->tmp);
    free (tr);
}

