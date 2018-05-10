/* Conjugate-gradient with shaping regularization, OpenMP parallelism*/
/*
  Copyright (C) 2018 University of Texas at Austin
  
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

#include <math.h>
#include <float.h>

#include <rsf.h>
#include "pblas.h"
/*^*/

static int np, nx, nr, nd;
static float *r, *sp, *sx, *sr, *gp, *gx, *gr;
static float eps, tol;
static bool verb, hasp0;

void sf_conjgradp_init(int np1     /* preconditioned size */, 
		      int nx1     /* model size */, 
		      int nd1     /* data size */, 
		      int nr1     /* residual size */, 
		      float eps1  /* scaling */,
		      float tol1  /* tolerance */, 
		      bool verb1  /* verbosity flag */, 
		      bool hasp01 /* if has initial model */) 
/*< solver constructor >*/
{
    np = np1; 
    nx = nx1;
    nr = nr1;
    nd = nd1;
    eps = eps1*eps1;
    tol = tol1;
    verb = verb1;
    hasp0 = hasp01;

    r = sf_floatalloc(nr);  
    sp = sf_floatalloc(np);
    gp = sf_floatalloc(np);
    sx = sf_floatalloc(nx);
    gx = sf_floatalloc(nx);
    sr = sf_floatalloc(nr);
    gr = sf_floatalloc(nr);
}

void sf_conjgradp_close(void) 
/*< Free allocated space >*/
{
    free (r);
    free (sp);
    free (gp);
    free (sx);
    free (gx);
    free (sr);
    free (gr);
}

void sf_conjgradp(sf_operator prec  /* data preconditioning */, 
		 sf_operator oper  /* linear operator */, 
		 sf_operator shape /* shaping operator */, 
		 float* p          /* preconditioned model */, 
		 float* x          /* estimated model */, 
		 float* dat        /* data */, 
		 int niter         /* number of iterations */) 
/*< Conjugate gradient solver with shaping >*/
{
    double gn, gnp, alpha, beta, g0, dg, r0;
    float *d=NULL;
    int i, iter;
    if (NULL != prec) {
	d = sf_floatalloc(nd); 
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i=0; i < nd; i++) {
	    d[i] = - dat[i];
	}
	prec(false,false,nd,nr,d,r);
    } else {
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i=0; i < nr; i++) {
	    r[i] = - dat[i];
	}
    }
    
    if (hasp0) { /* initial p */
	shape(false,false,np,nx,p,x);
	if (NULL != prec) {
	    oper(false,false,nx,nd,x,d);
	    prec(false,true,nd,nr,d,r);
	} else {
	    oper(false,true,nx,nr,x,r);
	}
    } else {
#ifdef _OPENMP
#pragma omp parallel
#endif
{
#ifdef _OPENMP
#pragma omp for
#endif
	for (i=0; i < np; i++) {
	    p[i] = 0.;
	}
#ifdef _OPENMP
#pragma omp for
#endif
	for (i=0; i < nx; i++) {
	    x[i] = 0.;
	}
}// end parallel
    } 
    
    dg = g0 = gnp = 0.;
//pblas call
    r0 = pblas_dsdot(nr,r,1,r,1);

    if (r0 == 0.) {
	if (verb) sf_warning("zero residual: r0=%g",r0);
	return;
    }

// timings
double ot1, ot2, st1, st2, dto, dts;

#ifdef _OPENMP
ot1 = omp_get_wtime();
   oper(true,true,nx,nr,gx,r);
ot2 = omp_get_wtime();
dto = ot2-ot1;
st1 = omp_get_wtime();
   shape(true,true,np,nx,gp,gx);
st2 = omp_get_wtime();
dts = st2-st1;
#endif

    for (iter=0; iter < niter; iter++) {
#ifdef _OPENMP
#pragma omp parallel
#endif
{
#ifdef _OPENMP
#pragma omp for
#endif
	for (i=0; i < np; i++) {
	    gp[i] = eps*p[i];
	}
#ifdef _OPENMP
#pragma omp for
#endif
	for (i=0; i < nx; i++) {
	    gx[i] = -eps*x[i];
	}
}
	if (NULL != prec) {
	    prec(true,false,nd,nr,d,r);
	    oper(true,true,nx,nd,gx,d);
	} else {
	    oper(true,true,nx,nr,gx,r);
	}
	shape(true,true,np,nx,gp,gx);
	shape(false,false,np,nx,gp,gx);

	if (NULL != prec) {
	    oper(false,false,nx,nd,gx,d);
	    prec(false,false,nd,nr,d,gr);
	} else {
	    oper(false,false,nx,nr,gx,gr);
	}

// pblas call
	gn = pblas_dsdot(np,gp,1,gp,1);

	if (iter==0) {
	    g0 = gn;
#ifdef _OPENMP
#pragma omp parallel
#endif
{
#ifdef _OPENMP
#pragma omp for
#endif
	    for (i=0; i < np; i++) {
		sp[i] = gp[i];
	    }
#ifdef _OPENMP
#pragma omp for
#endif
	    for (i=0; i < nx; i++) {
		sx[i] = gx[i];
	    }
#ifdef _OPENMP
#pragma omp for
#endif
	    for (i=0; i < nr; i++) {
		sr[i] = gr[i];
	    }
} // end parallel
	} else {
	    alpha = gn / gnp;
	    dg = gn / g0;

	    if (alpha < tol || dg < tol) {
		if (verb) 
		    sf_warning(
			"convergence in %d iterations, alpha=%g, gd=%g",
			iter,alpha,dg);
#ifdef _OPENMP
sf_warning("One Shaping (triangle smoothing) Operation: %g \n",dts);
sf_warning("One Weight Operation : %g \n",dto);

#endif
		break;
	    }

// pblas calls
#ifdef _OPENMP
#pragma omp parallel
#endif
{
	    pblas_saxpyA(np,alpha,sp,1,gp,1);
	    pblas_sswapA(np,sp,1,gp,1);

	    pblas_saxpyA(nx,alpha,sx,1,gx,1);
	    pblas_sswapA(nx,sx,1,gx,1);

	    pblas_saxpyA(nr,alpha,sr,1,gr,1);
	    pblas_sswapA(nr,sr,1,gr,1);
}
	}

// pblas calls
	beta = pblas_dsdot(nr,sr,1,sr,1) + eps*(pblas_dsdot(np,sp,1,sp,1) - pblas_dsdot(nx,sx,1,sx,1));
	
	if (verb) sf_warning("iteration %d res: %f grad: %f",
			     iter,pblas_snrm2(nr,r,1)/r0,dg);

	alpha = - gn / beta;
// pblas calls
#ifdef _OPENMP
#pragma omp parallel
#endif
{
	pblas_saxpyA(np,alpha,sp,1,p,1);
	pblas_saxpyA(nx,alpha,sx,1,x,1);
	pblas_saxpyA(nr,alpha,sr,1,r,1);
}
	gnp = gn;

    }

    if (NULL != prec) free (d);

}

// not parallelized
void sf_conjgradp_adj(bool adj /* adjoint flag */,
		     sf_operator oper  /* linear operator */, 
		     sf_operator shape /* shaping operator */, 
		     float* p          /* preconditioned model */, 
		     float* x          /* estimated model */, 
		     float* dat        /* data */, 
		     int niter         /* number of iterations */) 
/*< Conjugate gradient solver with shaping and its adjoint. >*/
{
    double gn, gnp, alpha, beta, g0, dg, r0;
    float *q, *sq, *gq, *y;
    int i, iter;

    q = sf_floatalloc(nx);
    y = sf_floatalloc(nx);
    sq = sf_floatalloc(nx);
    gq = sf_floatalloc(nx);
    
    if (adj) { /* dat -> x */
	for (i=0; i < nr; i++) {
	    r[i] = - dat[i];
	}
	oper(true,false,nx,nr,q,r);
    } else { /* x -> dat */
	for (i=0; i < nx; i++) {
	    q[i] = -x[i];
	}
    }
    
    for (i=0; i < np; i++) {
	p[i] = 0.;
    }
    for (i=0; i < nx; i++) {
	y[i] = 0.;
    }
    
    dg = g0 = gnp = 0.;
    r0 = cblas_snrm2(nx,q,1);
    if (r0 == 0.) {
	if (verb) sf_warning("zero residual: r0=%g",r0);
	return;
    }

    for (iter=0; iter < niter; iter++) {
	for (i=0; i < np; i++) {
	    gp[i] = eps*p[i];
	}
	for (i=0; i < nx; i++) {
	    gx[i] = q[i]-eps*y[i];
	}

	shape(true,true,np,nx,gp,gx);
	shape(false,false,np,nx,gp,gx);

	oper(false,false,nx,nr,gx,gr);
	oper(true,false, nx,nr,gq,gr);

	gn = cblas_dsdot(np,gp,1,gp,1);

	if (iter==0) {
	    g0 = gn;

	    for (i=0; i < np; i++) {
		sp[i] = gp[i];
	    }
	    for (i=0; i < nx; i++) {
		sx[i] = gx[i];
		sq[i] = gq[i];
	    }
	    for (i=0; i < nr; i++) {
		sr[i] = gr[i];
	    }
	} else {
	    alpha = gn / gnp;
	    dg = gn / g0;

	    if (alpha < tol || dg < tol) {
		if (verb) 
		    sf_warning(
			"convergence in %d iterations, alpha=%g, gd=%g",
			iter,alpha,dg);

		break;
	    }

	    cblas_saxpy(np,alpha,sp,1,gp,1);
	    cblas_sswap(np,sp,1,gp,1);

	    cblas_saxpy(nx,alpha,sx,1,gx,1);
	    cblas_sswap(nx,sx,1,gx,1);

	    cblas_saxpy(nx,alpha,sq,1,gq,1);
	    cblas_sswap(nx,sq,1,gq,1);

	    cblas_saxpy(nr,alpha,sr,1,gr,1);
	    cblas_sswap(nr,sr,1,gr,1);
	}

	beta = cblas_dsdot(nr,sr,1,sr,1) + eps*(cblas_dsdot(np,sp,1,sp,1) - cblas_dsdot(nx,sx,1,sx,1));
	
	if (verb) sf_warning("iteration %d res: %f grad: %f",
			     iter,cblas_snrm2(nx,q,1)/r0,dg);

	alpha = - gn / beta;
	
	cblas_saxpy(np,alpha,sp,1,p,1);
	cblas_saxpy(nx,alpha,sx,1,y,1);
	cblas_saxpy(nx,alpha,sq,1,q,1);
	
	gnp = gn;
    }

    if (adj) { /* dat -> x */
	for (i=0; i < nx; i++) {
	    x[i] = y[i];
	}
    } else { /* x -> dat */
	oper(false,false,nx,nd,y,dat);
    }   

    free(q);
    free(y);
    free(sq);
    free(gq);
}


/* 	$Id$	 */
