/* Conjugate-gradient with preconditioning. */
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

#include "conjprec.h"
#include "alloc.h"
#include "error.h"
#include "blas.h"

#include "_bool.h"
#include "_solver.h"
/*^*/

#ifndef _sf_conjprec_h

typedef void (*sf_operator2)(int,float*);
/*^*/

#endif

static int nx, nr;
static float *r, *sp, *sx, *sr, *gp, *gx, *gr;
static float eps, tol;
static bool verb, hasp0;

void sf_conjprec_init(int nx1     /* preconditioned size */, 
		      int nr1     /* residual size */, 
		      float eps1  /* scaling */,
		      float tol1  /* tolerance */, 
		      bool verb1  /* verbosity flag */, 
		      bool hasp01 /* if has initial model */)
/*< solver constructor >*/ 
{
    nx = nx1;
    nr = nr1;
    eps = eps1*eps1;
    tol = tol1;
    verb = verb1;
    hasp0 = hasp01;

    r = sf_floatalloc(nr);    
    sp = sf_floatalloc(nr);
    gp = sf_floatalloc(nr);
    sx = sf_floatalloc(nx);
    gx = sf_floatalloc(nx);
    sr = sf_floatalloc(nr);
    gr = sf_floatalloc(nr);
}

void sf_conjprec_close(void)
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

void sf_conjprec(sf_operator oper  /* linear operator */, 
		 sf_operator2 prec /* preconditioning */, 
		 float* p          /* preconditioned */, 
		 float* x          /* model */, 
		 const float* dat  /* data */, 
		 int niter         /* number of iterations */)
/*< Conjugate gradient solver with preconditioning >*/
{
    double gn, gnp, alpha, beta, g0, dg, r0, b0;
    int i, iter;
    
    for (i=0; i < nr; i++) {
	r[i] = - dat[i];
    }
    
    if (hasp0) { /* initial p */
	oper(true,false,nx,nr,x,p);
	prec(nx,x);
	oper(false,true,nx,nr,x,r);
    } else {
	for (i=0; i < nr; i++) {
	    p[i] = 0.;
	}
	for (i=0; i < nx; i++) {
	    x[i] = 0.;
	}
    } 
    
    dg = b0 = g0 = gnp = 0.;
    r0 = verb? cblas_dsdot(nr,r,1,r,1): 0.;

    for (iter=0; iter < niter; iter++) {
	for (i=0; i < nr; i++) {
	    gp[i] = eps*p[i] + r[i];
	}

	oper(true,false,nx,nr,gx,gp);
	prec(nx,gx);
	oper(false,false,nx,nr,gx,gr);

	gn = cblas_dsdot(nr,gp,1,gp,1);

	if (iter==0) {
	    g0 = gn;

	    for (i=0; i < nr; i++) {
		alpha = (double) gp[i];
		b0 += alpha*(gr[i] + eps*alpha);
	    }

	    for (i=0; i < nr; i++) {
		sp[i] = gp[i];
		sr[i] = gr[i];
	    }
	    for (i=0; i < nx; i++) {
		sx[i] = gx[i];
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

	    cblas_saxpy(nr,alpha,sp,1,gp,1);
	    cblas_sswap(nr,sp,1,gp,1);

	    cblas_saxpy(nr,alpha,sr,1,gr,1);
	    cblas_sswap(nr,sr,1,gr,1);

	    cblas_saxpy(nx,alpha,sx,1,gx,1);
	    cblas_sswap(nx,sx,1,gx,1);
	}

	beta = 0.;
	for (i=0; i < nr; i++) {
	    alpha = (double) sp[i];
	    beta += alpha*(sr[i] + eps*alpha);
	}

	/*
	if (beta/b0 < tol) {
	    if (verb) 
		sf_warning("convergence in %d iterations, beta=%g",iter,beta);
	    break;
	}
	*/
	
	if (verb) sf_warning("iteration %d res: %f grad: %f",
			     iter,cblas_sdot(nr,r,1,r,1)/r0,dg);

	alpha = - gn / beta;

	cblas_saxpy(nr,alpha,sp,1,p,1);
	cblas_saxpy(nx,alpha,sx,1,x,1);
	cblas_saxpy(nr,alpha,sr,1,r,1);

	gnp = gn;
    }
}
