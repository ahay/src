/* Conjugate-gradient with shaping regularization for complex numbers. */
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


#include <math.h>

#include "cconjgrad.h"
#include "alloc.h"
#include "error.h"
#include "komplex.h"

#include "_bool.h"
#include "c99.h"
#include "_solver.h"
/*^*/

static int np, nx, nr, nd;
static sf_complex *r, *d, *sp, *sx, *sr, *gp, *gx, *gr;
static float eps, tol;
static bool verb, hasp0;

static double norm (int n, const sf_complex* x) 
/* double-precision L2 norm of a complex number */
{
    double prod, xi, yi;
    int i;

    prod = 0.;
    for (i = 0; i < n; i++) {
	xi = (double) crealf(x[i]);
	yi = (double) cimagf(x[i]);
	prod += xi*xi + yi*yi;
    }
    return prod;
}

void sf_cconjgrad_init(int np1     /* preconditioned size */, 
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

    r = sf_complexalloc(nr);  
    d = sf_complexalloc(nd); 
    sp = sf_complexalloc(np);
    gp = sf_complexalloc(np);
    sx = sf_complexalloc(nx);
    gx = sf_complexalloc(nx);
    sr = sf_complexalloc(nr);
    gr = sf_complexalloc(nr);
}

void sf_cconjgrad_close(void) 
/*< Free allocated space >*/
{
    free (r);
    free (d);
    free (sp);
    free (gp);
    free (sx);
    free (gx);
    free (sr);
    free (gr);
}

void sf_cconjgrad(sf_coperator prec     /* data preconditioning */, 
		  sf_coperator oper     /* linear operator */, 
		  sf_coperator shape    /* shaping operator */, 
		  sf_complex* p         /* preconditioned model */, 
		  sf_complex* x         /* estimated model */, 
		  const sf_complex* dat /* data */, 
		  int niter             /* number of iterations */)
/*< Conjugate gradient solver with shaping >*/
{
    double gn, gnp, alpha, beta, g0, dg, r0, b0;
    int i, iter;
    
    if (NULL != prec) {
	for (i=0; i < nd; i++) {
#ifdef SF_HAS_COMPLEX_H
	    d[i] = - dat[i];
#else
	    d[i] = sf_cneg(dat[i]);
#endif
	}
	prec(false,false,nd,nr,d,r);
    } else {
	for (i=0; i < nr; i++) {
#ifdef SF_HAS_COMPLEX_H
	    r[i] = - dat[i];
#else
	    r[i] = sf_cneg(dat[i]);
#endif
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
	for (i=0; i < np; i++) {
	    p[i] = sf_cmplx(0.0,0.0);
	}
	for (i=0; i < nx; i++) {
	    x[i] = sf_cmplx(0.0,0.0);
	}
    } 
    
    dg = g0 = b0 = gnp = 0.;
    r0 = verb? norm(nr,r): 0.;

    for (iter=0; iter < niter; iter++) {
	for (i=0; i < np; i++) {
#ifdef SF_HAS_COMPLEX_H
	    gp[i] = eps*p[i];
#else
	    gp[i] = sf_crmul(p[i],eps);
#endif
	}
	for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H	    
	    gx[i] = -eps*x[i];
#else
	    gx[i] = sf_crmul(x[i],-eps);
#endif
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

	gn = norm(np,gp);

	if (iter==0) {
	    g0 = gn;
	    b0 = fabs(gn + eps*(norm(nr,gr)-norm(nx,gx)));

	    for (i=0; i < np; i++) {
		sp[i] = gp[i];
	    }
	    for (i=0; i < nx; i++) {
		sx[i] = gx[i];
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

	    for (i=0; i < np; i++) {
#ifdef SF_HAS_COMPLEX_H	 
		sp[i] = gp[i] + alpha * sp[i];
#else
		sp[i] = sf_cadd(gp[i],sf_crmul(sp[i],alpha));
#endif
	    }
	    for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H	 		
		sx[i] = gx[i] + alpha * sx[i];
#else
		sx[i] = sf_cadd(gx[i],sf_crmul(sx[i],alpha));
#endif
	    }
	    for (i=0; i < nr; i++) {
#ifdef SF_HAS_COMPLEX_H	 		
		sr[i] = gr[i] + alpha * sr[i];
#else
		sr[i] = sf_cadd(gr[i],sf_crmul(sr[i],alpha));
#endif
	    }
	}

	beta = norm(nr,sr) + eps*(norm(np,sp) - norm(nx,sx));

	/*
	if (beta/b0 < tol) {
	    if (verb) 
		sf_warning("convergence in %d iterations, beta=%g",iter,beta);
	    break;
	}
	*/
	
	if (verb) sf_warning("iteration %d res: %f grad: %f",
			     iter,norm(nr,r)/r0,dg);

	alpha = - gn / beta;

	for (i=0; i < np; i++) {
#ifdef SF_HAS_COMPLEX_H	 
	    p[i] += alpha * sp[i];
#else
	    p[i] = sf_cadd(p[i],sf_crmul(sp[i],alpha));
#endif
	}

	for (i=0; i < nx; i++) {
#ifdef SF_HAS_COMPLEX_H
	    x[i] += alpha * sx[i];
#else
	    x[i] = sf_cadd(x[i],sf_crmul(sx[i],alpha));
#endif
	}

	for (i=0; i < nr; i++) {
#ifdef SF_HAS_COMPLEX_H
	    r[i] += alpha * sr[i];
#else
	    r[i] = sf_cadd(r[i],sf_crmul(sr[i],alpha));
#endif
	}

	gnp = gn;
    }
}

/* 	$Id: cconjgrad.c 7107 2011-04-10 02:04:14Z ivlad $	 */
