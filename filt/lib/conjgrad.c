#include <math.h>

#include "conjgrad.h"
#include "alloc.h"
#include "c99.h"
#include "error.h"

static int np, nx, nr, nd;
static float *r, *sp, *sx, *sr, *gp, *gx, *gr;
static float eps, tol;
static bool verb, hasp0;

static double norm (int n, const float* x) {
    double prod, xi;
    int i;

    prod = 0.;
    for (i = 0; i < n; i++) {
	xi = (double) x[i];
	prod += xi*xi;
    }
    return prod;
}

void sf_conjgrad_init(int np1, int nx1, int nd1, int nr1, float eps1,
		      float tol1, bool verb1, bool hasp01) 
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

void sf_conjgrad_close(void) 
{
    free (r);
    free (sp);
    free (gp);
    free (sx);
    free (gx);
    free (sr);
    free (gr);
}

/* if prec != NULL, destroys dat */   
void sf_conjgrad(sf_operator prec, sf_operator oper, sf_operator shape, 
		 float* p, float* x, float* dat, int niter) 
{
    double gn, gnp, alpha, beta, g0, dg, r0, b0;
    int i, iter;
    
    if (NULL != prec) {
	prec(false,false,nd,nr,dat,r);
	for (i=0; i < nr; i++) {
	    r[i] = - r[i];
	}
    } else {
	for (i=0; i < nr; i++) {
	    r[i] = - dat[i];
	}
    }
    
    if (hasp0) { /* initial p */
	shape(false,false,np,nx,p,x);
	if (NULL != prec) {
	    oper(false,false,nx,nd,x,dat);
	    prec(false,true,nd,nr,dat,r);
	} else {
	    oper(false,true,nx,nr,x,r);
	}
    } else {
	for (i=0; i < np; i++) {
	    p[i] = 0.;
	}
	for (i=0; i < nx; i++) {
	    x[i] = 0.;
	}
    } 
    
    dg = g0 = b0 = gnp = 0.;
    r0 = verb? norm(nr,r): 0.;

    for (iter=0; iter < niter; iter++) {
	for (i=0; i < np; i++) {
	    gp[i] = eps*p[i];
	}
	for (i=0; i < nx; i++) {
	    gx[i] = -eps*x[i];
	}

	if (NULL != prec) {
	    prec(true,false,nd,nr,dat,r);
	    oper(true,true,nx,nd,gx,dat);
	} else {
	    oper(true,true,nx,nr,gx,r);
	}

	shape(true,true,np,nx,gp,gx);
	shape(false,false,np,nx,gp,gx);

	if (NULL != prec) {
	    oper(false,false,nx,nd,gx,dat);
	    prec(false,false,nd,nr,dat,gr);
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
		sp[i] = gp[i] + alpha * sp[i];
	    }
	    for (i=0; i < nx; i++) {
		sx[i] = gx[i] + alpha * sx[i];
	    }
	    for (i=0; i < nr; i++) {
		sr[i] = gr[i] + alpha * sr[i];
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
	    p[i] += alpha * sp[i];
	}

	for (i=0; i < nx; i++) {
	    x[i] += alpha * sx[i];
	}

	for (i=0; i < nr; i++) {
	    r[i] += alpha * sr[i];
	}

	gnp = gn;
    }
}
