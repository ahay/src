#include "conjprec.h"
#include "alloc.h"
#include "c99.h"
#include "error.h"

static int nx, nr;
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

void sf_conjprec_init(int nx1, int nr1, float eps1,
		      float tol1, bool verb1, bool hasp01) 
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
{
    free (r);
    free (sp);
    free (gp);
    free (sx);
    free (gx);
    free (sr);
    free (gr);
}
   
void sf_conjprec(sf_operator oper, sf_operator2 prec, 
		 float* p, float* x, const float* dat, int niter)
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
    r0 = verb? norm(nr,r): 0.;

    for (iter=0; iter < niter; iter++) {
	for (i=0; i < nr; i++) {
	    gp[i] = eps*p[i] + r[i];
	}

	oper(true,false,nx,nr,gx,gp);
	prec(nx,gx);
	oper(false,false,nx,nr,gx,gr);

	gn = norm(nr,gp);

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

	    for (i=0; i < nr; i++) {
		sp[i] = gp[i] + alpha * sp[i];
		sr[i] = gr[i] + alpha * sr[i];
	    }
	    for (i=0; i < nx; i++) {
		sx[i] = gx[i] + alpha * sx[i];
	    }
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
			     iter,norm(nr,r)/r0,dg);

	alpha = - gn / beta;

	for (i=0; i < nr; i++) {
	    p[i] += alpha * sp[i];
	    r[i] += alpha * sr[i];
	}

	for (i=0; i < nx; i++) {
	    x[i] += alpha * sx[i];
	}

	gnp = gn;
    }
}
