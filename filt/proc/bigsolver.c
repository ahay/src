#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "chain.h"
#include "bigsolver.h"

typedef void (*weight)(int,const float*,float*);

static const float TOLERANCE=1.e-12;

static float norm (int n, const float* x);

/* solver_prec
   -------------
   Generic preconditioned linear solver.
   Solves
   oper{x} =~ dat
   eps p   =~ 0
   where x = prec{p}

   oper - operator function
   solv - solver step
   prec - preconditioning operator
   nprec - size of p (precondioned model)
   nx - model size
   ny - data size
   x[nx] - model
   dat[ny] - data
   niter - number of iterations
   eps - scaling factor */
void solver_prec (operator oper, solverstep solv, operator prec, 
		  int nprec, int nx, int ny, float* x, const float* dat, 
		  int niter, float eps, ...) {
    va_list args;
    char* par;
    float* wt = NULL;
    weight wght = NULL;
    float* x0 = NULL;
    operator nloper = NULL;
    float* mwt = NULL;
    bool verb = false;
    bool* known = NULL;
    int nmem = 0;
    int nfreq = 0;
    float** xmov = NULL;
    float** rmov = NULL;
    float* err = NULL;
    float* res = NULL;
    float* xp = NULL;
    float* wht = NULL;
    float *p, *g, *rr, *gg, *tp = NULL, *td = NULL;
    int i, iter, dprr, dppd, dppm, dpgm;
    bool forget = false;

    va_start (args, eps);
    for (;;) {
	par = va_arg (args, char *);
	if      (!strcmp (par,"end")) {break;}
	else if (!strcmp (par,"wt"))      
	{                    wt = va_arg (args, float*);}
	else if (!strcmp (par,"wght"))      
	{                    wght = va_arg (args, weight);}
	else if (!strcmp (par,"x0"))      
	{                    x0 = va_arg (args, float*);}
	else if (!strcmp (par,"nloper"))      
	{                    nloper = va_arg (args, operator);}
	else if (!strcmp (par,"mwt"))      
	{                    mwt = va_arg (args, float*);}
	else if (!strcmp (par,"verb"))      
	{                    verb = va_arg (args, int);}    
	else if (!strcmp (par,"known"))      
	{                    known = va_arg (args, bool*);}  
	else if (!strcmp (par,"nmem"))      
	{                    nmem = va_arg (args, int);}
	else if (!strcmp (par,"nfreq"))      
	{                    nfreq = va_arg (args, int);}
	else if (!strcmp (par,"xmov"))      
	{                    xmov = va_arg (args, float**);}
	else if (!strcmp (par,"rmov"))      
	{                    rmov = va_arg (args, float**);}
	else if (!strcmp (par,"err"))      
	{                    err = va_arg (args, float*);}
	else if (!strcmp (par,"res"))      
	{                    res = va_arg (args, float*);}
	else if (!strcmp (par,"xp"))      
	{                    xp = va_arg (args, float*);}
	else 
	{ exit (1);}
    }
    va_end (args);
  
    p = sf_floatalloc (ny+nprec);
    g = sf_floatalloc (ny+nprec);
    rr = sf_floatalloc (ny);
    gg = sf_floatalloc (ny);
    for (i=0; i < ny; i++) {
	rr[i] = -dat[i];
	p[i+nprec] = 0.0;
    }

    if (wt != NULL || wght != NULL) {
	td = sf_floatalloc (ny);
	if (wt != NULL) {
	    wht = wt;
	} else {
	    wht = sf_floatalloc (ny);
	    for (i=0; i < ny; i++) {
		wht[i] = 1.0;
	    }
	} 
    }

    if (mwt != NULL) tp = sf_floatalloc (nprec);

    if (x0 != NULL) {
	for (i=0; i < nprec; i++) {
	    p[i] = x0[i]; 
	}
	if (nloper != NULL) {
	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
		    tp[i] = p[i]*mwt[i];
		}
		chain (nloper, prec, false, true, nprec, ny, nx, tp, rr, x);
	    } else { 
		chain (nloper, prec, false, true, nprec, ny, nx,  p, rr, x);
	    }
	} else {
	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
		    tp[i] = p[i]*mwt[i];
		}
		chain (  oper, prec, false, true, nprec, ny, nx, tp, rr, x);
	    } else { 
		chain (  oper, prec, false, true, nprec, ny, nx,  p, rr, x);
	    }
	}
    } else {
	for (i=0; i < nprec; i++) {
	    p[i] = 0.0; 
	}
    }

    for (iter = 0; iter < niter; iter++) {
	if (nmem != 0) {
	    forget = (iter > nmem)? true: false;
	}
	if (wght != NULL && forget == true) {
	    wght (ny, rr, wht);
	}
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		rr[i] = eps*p[i+nprec] + wht[i]*rr[i];
		td[i] = rr[i]*wht[i];
	    } 
	    chain (oper, prec, true, false, nprec, ny, nx, g, td, x); 
	} else {
	    chain (oper, prec, true, false, nprec, ny, nx, g, rr, x);
	}
	if (mwt != NULL) {
	    for (i=0; i < nprec; i++) {
		g[i] = g[i]*mwt[i];
	    }
	}
	for (i=0; i < ny; i++) {
	    g[i+nprec] = eps*rr[i];
	}
	if (known != NULL) {
	    for (i=0; i < nprec; i++) {
		if (known[i] == true) {
		    g[i] = 0.0;
		} 
	    }
	}

	if (mwt != NULL) {
	    for (i=0; i < nprec; i++) {
		tp[i] = g[i]*mwt[i];
	    }
	    chain (oper, prec, false, false, nprec, ny, nx, tp, gg, x);
	} else {
	    chain (oper, prec, false, false, nprec, ny, nx,  g, gg, x);
	}
	for (i=0; i < ny; i++) {
	    if (wht != NULL) {
		gg[i] = eps*g[i+nprec] + wht[i]*gg[i];
	    } else {
		gg[i] = eps*g[i+nprec] + gg[i];
	    }
	}
	if (forget == true && nfreq != 0) {  /* periodic restart */
	    forget = (iter%nfreq == 0)? true: false;
	} 
	if (verb == true) {
	    dprr = (int) norm (ny, rr);
	    dppd = (int) norm (ny, p+nprec);
	    dppm = (int) norm (nprec, p);
	    dpgm = (int) norm (nprec, g);
	    fprintf(stderr, 
		    "iteration %d res %d prec dat %d prec mod %d grad %d\n", 
		    iter, dprr, dppd, dppm, dpgm);
	}
    
	solv (forget, nprec+ny, ny, p, g, rr, gg);
	forget = false;

	if (nloper != NULL) {
	    for (i=0; i < ny; i++) {
		rr[i] = eps*p[i+nprec] - dat[i];
	    }
	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
		    tp[i] = p[i]*mwt[i];
		}
		chain (nloper, prec, false, true, nprec, ny, nx, tp, rr, x);
	    } else { 
		chain (nloper, prec, false, true, nprec, ny, nx,  p, rr, x);
	    }
	} else if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		rr[i] = -dat[i];
	    }
	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
		    tp[i] = p[i]*mwt[i];
		}
		chain (  oper, prec, false, true, nprec, ny, nx, tp, rr, x);
	    } else { 
		chain (  oper, prec, false, true, nprec, ny, nx,  p, rr, x);
	    }	
	} else if (xmov != NULL || iter == niter-1) {
	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
		    tp[i] = p[i]*mwt[i];
		}
		prec (false, false, nprec, nx, tp, x);
	    } else {
		prec (false, false, nprec, nx,  p, x);
	    }
	}
	if (xmov != NULL) {
	    for (i=0; i < nx; i++) {
		xmov[iter][i] =  x[i];
	    }
	}
	if (rmov != NULL) {
	    for (i=0; i < ny; i++) {
		rmov[iter][i] =  p[i+nprec] * eps;
	    }
	}
    
	if (err != NULL) {
	    err[iter] = norm(ny, rr);
	}
    }

    if (xp != NULL) {
	for (i=0; i < nprec; i++) {
	    xp[i] = p[i];
	}
    }
    if (res != NULL) {
	for (i=0; i < ny; i++) {
	    res[i] = rr[i];
	}
    }
  
    free (p);
    free (g);
    free (rr);
    free (gg);

    if (wht != NULL) {
	free (td);
	if (wt == NULL) {
	    free (wht);
	}
    }

    if (mwt != NULL) {
	free (tp);
    }

}

/* solver_reg
   ----------
   Generic regularized linear solver.
   Solves
   oper{x}    =~ dat
   eps reg{x} =~ 0

   oper - operator function
   solv - solver step
   reg - regularization operator
   nreg - size of reg{x}
   nx - model size
   ny - data size
   x[nx] - model
   dat[ny] - data
   niter - number of iterations
   eps - scaling factor */
void solver_reg (operator oper, solverstep solv, operator reg, 
		 int nreg, int nx, int ny, float* x, const float* dat, 
		 int niter, float eps, ...) {

    va_list args;
    char* par;
    float* wt = NULL;
    weight wght = NULL;
    float* x0 = NULL;
    operator nloper = NULL;
    operator nlreg = NULL;
    bool verb = false;
    bool* known = NULL;
    int nmem = 0;
    int nfreq = 0;
    float** xmov = NULL;
    float** rmov = NULL;
    float* err = NULL;
    float* res = NULL;
    float* resm = NULL;
    float* wht = NULL;
    float *g, *rr, *gg, *tr, *td = NULL;
    int i, iter, dprd, dprm, dpx, dpg;
    bool forget = false;

    va_start (args, eps);
    for (;;) {
	par = va_arg (args, char *);
	if      (!strcmp (par,"end")) {break;}
	else if (!strcmp (par,"wt"))      
	{                    wt = va_arg (args, float*);}
	else if (!strcmp (par,"wght"))      
	{                    wght = va_arg (args, weight);}
	else if (!strcmp (par,"x0"))      
	{                    x0 = va_arg (args, float*);}
	else if (!strcmp (par,"nloper"))      
	{                    nloper = va_arg (args, operator);}
	else if (!strcmp (par,"nlreg"))      
	{                    nlreg = va_arg (args, operator);}
	else if (!strcmp (par,"verb"))      
	{                    verb = va_arg (args, int);}    
	else if (!strcmp (par,"known"))      
	{                    known = va_arg (args, bool*);}  
	else if (!strcmp (par,"nmem"))      
	{                    nmem = va_arg (args, int);}
	else if (!strcmp (par,"nfreq"))      
	{                    nfreq = va_arg (args, int);}
	else if (!strcmp (par,"xmov"))      
	{                    xmov = va_arg (args, float**);}
	else if (!strcmp (par,"rmov"))      
	{                    rmov = va_arg (args, float**);}
	else if (!strcmp (par,"err"))      
	{                    err = va_arg (args, float*);}
	else if (!strcmp (par,"res"))      
	{                    res = va_arg (args, float*);}
	else if (!strcmp (par,"resm"))      
	{                    resm = va_arg (args, float*);}
	else 
	{ exit (1);}
    }
    va_end (args);
 
    g =  sf_floatalloc (nx);
    tr = sf_floatalloc (nreg);
    rr = sf_floatalloc (ny+nreg);
    gg = sf_floatalloc (ny+nreg);

    if (wt != NULL || wght != NULL) {
	td = sf_floatalloc (ny);
	if (wt != NULL) {
	    wht = wt;
	} else {
	    wht = sf_floatalloc (ny);
	    for (i=0; i < ny; i++) {
		wht[i] = 1.0;
	    }
	} 
    }

    for (i=0; i < ny; i++) {
	rr[i] = - dat[i];
    }
    if (x0 != NULL) {
	for (i=0; i < nx; i++) {
	    x[i] = x0[i];
	} 
	if (nloper != NULL) {
	    nloper (false, true, nx, ny, x, rr); 
	} else {
	    oper (false, true, nx, ny, x, rr); 
	}
	if (nlreg != NULL) {
	    nlreg  (false, false, nx, ny, x, rr+ny);
	} else {
	    reg  (false, false, nx, ny, x, rr+ny);            
	}
	for (i=0; i < nreg; i++) {
	    rr[i+ny] *= eps;
	}
    } else {
	for (i=0; i < nx; i++) {
	    x[i] = 0.0;
	} 
	for (i=0; i < nreg; i++) {
	    rr[i+ny] = 0.0;
	}
    }

    for (iter=0; iter < niter; iter++) {
	if ( nmem != 0) {  /* restart */
	    forget = (iter > nmem)? true : false;
	}
	if (wght != NULL && forget == true) {
	    wght (ny, rr, wht);
	}
	for (i=0; i < nreg; i++) {
	    tr[i] = rr[i+ny]*eps;
	}
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		rr[i] *= wht[i];
		td[i] = rr[i]*wht[i];
	    }
      
	    array (oper, reg, true, false, nx, ny, nreg, g, td, tr);
	} else {
	    array (oper, reg, true, false, nx, ny, nreg, g, rr, tr);
	} 
	if (known != NULL) {
	    for (i=0; i < nx; i++) {
		if (known[i] == true) g[i] = 0.0;
	    }
	} 
	array (oper, reg, false, false, nx, ny, nreg, g, gg, gg+ny);
	for (i=0; i < nreg; i++) {
	    gg[i+ny] *= eps;
	}
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		gg[i] *= wht[i];
	    }
	}
 
	if (forget == true && nfreq != 0) { /* periodic restart */
	    forget = (iter%nfreq == 0)? true : false; 
	}

	if (verb == true) {
	    dprd = (int) norm (ny, rr);
	    dprm = (int) norm (nreg, rr+ny);
	    dpx  = (int) norm (nx, x);
	    dpg  = (int) norm (nx, g);
	    fprintf (stderr, 
		     "iteration %d res dat %d res mod %d mod %d grad %d\n",
		     iter, dprd, dprm, dpx, dpg);
	}

	solv (forget, nx, ny+nreg, x, g, rr, gg);
	forget = false;

	if (nloper != NULL) {
	    for (i=0; i < ny; i++) {
		rr[i] = -dat[i]; 
	    }
	    nloper (false, true, nx, ny, x, rr);
	}
	if (nlreg != NULL) {
	    nlreg  (false, false, nx, nreg, x, rr+ny); 
	    for (i=0; i < nreg; i++) {
		rr[i+ny] *= eps; 
	    }
	}
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		rr[i] = -dat[i]; 
	    }
	    oper (false, true, nx, ny, x, rr);
	}  
	if (xmov != NULL) {
	    for (i=0; i < nx; i++) {
		xmov[iter][i] =  x[i];
	    }
	}
	if (rmov != NULL) {
	    for (i=0; i < ny; i++) {
		rmov[iter][i] =  rr[i];
	    }
	}
    
	if (err != NULL) {
	    err[iter] = norm(ny, rr);
	}
    }

    if (resm != NULL) {
	for (i=0; i < nreg; i++) {
	    resm[i] = rr[i+ny];
	}
    }
    if (res != NULL) {
	for (i=0; i < ny; i++) {
	    res[i] = rr[i];
	}
    }

    free (tr);
    free (g); 
    free (rr);
    free (gg);

    if (wht != NULL) {
	free (td);
	if (wt == NULL) {
	    free (wht);
	}
    }
}

/* solver
   ------
   Generic linear solver.
   Solves
   oper{x} =~ dat

   oper - operator function
   solv - solver step
   nx - model size
   ny - data size
   x[nx] - model
   dat[ny] - data
   niter - number of iterations */
void solver (operator oper, solverstep solv, int nx, int ny, 
	     float* x, const float* dat, int niter, ...) {

    va_list args;
    char* par;
    float* wt = NULL;
    weight wght = NULL;
    float* x0 = NULL;
    operator nloper = NULL;
    bool verb = false;
    bool* known = NULL;
    int nmem = 0;
    int nfreq = 0;
    float** xmov = NULL;
    float** rmov = NULL;
    float* err = NULL;
    float* res = NULL;
    float* wht = NULL;
    float *g, *rr, *gg, *td = NULL;
    float dpr, dpg, dpr0, dpg0;
    int i, iter; 
    bool forget = false;

    va_start (args, niter);
    for (;;) {
	par = va_arg (args, char *);
	if      (!strcmp (par,"end")) {break;}
	else if (!strcmp (par,"wt"))      
	{                    wt = va_arg (args, float*);}
	else if (!strcmp (par,"wght"))      
	{                    wght = va_arg (args, weight);}
	else if (!strcmp (par,"x0"))      
	{                    x0 = va_arg (args, float*);}
	else if (!strcmp (par,"nloper"))      
	{                    nloper = va_arg (args, operator);}
	else if (!strcmp (par,"verb"))      
	{                    verb = va_arg (args, int);}    
	else if (!strcmp (par,"known"))      
	{                    known = va_arg (args, bool*);}  
	else if (!strcmp (par,"nmem"))      
	{                    nmem = va_arg (args, int);}
	else if (!strcmp (par,"nfreq"))      
	{                    nfreq = va_arg (args, int);}
	else if (!strcmp (par,"xmov"))      
	{                    xmov = va_arg (args, float**);}
	else if (!strcmp (par,"rmov"))      
	{                    rmov = va_arg (args, float**);}
	else if (!strcmp (par,"err"))      
	{                    err = va_arg (args, float*);}
	else if (!strcmp (par,"res"))      
	{                    res = va_arg (args, float*);}
	else 
	{ fprintf(stderr,"solver: unknown argument %s\n",par);
	exit (1);}
    }
    va_end (args);
 
    g =  sf_floatalloc (nx);
    rr = sf_floatalloc (ny);
    gg = sf_floatalloc (ny);

    if (wt != NULL || wght != NULL) {
	td = sf_floatalloc (ny);
	if (wt != NULL) {
	    wht = wt;
	} else {
	    wht = sf_floatalloc (ny);
	    for (i=0; i < ny; i++) {
		wht[i] = 1.0;
	    }
	} 
    }

    for (i=0; i < ny; i++) {
	rr[i] = - dat[i];
    }
    if (x0 != NULL) {
	for (i=0; i < nx; i++) {
	    x[i] = x0[i];
	} 
	if (nloper != NULL) {
	    nloper (false, true, nx, ny, x, rr); 
	} else {
	    oper (false, true, nx, ny, x, rr); 
	}
    } else {
	for (i=0; i < nx; i++) {
	    x[i] = 0.0;
	} 
    }

    dpr0 = norm(ny, rr);
    dpg0 = 1.;

    for (iter=0; iter < niter; iter++) {
	if ( nmem != 0) {  /* restart */
	    forget = (iter > nmem)? true : false;
	}
	if (wght != NULL && forget == true) {
	    wght (ny, rr, wht);
	}
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		rr[i] *= wht[i];
		td[i] = rr[i]*wht[i];
	    }
      
	    oper (true, false, nx, ny, g, td);
	} else {
	    oper (true, false, nx, ny, g, rr);
	} 
	if (known != NULL) {
	    for (i=0; i < nx; i++) {
		if (known[i] == true) g[i] = 0.0;
	    }
	} 
	oper (false, false, nx, ny, g, gg);
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		gg[i] *= wht[i];
	    }
	}
 
	if (forget == true && nfreq != 0) { /* periodic restart */
	    forget = (iter%nfreq == 0)? true : false; 
	}


	if (iter == 0) {
	    dpg0  = norm (nx, g);
	    dpr = 1.;
	    dpg = 1.;
	} else {
	    dpr = norm (ny, rr)/dpr0;
	    dpg = norm (nx, g)/dpg0;
	}    

	if (dpr < TOLERANCE || dpg < TOLERANCE) {
	    if (verb == true) 
		fprintf(stderr,"convergence in %d iterations\n",iter+1);
	    break;
	}

	if (verb == true) {
	    fprintf (stderr, 
		     "iteration %d res %f mod %f grad %f\n",
		     iter+1, dpr, norm (nx, x), dpg);
	}

	solv (forget, nx, ny, x, g, rr, gg);
	forget = false;

	if (nloper != NULL) {
	    for (i=0; i < ny; i++) {
		rr[i] = -dat[i]; 
	    }
	    nloper (false, true, nx, ny, x, rr);
	}
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		rr[i] = -dat[i]; 
	    }
	    oper (false, true, nx, ny, x, rr);
	}  
	if (xmov != NULL) {
	    for (i=0; i < nx; i++) {
		xmov[iter][i] =  x[i];
	    }
	}
	if (rmov != NULL) {
	    for (i=0; i < ny; i++) {
		rmov[iter][i] =  rr[i];
	    }
	}
    
	if (err != NULL) {
	    err[iter] = norm(ny, rr);
	}
    }

    if (res != NULL) {
	for (i=0; i < ny; i++) {
	    res[i] = rr[i];
	}
    }
  
    free (g);
    free (rr);
    free (gg);

    if (wht != NULL) {
	free (td);
	if (wt == NULL) {
	    free (wht);
	}
    }
}
  
static float norm (int n, const float* x) {
    int i;
    float xn;

    xn = 0.0;
    for (i=0; i < n; i++) {
	xn += x[i]*x[i];
    }
    return xn;
}

/* 	$Id: bigsolver.c,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
