/* Solver functions for iterative least-squares optimization. */
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

#include <stdarg.h>
#include <string.h>
#include <stdlib.h>

#include "_bool.h"
#include "c99.h"
/*^*/

#include "_solver.h"
/*^*/

#include "bigsolver.h"
#include "chain.h"
#include "alloc.h"
#include "error.h"
#include "komplex.h"
#include "blas.h"

static const float TOLERANCE=1.e-12;

void sf_solver_prec (sf_operator oper   /* linear operator */, 
		     sf_solverstep solv /* stepping function */, 
		     sf_operator prec   /* preconditioning operator */, 
		     int nprec          /* size of p */, 
		     int nx             /* size of x */, 
		     int ny             /* size of dat */, 
		     float* x           /* estimated model */, 
		     const float* dat   /* data */, 
		     int niter          /* number of iterations */, 
		     float eps          /* regularization parameter */, 
		     ...                /* variable number of arguments */) 
/*< Generic preconditioned linear solver.
 ---
 Solves
 oper{x} =~ dat
 eps p   =~ 0
 where x = prec{p}
 ---
 The last parameter in the call to this function should be "end".
 Example: 
 ---
 sf_solver_prec (oper_lop,sf_cgstep,prec_lop,
 np,nx,ny,x,y,100,1.0,"x0",x0,"end");
 ---
 Parameters in ...:
 ... 
 "wt":     float*:         weight      
 "wght":   sf_weight wght: weighting function
 "x0":     float*:         initial model
 "nloper": sf_operator:    nonlinear operator  
 "mwt":    float*:         model weight
 "verb":   bool:           verbosity flag
 "known":  bool*:          known model mask
 "nmem":   int:            iteration memory
 "nfreq":  int:            periodic restart
 "xmov":   float**:        model iteration
 "rmov":   float**:        residual iteration
 "err":    float*:         final error
 "res":    float*:         final residual
 "xp":     float*:         preconditioned model
 >*/
{
    va_list args;
    char* par;
    float* wt = NULL;
    sf_weight wght = NULL;
    float* x0 = NULL;
    sf_operator nloper = NULL;
    float* mwt = NULL;
    bool verb = false;
    bool* known = NULL;
    int nmem = -1;
    int nfreq = 0;
    float** xmov = NULL;
    float** rmov = NULL;
    float* err = NULL;
    float* res = NULL;
    float* xp = NULL;
    float* wht = NULL;
    float *p, *g, *rr, *gg, *tp = NULL, *td = NULL;
    int i, iter;
    double dprr, dppd, dppm, dpgm, dprr0=1., dpgm0=1.;
    bool forget = false;

    va_start (args, eps);
    for (;;) {
	par = va_arg (args, char *);
	if      (0 == strcmp (par,"end")) {break;}
	else if (0 == strcmp (par,"wt"))      
	{                    wt = va_arg (args, float*);}
	else if (0 == strcmp (par,"wght"))      
	{                    wght = va_arg (args, sf_weight);}
	else if (0 == strcmp (par,"x0"))      
	{                    x0 = va_arg (args, float*);}
	else if (0 == strcmp (par,"nloper"))      
	{                    nloper = va_arg (args, sf_operator);}
	else if (0 == strcmp (par,"mwt"))      
	{                    mwt = va_arg (args, float*);}
	else if (0 == strcmp (par,"verb"))      
	{                    verb = (bool) va_arg (args, int);}    
	else if (0 == strcmp (par,"known"))      
	{                    known = va_arg (args, bool*);}  
	else if (0 == strcmp (par,"nmem"))      
	{                    nmem = va_arg (args, int);}
	else if (0 == strcmp (par,"nfreq"))      
	{                    nfreq = va_arg (args, int);}
	else if (0 == strcmp (par,"xmov"))      
	{                    xmov = va_arg (args, float**);}
	else if (0 == strcmp (par,"rmov"))      
	{                    rmov = va_arg (args, float**);}
	else if (0 == strcmp (par,"err"))      
	{                    err = va_arg (args, float*);}
	else if (0 == strcmp (par,"res"))      
	{                    res = va_arg (args, float*);}
	else if (0 == strcmp (par,"xp"))      
	{                    xp = va_arg (args, float*);}
	else 
	{ sf_error("%s: unknown parameter %s",__FILE__,par);}
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
		sf_chain (nloper, prec, false, true, nprec, ny, nx, tp, rr, x);
	    } else { 
		sf_chain (nloper, prec, false, true, nprec, ny, nx,  p, rr, x);
	    }
	} else {
	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
		    tp[i] = p[i]*mwt[i];
		}
		sf_chain (  oper, prec, false, true, nprec, ny, nx, tp, rr, x);
	    } else { 
		sf_chain (  oper, prec, false, true, nprec, ny, nx,  p, rr, x);
	    }
	}
    } else {
	for (i=0; i < nprec; i++) {
	    p[i] = 0.0; 
	}
    }

    for (iter = 0; iter < niter; iter++) {
	if (nmem >= 0) {
	    forget = (bool) (iter >= nmem);
	}
	if (wght != NULL && forget) {
	    wght (ny, rr, wht);
	}
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		rr[i] = eps*p[i+nprec] + wht[i]*rr[i];
		td[i] = rr[i]*wht[i];
	    } 
	    sf_chain (oper, prec, true, false, nprec, ny, nx, g, td, x); 
	} else {
	    sf_chain (oper, prec, true, false, nprec, ny, nx, g, rr, x);
	}
	if (mwt != NULL) {
	    for (i=0; i < nprec; i++) {
		g[i] *= mwt[i];
	    }
	}
	for (i=0; i < ny; i++) {
	    g[i+nprec] = eps*rr[i];
	}
	if (known != NULL) {
	    for (i=0; i < nprec; i++) {
		if (known[i]) {
		    g[i] = 0.0;
		} 
	    }
	}

	if (mwt != NULL) {
	    for (i=0; i < nprec; i++) {
		tp[i] = g[i]*mwt[i];
	    }
	    sf_chain (oper, prec, false, false, nprec, ny, nx, tp, gg, x);
	} else {
	    sf_chain (oper, prec, false, false, nprec, ny, nx,  g, gg, x);
	}
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		gg[i] *= wht[i];
	    }
	}
	cblas_saxpy(ny,eps,g+nprec,1,gg,1);

	if (forget && nfreq != 0) {  /* periodic restart */
	    forget = (bool) (0 == (iter+1)%nfreq);
	} 
	
	if (iter == 0) {
	    dprr0 = cblas_snrm2 (ny, rr, 1);
	    dpgm0 = cblas_snrm2 (nprec, g, 1);
	    dprr = 1.;
	    dpgm = 1.;
	} else {
	    dprr = cblas_snrm2 (ny, rr, 1)/dprr0;
	    dpgm = cblas_snrm2 (nprec, g, 1)/dpgm0;
	}
	dppd = cblas_snrm2 (ny, p+nprec, 1);
	dppm = cblas_snrm2 (nprec, p, 1);

	if (verb) 
	    sf_warning("iteration %d res %g prec dat %g prec mod %g grad %g", 
		       iter, dprr, dppd, dppm, dpgm);
	
	if (dprr < TOLERANCE || dpgm < TOLERANCE) {
	    if (verb) 
		sf_warning("convergence in %d iterations",iter+1);

	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
		    tp[i] = p[i]*mwt[i];
		}
		prec (false, false, nprec, nx, tp, x);
	    } else {
		prec (false, false, nprec, nx,  p, x);
	    }

	    break;
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
		sf_chain (nloper, prec, false, true, nprec, ny, nx, tp, rr, x);
	    } else { 
		sf_chain (nloper, prec, false, true, nprec, ny, nx,  p, rr, x);
	    }
	} else if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		rr[i] = -dat[i];
	    }
	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
		    tp[i] = p[i]*mwt[i];
		}
		sf_chain (  oper, prec, false, true, nprec, ny, nx, tp, rr, x);
	    } else { 
		sf_chain (  oper, prec, false, true, nprec, ny, nx,  p, rr, x);
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
	if (err != NULL) err[iter] = cblas_snrm2(ny, rr, 1);
    } /* iter */

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

    for (; iter < niter; iter++) {
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
	if (err != NULL) err[iter] = cblas_snrm2(ny, rr, 1);
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

void sf_csolver_prec (sf_coperator oper   /* linear operator */, 
		      sf_csolverstep solv /* stepping function */, 
		      sf_coperator prec   /* preconditioning operator */, 
		      int nprec          /* size of p */, 
		      int nx             /* size of x */, 
		      int ny             /* size of dat */, 
		      sf_complex* x           /* estimated model */, 
		      const sf_complex* dat   /* data */, 
		      int niter          /* number of iterations */, 
		      float eps          /* regularization parameter */, 
		     ...                /* variable number of arguments */) 
/*< Generic preconditioned linear solver.
 ---
 Solves
 oper{x} =~ dat
 eps p   =~ 0
 where x = prec{p}
 ---
 The last parameter in the call to this function should be "end".
 Example: 
 ---
 sf_csolver_prec (oper_lop,sf_cgstep,prec_lop,
 np,nx,ny,x,y,100,1.0,"x0",x0,"end");
 ---
 Parameters in ...:
 ... 
 "wt":     float*:         weight      
 "wght":   sf_cweight wght: weighting function
 "x0":     sf_complex*:         initial model
 "nloper": sf_coperator:   nonlinear operator  
 "mwt":    float*:         model weight
 "verb":   bool:           verbosity flag
 "known":  bool*:          known model mask
 "nmem":   int:            iteration memory
 "nfreq":  int:            periodic restart
 "xmov":   sf_complex**:        model iteration
 "rmov":   sf_complex**:        residual iteration
 "err":    sf_float*:         final error
 "res":    sf_complex*:         final residual
 "xp":     sf_complex*:         preconditioned model
 >*/
{
    va_list args;
    char* par;
    float* wt = NULL;
    sf_cweight wght = NULL;
    sf_complex* x0 = NULL;
    sf_coperator nloper = NULL;
    float* mwt = NULL;
    bool verb = false;
    bool* known = NULL;
    int nmem = -1;
    int nfreq = 0;
    sf_complex** xmov = NULL;
    sf_complex** rmov = NULL;
    float* err = NULL;
    sf_complex* res = NULL;
    sf_complex* xp = NULL;
    float* wht = NULL;
    sf_complex *p, *g, *rr, *gg, *tp = NULL, *td = NULL;
    int i, iter;
    double dprr, dppd, dppm, dpgm, dprr0=1., dpgm0=1.;
    bool forget = false;

    va_start (args, eps);
    for (;;) {
	par = va_arg (args, char *);
	if      (0 == strcmp (par,"end")) {break;}
	else if (0 == strcmp (par,"wt"))      
	{                    wt = va_arg (args, float*);}
	else if (0 == strcmp (par,"wght"))      
	{                    wght = va_arg (args, sf_cweight);}
	else if (0 == strcmp (par,"x0"))      
	{                    x0 = va_arg (args, sf_complex*);}
	else if (0 == strcmp (par,"nloper"))      
	{                    nloper = va_arg (args, sf_coperator);}
	else if (0 == strcmp (par,"mwt"))      
	{                    mwt = va_arg (args, float*);}
	else if (0 == strcmp (par,"verb"))      
	{                    verb = (bool) va_arg (args, int);}    
	else if (0 == strcmp (par,"known"))      
	{                    known = va_arg (args, bool*);}  
	else if (0 == strcmp (par,"nmem"))      
	{                    nmem = va_arg (args, int);}
	else if (0 == strcmp (par,"nfreq"))      
	{                    nfreq = va_arg (args, int);}
	else if (0 == strcmp (par,"xmov"))      
	{                    xmov = va_arg (args, sf_complex**);}
	else if (0 == strcmp (par,"rmov"))      
	{                    rmov = va_arg (args, sf_complex**);}
	else if (0 == strcmp (par,"err"))      
	{                    err = va_arg (args, float*);}
	else if (0 == strcmp (par,"res"))      
	{                    res = va_arg (args, sf_complex*);}
	else if (0 == strcmp (par,"xp"))      
	{                    xp = va_arg (args, sf_complex*);}
	else 
	{ sf_error("%s: unknown parameter %s",__FILE__,par);}
    }
    va_end (args);
  
    p = sf_complexalloc (ny+nprec);
    g = sf_complexalloc (ny+nprec);
    rr = sf_complexalloc (ny);
    gg = sf_complexalloc (ny);
    for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
	rr[i] = -dat[i];
#else
	rr[i] = sf_cneg(dat[i]);
#endif	
	p[i+nprec] = sf_cmplx(0.0,0.0);
    }

    if (wt != NULL || wght != NULL) {
	td = sf_complexalloc (ny);
	if (wt != NULL) {
	    wht = wt;
	} else {
	    wht = sf_floatalloc (ny);
	    for (i=0; i < ny; i++) {
		wht[i] = 1.0;
	    }
	} 
    }

    if (mwt != NULL) tp = sf_complexalloc (nprec);

    if (x0 != NULL) {
	for (i=0; i < nprec; i++) {
	    p[i] = x0[i]; 
	}
	if (nloper != NULL) {
	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
#ifdef SF_HAS_COMPLEX_H
		    tp[i] = p[i]*mwt[i];
#else
		    tp[i] = sf_crmul(p[i],mwt[i]);
#endif
		}
		sf_cchain (nloper, prec, false, true, nprec, ny, nx, tp, rr, x);
	    } else { 
		sf_cchain (nloper, prec, false, true, nprec, ny, nx,  p, rr, x);
	    }
	} else {
	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
#ifdef SF_HAS_COMPLEX_H
		    tp[i] = p[i]*mwt[i];
#else
		    tp[i] = sf_crmul(p[i],mwt[i]);
#endif
		}
		sf_cchain (  oper, prec, false, true, nprec, ny, nx, tp, rr, x);
	    } else { 
		sf_cchain (  oper, prec, false, true, nprec, ny, nx,  p, rr, x);
	    }
	}
    } else {
	for (i=0; i < nprec; i++) {
	    p[i] = sf_cmplx(0.0,0.0); 
	}
    }

    for (iter = 0; iter < niter; iter++) {
	if (nmem >= 0) {
	    forget = (bool) (iter >= nmem);
	}
	if (wght != NULL && forget) {
	    wght (ny, rr, wht);
	}
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
		rr[i] = eps*p[i+nprec] + wht[i]*rr[i];
		td[i] = rr[i]*wht[i];
#else
		rr[i] = sf_cadd(sf_crmul(p[i+nprec],eps),
				sf_crmul(rr[i],wht[i]));
		td[i] = sf_crmul(rr[i],wht[i]);
#endif
	    } 
	    sf_cchain (oper, prec, true, false, nprec, ny, nx, g, td, x); 
	} else {
	    sf_cchain (oper, prec, true, false, nprec, ny, nx, g, rr, x);
	}
	if (mwt != NULL) {
	    for (i=0; i < nprec; i++) {
#ifdef SF_HAS_COMPLEX_H
		g[i] *= mwt[i];
#else
		g[i] = sf_crmul(g[i],mwt[i]);
#endif
	    }
	}
	for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
	    g[i+nprec] = eps*rr[i];
#else
	    g[i+nprec] = sf_crmul(rr[i],eps);
#endif
	}
	if (known != NULL) {
	    for (i=0; i < nprec; i++) {
		if (known[i]) {
		    g[i] = sf_cmplx(0.0,0.0);
		} 
	    }
	}

	if (mwt != NULL) {
	    for (i=0; i < nprec; i++) {
#ifdef SF_HAS_COMPLEX_H
		tp[i] = g[i]*mwt[i];
#else
		tp[i] = sf_crmul(g[i],mwt[i]);
#endif
	    }
	    sf_cchain (oper, prec, false, false, nprec, ny, nx, tp, gg, x);
	} else {
	    sf_cchain (oper, prec, false, false, nprec, ny, nx,  g, gg, x);
	}
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
		gg[i] *= wht[i];
#else
		gg[i] = sf_crmul(gg[i],wht[i]);
#endif
	    }
	}
        /*	cblas_saxpy(ny,eps,g+nprec,1,gg,1); */

	for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
	    gg[i] += eps * g[i+nprec];
#else
	    gg[i] = sf_cadd(gg[i],sf_crmul(g[i+nprec],eps));
#endif
	}

	if (forget && nfreq != 0) {  /* periodic restart */
	    forget = (bool) (0 == (iter+1)%nfreq);
	} 
	
	if (iter == 0) {
	    dprr0 = cblas_scnrm2 (ny, rr, 1);
	    dpgm0 = cblas_scnrm2 (nprec, g, 1);
	    dprr = 1.;
	    dpgm = 1.;
	} else {
	    dprr = cblas_scnrm2 (ny, rr, 1)/dprr0;
	    dpgm = cblas_scnrm2 (nprec, g, 1)/dpgm0;
	}
	dppd = cblas_scnrm2 (ny, p+nprec, 1);
	dppm = cblas_scnrm2 (nprec, p, 1);

	if (verb) 
	    sf_warning("iteration %d res %g prec dat %g prec mod %g grad %g", 
		       iter, dprr, dppd, dppm, dpgm);
	
	if (dprr < TOLERANCE || dpgm < TOLERANCE) {
	    if (verb) 
		sf_warning("convergence in %d iterations",iter+1);

	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
#ifdef SF_HAS_COMPLEX_H
		    tp[i] = p[i]*mwt[i];
#else
		    tp[i] = sf_crmul(p[i],mwt[i]);
#endif
		}
		prec (false, false, nprec, nx, tp, x);
	    } else {
		prec (false, false, nprec, nx,  p, x);
	    }

	    break;
	}

	solv (forget, nprec+ny, ny, p, g, rr, gg);
	forget = false;

	if (nloper != NULL) {
	    for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
		rr[i] = eps*p[i+nprec] - dat[i];
#else
		rr[i] = sf_cadd(sf_crmul(p[i+nprec],eps),sf_cneg(dat[i]));
#endif
	    }
	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
#ifdef SF_HAS_COMPLEX_H
		    tp[i] = p[i]*mwt[i];
#else
		    tp[i] = sf_crmul(p[i],mwt[i]);
#endif
		}
		sf_cchain (nloper, prec, false, true, nprec, ny, nx, tp, rr, x);
	    } else { 
		sf_cchain (nloper, prec, false, true, nprec, ny, nx,  p, rr, x);
	    }
	} else if (wht != NULL) {
	    for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
		rr[i] = -dat[i];
#else
		rr[i] = sf_cneg(dat[i]);
#endif
	    }
	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
#ifdef SF_HAS_COMPLEX_H
		    tp[i] = p[i]*mwt[i];
#else
		    tp[i] = sf_crmul(p[i],mwt[i]);
#endif
		}
		sf_cchain (  oper, prec, false, true, nprec, ny, nx, tp, rr, x);
	    } else { 
		sf_cchain (  oper, prec, false, true, nprec, ny, nx,  p, rr, x);
	    }	
	} else if (xmov != NULL || iter == niter-1) {
	    if (mwt != NULL) {
		for (i=0; i < nprec; i++) {
#ifdef SF_HAS_COMPLEX_H
		    tp[i] = p[i]*mwt[i];
#else
		    tp[i] = sf_crmul(p[i],mwt[i]);
#endif
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
#ifdef SF_HAS_COMPLEX_H
		rmov[iter][i] =  p[i+nprec] * eps;
#else
		rmov[iter][i] =  sf_crmul(p[i+nprec],eps);
#endif
	    }
	}
	if (err != NULL) err[iter] = cblas_scnrm2(ny, rr, 1);
    } /* iter */

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

    for (; iter < niter; iter++) {
	if (xmov != NULL) {
	    for (i=0; i < nx; i++) {
		xmov[iter][i] =  x[i];
	    }
	}
	if (rmov != NULL) {
	    for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
		rmov[iter][i] =  p[i+nprec] * eps;
#else
		rmov[iter][i] =  sf_crmul(p[i+nprec],eps);
#endif
	    }
	}    
	if (err != NULL) err[iter] = cblas_scnrm2(ny, rr, 1);
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

void sf_solver_reg (sf_operator oper   /* linear operator */, 
		    sf_solverstep solv /* stepping function */,
		    sf_operator reg    /* regularization operator */, 
		    int nreg           /* size of reg{x} */, 
		    int nx             /* size of x */, 
		    int ny             /* size of dat */, 
		    float* x           /* estimated model */, 
		    const float* dat   /* data */, 
		    int niter          /* number of iterations */, 
		    float eps          /* regularization parameter */, 
		    ...                /* variable number of arguments */) 
/*< Generic regularized linear solver.
  ---
  Solves
  oper{x}    =~ dat
  eps reg{x} =~ 0
  ---
  The last parameter in the call to this function should be "end".
  Example: 
  ---
  sf_solver_reg (oper_lop,sf_cgstep,reg_lop,
  np,nx,ny,x,y,100,1.0,"x0",x0,"end");
  ---
  Parameters in ...:
  
  "wt":     float*:         weight      
  "wght":   sf_weight wght: weighting function
  "x0":     float*:         initial model
  "nloper": sf_operator:    nonlinear operator  
  "nlreg":  sf_operator:    nonlinear regularization operator
  "verb":   bool:           verbosity flag
  "known":  bool*:          known model mask
  "nmem":   int:            iteration memory
  "nfreq":  int:            periodic restart
  "xmov":   float**:        model iteration
  "rmov":   float**:        residual iteration
  "err":    float*:         final error
  "res":    float*:         final residual
  "resm":   float*:         final model residual
  >*/
{

    va_list args;
    char* par;
    float* wt = NULL;
    sf_weight wght = NULL;
    float* x0 = NULL;
    sf_operator nloper = NULL;
    sf_operator nlreg = NULL;
    bool verb = false;
    bool* known = NULL;
    int nmem = -1;
    int nfreq = 0;
    float** xmov = NULL;
    float** rmov = NULL;
    float* err = NULL;
    float* res = NULL;
    float* resm = NULL;
    float* wht = NULL;
    float *g, *rr, *gg, *tr, *td = NULL;
    float dpr, dpg, dpr0, dpg0;
    int i, iter;
    bool forget = false;

    va_start (args, eps);
    for (;;) {
	par = va_arg (args, char *);
	if      (0 == strcmp (par,"end")) {break;}
	else if (0 == strcmp (par,"wt"))      
	{                    wt = va_arg (args, float*);}
	else if (0 == strcmp (par,"wght"))      
	{                    wght = va_arg (args, sf_weight);}
	else if (0 == strcmp (par,"x0"))      
	{                    x0 = va_arg (args, float*);}
	else if (0 == strcmp (par,"nloper"))      
	{                    nloper = va_arg (args, sf_operator);}
	else if (0 == strcmp (par,"nlreg"))      
	{                    nlreg = va_arg (args, sf_operator);}
	else if (0 == strcmp (par,"verb"))      
	{                    verb = (bool) va_arg (args, int);}    
	else if (0 == strcmp (par,"known"))      
	{                    known = va_arg (args, bool*);}  
	else if (0 == strcmp (par,"nmem"))      
	{                    nmem = va_arg (args, int);}
	else if (0 == strcmp (par,"nfreq"))      
	{                    nfreq = va_arg (args, int);}
	else if (0 == strcmp (par,"xmov"))      
	{                    xmov = va_arg (args, float**);}
	else if (0 == strcmp (par,"rmov"))      
	{                    rmov = va_arg (args, float**);}
	else if (0 == strcmp (par,"err"))      
	{                    err = va_arg (args, float*);}
	else if (0 == strcmp (par,"res"))      
	{                    res = va_arg (args, float*);}
	else if (0 == strcmp (par,"resm"))      
	{                    resm = va_arg (args, float*);}
	else 
	{ sf_error ("%s: unknown parameter %s",__FILE__,par);}
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
	    nlreg  (false, false, nx, nreg, x, rr+ny);
	} else {
	    reg  (false, false, nx, nreg, x, rr+ny);            
	}
	cblas_sscal(nreg,eps,rr+ny,1);
    } else {
	for (i=0; i < nx; i++) {
	    x[i] = 0.0;
	} 
	for (i=0; i < nreg; i++) {
	    rr[i+ny] = 0.0;
	}
    }

    dpr0 = cblas_snrm2(ny, rr, 1);
    dpg0 = 1.;

    for (iter=0; iter < niter; iter++) {
	if ( nmem >= 0) {  /* restart */
	    forget = (bool) (iter >= nmem);
	}
	if (wght != NULL && forget) {
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
      
	    sf_array (oper, reg, true, false, nx, ny, nreg, g, td, tr);
	} else {
	    sf_array (oper, reg, true, false, nx, ny, nreg, g, rr, tr);
	} 
	if (known != NULL) {
	    for (i=0; i < nx; i++) {
		if (known[i]) g[i] = 0.0;
	    }
	} 
	sf_array (oper, reg, false, false, nx, ny, nreg, g, gg, gg+ny);

	cblas_sscal(nreg,eps,gg+ny,1);

	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		gg[i] *= wht[i];
	    }
	}
 
	if (forget && nfreq != 0) { /* periodic restart */
	    forget = (bool) (0 == (iter+1)%nfreq);
	}

	if (iter == 0) {
	    dpg0  = cblas_snrm2 (nx, g, 1);
	    dpr = 1.;
	    dpg = 1.;
	} else {
	    dpr = cblas_snrm2 (ny, rr, 1)/dpr0;	    
	    dpg = cblas_snrm2 (nx, g , 1)/dpg0;
	} 

	if (verb) 
	    sf_warning ("iteration %d res dat %f res mod %f mod %f grad %f",
			iter, dpr, cblas_snrm2 (nreg, rr+ny, 1), cblas_snrm2 (nx, x, 1), dpg);

	if (dpr < TOLERANCE || dpg < TOLERANCE) {
	    if (verb) 
		sf_warning("convergence in %d iterations",iter+1);
	    break;
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

	    cblas_sscal(nreg,eps,rr+ny,1);
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
	    err[iter] = cblas_snrm2(ny, rr, 1);
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

void sf_csolver_reg (sf_coperator oper   /* linear operator */, 
		     sf_csolverstep solv /* stepping function */,
		     sf_coperator reg    /* regularization operator */, 
		     int nreg            /* size of reg{x} */, 
		     int nx              /* size of x */, 
		     int ny              /* size of dat */, 
		     sf_complex* x           /* estimated model */, 
		     const sf_complex* dat   /* data */, 
		     int niter          /* number of iterations */, 
		     float eps          /* regularization parameter */, 
		    ...                /* variable number of arguments */) 
/*< Generic regularized linear solver.
  ---
  Solves
  oper{x}    =~ dat
  eps reg{x} =~ 0
  ---
  The last parameter in the call to this function should be "end".
  Example: 
  ---
  sf_solver_reg (oper_lop,sf_cgstep,reg_lop,
  np,nx,ny,x,y,100,1.0,"x0",x0,"end");
  ---
  Parameters in ...:
  
  "x0":     sf_complex*:         initial model
  "nloper": sf_operator:    nonlinear operator  
  "nlreg":  sf_operator:    nonlinear regularization operator
  "verb":   bool:           verbosity flag
  "known":  bool*:          known model mask
  "nmem":   int:            iteration memory
  "nfreq":  int:            periodic restart
  "xmov":   sf_complex**:        model iteration
  "rmov":   sf_complex**:        residual iteration
  "err":    float*:         final error
  "res":    sf_complex*:         final residual
  "resm":   sf_complex*:         final model residual
  >*/
{

    va_list args;
    char* par;
    sf_complex* x0 = NULL;
    sf_coperator nloper = NULL;
    sf_coperator nlreg = NULL;
    bool verb = false;
    bool* known = NULL;
    int nmem = -1;
    int nfreq = 0;
    sf_complex** xmov = NULL;
    sf_complex** rmov = NULL;
    float* err = NULL;
    sf_complex* res = NULL;
    sf_complex* resm = NULL;
    sf_complex *g, *rr, *gg, *tr;
    float dpr, dpg, dpr0, dpg0;
    int i, iter;
    bool forget = false;

    va_start (args, eps);
    for (;;) {
	par = va_arg (args, char *);
	if      (0 == strcmp (par,"end")) {break;}
	else if (0 == strcmp (par,"x0"))      
	{                    x0 = va_arg (args, sf_complex*);}
	else if (0 == strcmp (par,"nloper"))      
	{                    nloper = va_arg (args, sf_coperator);}
	else if (0 == strcmp (par,"nlreg"))      
	{                    nlreg = va_arg (args, sf_coperator);}
	else if (0 == strcmp (par,"verb"))      
	{                    verb = (bool) va_arg (args, int);}    
	else if (0 == strcmp (par,"known"))      
	{                    known = va_arg (args, bool*);}  
	else if (0 == strcmp (par,"nmem"))      
	{                    nmem = va_arg (args, int);}
	else if (0 == strcmp (par,"nfreq"))      
	{                    nfreq = va_arg (args, int);}
	else if (0 == strcmp (par,"xmov"))      
	{                    xmov = va_arg (args, sf_complex**);}
	else if (0 == strcmp (par,"rmov"))      
	{                    rmov = va_arg (args, sf_complex**);}
	else if (0 == strcmp (par,"err"))      
	{                    err = va_arg (args, float*);}
	else if (0 == strcmp (par,"res"))      
	{                    res = va_arg (args, sf_complex*);}
	else if (0 == strcmp (par,"resm"))      
	{                    resm = va_arg (args, sf_complex*);}
	else 
	{ sf_error ("%s: unknown parameter %s",__FILE__,par);}
    }
    va_end (args);
 
    g =  sf_complexalloc (nx);
    tr = sf_complexalloc (nreg);
    rr = sf_complexalloc (ny+nreg);
    gg = sf_complexalloc (ny+nreg);

    for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
	rr[i] = - dat[i];
#else
	rr[i] = sf_cneg(dat[i]);
#endif
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
	    nlreg  (false, false, nx, nreg, x, rr+ny);
	} else {
	    reg  (false, false, nx, nreg, x, rr+ny);            
	}
	cblas_csscal(nreg,eps,rr+ny,1);
    } else {
	for (i=0; i < nx; i++) {
	    x[i] = sf_cmplx(0.0,0.0);
	} 
	for (i=0; i < nreg; i++) {
	    rr[i+ny] = sf_cmplx(0.0,0.0);
	}
    }

    dpr0 = cblas_scnrm2(ny, rr, 1);
    dpg0 = 1.;

    for (iter=0; iter < niter; iter++) {
	if ( nmem >= 0) {  /* restart */
	    forget = (bool) (iter >= nmem);
	}
	for (i=0; i < nreg; i++) {
#ifdef SF_HAS_COMPLEX_H
	    tr[i] = rr[i+ny]*eps;
#else
	    tr[i] = sf_crmul(rr[i+ny],eps);
#endif
	}
	sf_carray (oper, reg, true, false, nx, ny, nreg, g, rr, tr);
	if (known != NULL) {
	    for (i=0; i < nx; i++) {
		if (known[i]) g[i] = sf_cmplx(0.0,0.0);
	    }
	} 
	sf_carray (oper, reg, false, false, nx, ny, nreg, g, gg, gg+ny);

	cblas_csscal(nreg,eps,gg+ny,1);

	if (forget && nfreq != 0) { /* periodic restart */
	    forget = (bool) (0 == (iter+1)%nfreq);
	}

	if (iter == 0) {
	    dpg0  = cblas_scnrm2 (nx, g, 1);
	    dpr = 1.;
	    dpg = 1.;
	} else {
	    dpr = cblas_scnrm2 (ny, rr, 1)/dpr0;	    
	    dpg = cblas_scnrm2 (nx, g , 1)/dpg0;
	} 

	if (verb) 
	    sf_warning ("iteration %d res dat %f res mod %f mod %f grad %f",
			iter, dpr, cblas_scnrm2 (nreg, rr+ny, 1), cblas_scnrm2 (nx, x, 1), dpg);

	if (dpr < TOLERANCE || dpg < TOLERANCE) {
	    if (verb) 
		sf_warning("convergence in %d iterations",iter+1);
	    break;
	}

	solv (forget, nx, ny+nreg, x, g, rr, gg);
	forget = false;

	if (nloper != NULL) {
	    for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
		rr[i] = -dat[i]; 
#else
		rr[i] = sf_cneg(dat[i]);
#endif
	    }
	    nloper (false, true, nx, ny, x, rr);
	}
	if (nlreg != NULL) {
	    nlreg  (false, false, nx, nreg, x, rr+ny); 

	    cblas_csscal(nreg,eps,rr+ny,1);
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
	    err[iter] = cblas_scnrm2(ny, rr, 1);
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
}


void sf_solver (sf_operator oper   /* linear operator */, 
		sf_solverstep solv /* stepping function */, 
		int nx             /* size of x */, 
		int ny             /* size of dat */, 
		float* x           /* estimated model */, 
		const float* dat   /* data */, 
		int niter          /* number of iterations */, 
		...                /* variable number of arguments */)
/*< Generic linear solver.
  ---
  Solves
  oper{x}    =~ dat
  ---
  The last parameter in the call to this function should be "end".
  Example: 
  ---
  sf_solver (oper_lop,sf_cgstep,nx,ny,x,y,100,"x0",x0,"end");
  ---
  Parameters in ...:
  ---
  "wt":     float*:         weight      
  "wght":   sf_weight wght: weighting function
  "x0":     float*:         initial model
  "nloper": sf_operator:    nonlinear operator
  "mwt":    float*:         model weight
  "verb":   bool:           verbosity flag
  "known":  bool*:          known model mask
  "nmem":   int:            iteration memory
  "nfreq":  int:            periodic restart
  "xmov":   float**:        model iteration
  "rmov":   float**:        residual iteration
  "err":    float*:         final error
  "res":    float*:         final residual
  >*/ 
{

    va_list args;
    char* par;
    float* wt = NULL;
    sf_weight wght = NULL;
    float* x0 = NULL;
    sf_operator nloper = NULL;
    float* mwt = NULL;
    bool verb = false;
    bool* known = NULL;
    int nmem = -1;
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
	if      (0 == strcmp (par,"end")) {break;}
	else if (0 == strcmp (par,"wt"))      
	{                    wt = va_arg (args, float*);}
	else if (0 == strcmp (par,"wght"))      
	{                    wght = va_arg (args, sf_weight);}
	else if (0 == strcmp (par,"x0"))      
	{                    x0 = va_arg (args, float*);}
	else if (0 == strcmp (par,"nloper"))      
	{                    nloper = va_arg (args, sf_operator);}
	else if (0 == strcmp (par,"mwt"))      
	{                    mwt = va_arg (args, float*);}
	else if (0 == strcmp (par,"verb"))      
	{                    verb = (bool) va_arg (args, int);}    
	else if (0 == strcmp (par,"known"))      
	{                    known = va_arg (args, bool*);}  
	else if (0 == strcmp (par,"nmem"))      
	{                    nmem = va_arg (args, int);}
	else if (0 == strcmp (par,"nfreq"))      
	{                    nfreq = va_arg (args, int);}
	else if (0 == strcmp (par,"xmov"))      
	{                    xmov = va_arg (args, float**);}
	else if (0 == strcmp (par,"rmov"))      
	{                    rmov = va_arg (args, float**);}
	else if (0 == strcmp (par,"err"))      
	{                    err = va_arg (args, float*);}
	else if (0 == strcmp (par,"res"))      
	{                    res = va_arg (args, float*);}
	else 
	{ sf_error("solver: unknown argument %s",par);}
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
	if (mwt != NULL) {
	    for (i=0; i < nx; i++) {
		x[i] *= mwt[i];
	    }
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

    dpr0 = cblas_snrm2(ny, rr, 1);
    dpg0 = 1.;

    for (iter=0; iter < niter; iter++) {
	if ( nmem >= 0) {  /* restart */
	    forget = (bool) (iter >= nmem);
	}
	if (wght != NULL && forget) {
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

	if (mwt != NULL) {
	    for (i=0; i < nx; i++) {
		g[i] *= mwt[i];
	    }
	}
	if (known != NULL) {
	    for (i=0; i < nx; i++) {
		if (known[i]) g[i] = 0.0;
	    }
	} 

	if (mwt != NULL) {
	    for (i=0; i < nx; i++) {
		g[i] *= mwt[i];
	    }
	}
	
	oper (false, false, nx, ny, g, gg);

	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		gg[i] *= wht[i];
	    }
	}
 
	if (forget && nfreq != 0) { /* periodic restart */
	    forget = (bool) (0 == (iter+1)%nfreq); 
	}


	if (iter == 0) {
	    dpg0  = cblas_snrm2 (nx, g, 1);
	    dpr = 1.;
	    dpg = 1.;
	} else {
	    dpr = cblas_snrm2 (ny, rr, 1)/dpr0;
	    dpg = cblas_snrm2 (nx, g , 1)/dpg0;
	}    

	if (verb) 
	    sf_warning ("iteration %d res %f mod %f grad %f",
			iter+1, dpr, cblas_snrm2 (nx, x, 1), dpg);
	

	if (dpr < TOLERANCE || dpg < TOLERANCE) {
	    if (verb) 
		sf_warning("convergence in %d iterations",iter+1);

	    if (mwt != NULL) {
		for (i=0; i < nx; i++) {
		    x[i] *= mwt[i];
		}
	    }
	    break;
	}

	solv (forget, nx, ny, x, g, rr, gg);
	forget = false;

	if (nloper != NULL) {
	    for (i=0; i < ny; i++) {
		rr[i] = -dat[i]; 
	    }

	    if (mwt != NULL) {
		for (i=0; i < nx; i++) {
		    x[i] *= mwt[i];
		}
	    }
	    nloper (false, true, nx, ny, x, rr);
	} else if (wht != NULL) {
	    for (i=0; i < ny; i++) {
		rr[i] = -dat[i]; 
	    }
	    if (mwt != NULL) {
		for (i=0; i < nx; i++) {
		    x[i] *= mwt[i];
		}
	    }
	    oper (false, true, nx, ny, x, rr);
	}  else if (mwt != NULL && (xmov != NULL || iter == niter-1)) {
	    for (i=0; i < nx; i++) {
		x[i] *= mwt[i];
	    }
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
	    err[iter] = cblas_snrm2(ny, rr, 1);
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

void sf_left_solver (sf_operator oper   /* linear operator */, 
		     sf_solverstep solv /* stepping function */, 
		     int nx             /* size of x and dat */, 
		     float* x           /* estimated model */, 
		     const float* dat   /* data */, 
		     int niter          /* number of iterations */, 
		     ...                /* variable number of arguments */)
/*< Generic linear solver for non-symmetric operators.
  ---
  Solves
  oper{x}    =~ dat
  ---
  The last parameter in the call to this function should be "end".
  Example: 
  ---
  sf_left_solver (oper_lop,sf_cdstep,nx,ny,x,y,100,"x0",x0,"end");
  ---
  Parameters in ...:
  ---
  "wt":     float*:         weight      
  "wght":   sf_weight wght: weighting function
  "x0":     float*:         initial model
  "nloper": sf_operator:    nonlinear operator
  "mwt":    float*:         model weight
  "verb":   bool:           verbosity flag
  "known":  bool*:          known model mask
  "nmem":   int:            iteration memory
  "nfreq":  int:            periodic restart
  "xmov":   float**:        model iteration
  "rmov":   float**:        residual iteration
  "err":    float*:         final error
  "res":    float*:         final residual
  >*/ 
{

    va_list args;
    char* par;
    float* wt = NULL;
    sf_weight wght = NULL;
    float* x0 = NULL;
    sf_operator nloper = NULL;
    float* mwt = NULL;
    bool verb = false;
    bool* known = NULL;
    int nmem = -1;
    int nfreq = 0;
    float** xmov = NULL;
    float** rmov = NULL;
    float* err = NULL;
    float* res = NULL;
    float* wht = NULL;
    float *g, *rr, *gg;
    float dpr, dpg, dpr0, dpg0;
    int i, iter; 
    bool forget = false;

    va_start (args, niter);
    for (;;) {
	par = va_arg (args, char *);
	if      (0 == strcmp (par,"end")) {break;}
	else if (0 == strcmp (par,"wt"))      
	{                    wt = va_arg (args, float*);}
	else if (0 == strcmp (par,"wght"))      
	{                    wght = va_arg (args, sf_weight);}
	else if (0 == strcmp (par,"x0"))      
	{                    x0 = va_arg (args, float*);}
	else if (0 == strcmp (par,"nloper"))      
	{                    nloper = va_arg (args, sf_operator);}
	else if (0 == strcmp (par,"mwt"))      
	{                    mwt = va_arg (args, float*);}
	else if (0 == strcmp (par,"verb"))      
	{                    verb = (bool) va_arg (args, int);}    
	else if (0 == strcmp (par,"known"))      
	{                    known = va_arg (args, bool*);}  
	else if (0 == strcmp (par,"nmem"))      
	{                    nmem = va_arg (args, int);}
	else if (0 == strcmp (par,"nfreq"))      
	{                    nfreq = va_arg (args, int);}
	else if (0 == strcmp (par,"xmov"))      
	{                    xmov = va_arg (args, float**);}
	else if (0 == strcmp (par,"rmov"))      
	{                    rmov = va_arg (args, float**);}
	else if (0 == strcmp (par,"err"))      
	{                    err = va_arg (args, float*);}
	else if (0 == strcmp (par,"res"))      
	{                    res = va_arg (args, float*);}
	else 
	{ sf_error("solver: unknown argument %s",par);}
    }
    va_end (args);
 
    g =  sf_floatalloc (nx);
    rr = sf_floatalloc (nx);
    gg = sf_floatalloc (nx);

    if (wt != NULL || wght != NULL) {
	if (wt != NULL) {
	    wht = wt;
	} else {
	    wht = sf_floatalloc (nx);
	    for (i=0; i < nx; i++) {
		wht[i] = 1.0;
	    }
	} 
    }

    for (i=0; i < nx; i++) {
	rr[i] = - dat[i];
    }
    if (x0 != NULL) {
	for (i=0; i < nx; i++) {
	    x[i] = x0[i];
	} 	
	if (mwt != NULL) {
	    for (i=0; i < nx; i++) {
		x[i] *= mwt[i];
	    }
	} 
	if (nloper != NULL) {
	    nloper (false, true, nx, nx, x, rr);
	} else {
	    oper (false, true, nx, nx, x, rr);
            
	}
    } else {
	for (i=0; i < nx; i++) {
	    x[i] = 0.0;
	} 
    }

    dpr0 = cblas_snrm2(nx, rr, 1);
    dpg0 = 1.;

    for (iter=0; iter < niter; iter++) {
	if ( nmem >= 0) {  /* restart */
	    forget = (bool) (iter >= nmem);
	}
	if (wght != NULL && forget) {
	    wght (nx, rr, wht);
	}
	if (wht != NULL) {
	    for (i=0; i < nx; i++) {
		rr[i] *= wht[i];
		g[i] = rr[i]*wht[i];
	    }
	} else {
	    for (i=0; i < nx; i++) {
		g[i] = rr[i];
	    }
	} 

	if (mwt != NULL) {
	    for (i=0; i < nx; i++) {
		g[i] *= mwt[i];
	    }
	}
	if (known != NULL) {
	    for (i=0; i < nx; i++) {
		if (known[i]) g[i] = 0.0;
	    }
	} 

	if (mwt != NULL) {
	    for (i=0; i < nx; i++) {
		g[i] *= mwt[i];
	    }
	}
	
	oper (false, false, nx, nx, g, gg);

	if (wht != NULL) {
	    for (i=0; i < nx; i++) {
		gg[i] *= wht[i];
	    }
	}
 
	if (forget && nfreq != 0) { /* periodic restart */
	    forget = (bool) (0 == (iter+1)%nfreq); 
	}


	if (iter == 0) {
	    dpg0  = cblas_snrm2 (nx, g, 1);
	    dpr = 1.;
	    dpg = 1.;
	} else {
	    dpr = cblas_snrm2 (nx, rr, 1)/dpr0;
	    dpg = cblas_snrm2 (nx, g , 1)/dpg0;
	}    

	if (verb) 
	    sf_warning ("iteration %d res %f mod %f grad %f",
			iter+1, dpr, cblas_snrm2 (nx, x, 1), dpg);
	

	if (dpr < TOLERANCE || dpg < TOLERANCE) {
	    if (verb) 
		sf_warning("convergence in %d iterations",iter+1);

	    if (mwt != NULL) {
		for (i=0; i < nx; i++) {
		    x[i] *= mwt[i];
		}
	    }
	    break;
	}

	solv (forget, nx, nx, x, g, rr, gg);
	forget = false;

	if (nloper != NULL) {
	    for (i=0; i < nx; i++) {
		rr[i] = -dat[i]; 
	    }

	    if (mwt != NULL) {
		for (i=0; i < nx; i++) {
		    x[i] *= mwt[i];
		}
	    }
	    nloper (false, true, nx, nx, x, rr);
	} else if (wht != NULL) {
	    for (i=0; i < nx; i++) {
		rr[i] = -dat[i]; 
	    }
	    if (mwt != NULL) {
		for (i=0; i < nx; i++) {
		    x[i] *= mwt[i];
		}
	    }
	    oper (false, true, nx, nx, x, rr);
	}  else if (mwt != NULL && (xmov != NULL || iter == niter-1)) {
	    for (i=0; i < nx; i++) {
		x[i] *= mwt[i];
	    }
	}

	if (xmov != NULL) {
	    for (i=0; i < nx; i++) {
		xmov[iter][i] =  x[i];
	    }
	}
	if (rmov != NULL) {
	    for (i=0; i < nx; i++) {
		rmov[iter][i] =  rr[i];
	    }
	}
    
	if (err != NULL) {
	    err[iter] = cblas_snrm2(nx, rr, 1);
	}
    }

    if (res != NULL) {
	for (i=0; i < nx; i++) {
	    res[i] = rr[i];
	}
    }
  
    free (g);
    free (rr);
    free (gg);

    if (wht != NULL && wt == NULL) free (wht);
}


void sf_csolver (sf_coperator oper        /* linear operator */, 
		 sf_csolverstep solv      /* stepping function */, 
		 int nx                   /* size of x */, 
		 int ny                   /* size of dat */, 
		 sf_complex* x            /* estimated model */, 
		 const sf_complex* dat    /* data */, 
		 int niter                /* number of iterations */, 
		 ...                      /* variable number of arguments */) 
/*< Generic linear solver for complex data.
  ---
  Solves
  oper{x}    =~ dat
  ---
  The last parameter in the call to this function should be "end".
  Example: 
  ---
  sf_csolver (oper_lop,sf_cgstep,nx,ny,x,y,100,"x0",x0,"end");
  ---
  Parameters in ...:
  ---
  "wt":     float*:          weight      
  "wght":   sf_cweight wght: weighting function
  "x0":     sf_complex*:  initial model
  "nloper": sf_coperator:    nonlinear operator  
  "verb":   bool:            verbosity flag
  "known":  bool*:           known model mask
  "nmem":   int:             iteration memory
  "nfreq":  int:             periodic restart
  "xmov":   sf_complex**: model iteration
  "rmov":   sf_complex**: residual iteration
  "err":    float*:  final error
  "res":    sf_complex*:  final residual
  >*/ 
{

    va_list args;
    char* par;
    float * wt = NULL;
    sf_cweight wght = NULL;
    sf_complex * x0 = NULL;
    sf_coperator nloper = NULL;
    bool verb = false;
    bool* known = NULL;
    int nmem = -1;
    int nfreq = 0;
    sf_complex ** xmov = NULL;
    sf_complex ** rmov = NULL;
    float * err = NULL;
    sf_complex * res = NULL;
    float * wht = NULL;
    sf_complex *g, *rr, *gg, *td = NULL;
    float dpr, dpg, dpr0, dpg0;
    int i, iter; 
    bool forget = false;

    va_start (args, niter);
    for (;;) {
	par = va_arg (args, char *);
	if      (0 == strcmp (par,"end")) {break;}
	else if (0 == strcmp (par,"wt"))      
	{                    wt = va_arg (args, float*);}
	else if (0 == strcmp (par,"wght"))      
	{                    wght = va_arg (args, sf_cweight);}
	else if (0 == strcmp (par,"x0"))      
	{                    x0 = va_arg (args, sf_complex *);}
	else if (0 == strcmp (par,"nloper"))      
	{                    nloper = va_arg (args, sf_coperator);}
	else if (0 == strcmp (par,"verb"))      
	{                    verb = (bool) va_arg (args, int);}    
	else if (0 == strcmp (par,"known"))      
	{                    known = va_arg (args, bool*);}  
	else if (0 == strcmp (par,"nmem"))      
	{                    nmem = va_arg (args, int);}
	else if (0 == strcmp (par,"nfreq"))      
	{                    nfreq = va_arg (args, int);}
	else if (0 == strcmp (par,"xmov"))      
	{                    xmov = va_arg (args, sf_complex **);}
	else if (0 == strcmp (par,"rmov"))      
	{                    rmov = va_arg (args, sf_complex **);}
	else if (0 == strcmp (par,"err"))      
	{                    err = va_arg (args, float *);}
	else if (0 == strcmp (par,"res"))      
	{                    res = va_arg (args, sf_complex *);}
	else 
	{ sf_error("solver: unknown argument %s",par);}
    }
    va_end (args);
 
    g =  sf_complexalloc (nx);
    rr = sf_complexalloc (ny);
    gg = sf_complexalloc (ny);

    if (wt != NULL || wght != NULL) {
	td = sf_complexalloc (ny);
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
#ifdef SF_HAS_COMPLEX_H
	rr[i] = - dat[i];
#else 
	rr[i] = sf_cneg(dat[i]);
#endif
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
	    x[i] = sf_cmplx(0.0,0.0);
	} 
    }

    dpr0 = cblas_scnrm2(ny, rr, 1);
    dpg0 = 1.;

    for (iter=0; iter < niter; iter++) {
	if ( nmem >= 0) {  /* restart */
	    forget = (bool) (iter >= nmem);
	}
	if (wght != NULL && forget) {
	    wght (ny, rr, wht);
	}
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
		rr[i] *= wht[i];
		td[i] = rr[i]*wht[i];
#else
		rr[i] = sf_crmul(rr[i],wht[i]);
		td[i] = sf_crmul(rr[i],wht[i]);
#endif
	    }
      
	    oper (true, false, nx, ny, g, td);
	} else {
	    oper (true, false, nx, ny, g, rr);
	} 
	if (known != NULL) {
	    for (i=0; i < nx; i++) {
		if (known[i]) g[i] = sf_cmplx(0.0,0.0);
	    }
	} 
	oper (false, false, nx, ny, g, gg);
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
		gg[i] *= wht[i];
#else
		gg[i] = sf_crmul(gg[i],wht[i]);
#endif
	    }
	}
 
	if (forget && nfreq != 0) { /* periodic restart */
	    forget = (bool) (0 == (iter+1)%nfreq); 
	}

	if (iter == 0) {
	    dpg0  = cblas_scnrm2 (nx, g, 1);
	    dpr = 1.;
	    dpg = 1.;
	} else {
	    dpr = cblas_scnrm2 (ny, rr, 1)/dpr0;
	    dpg = cblas_scnrm2 (nx, g , 1)/dpg0;
	}    

	if (verb) {
	    sf_warning ("iteration %d res %f mod %f grad %f",
			iter+1, dpr, cblas_scnrm2 (nx, x, 1), dpg);
	}

	if (dpr < TOLERANCE || dpg < TOLERANCE) {
	    if (verb) 
		sf_warning("convergence in %d iterations",iter+1);
	    break;
	}

	solv (forget, nx, ny, x, g, rr, gg);
	forget = false;

	if (nloper != NULL) {
	    for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
		rr[i] = -dat[i]; 
#else
		rr[i] = sf_cneg(dat[i]);
#endif
	    }
	    nloper (false, true, nx, ny, x, rr);
	}
	if (wht != NULL) {
	    for (i=0; i < ny; i++) {
#ifdef SF_HAS_COMPLEX_H
		rr[i] = -dat[i]; 
#else
		rr[i] = sf_cneg(dat[i]);
#endif
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
	    err[iter] = cblas_scnrm2(ny, rr, 1);
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
  
/* 	$Id: bigsolver.c 7717 2011-10-04 04:29:03Z psava $	 */
