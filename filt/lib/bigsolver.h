#ifndef _sf_bigsolver_h
#define _sf_bigsolver_h

#include "c99.h"

typedef void (*sf_operator)(bool,bool,int,int,float*,float*);
typedef void (*sf_solverstep)(bool,int,int,float*,
			   const float*,float*,const float*);

typedef void (*sf_coperator)(bool,bool,int,int,float complex*,float complex*);
typedef void (*sf_csolverstep)(bool,int,int,float complex*,
			       const float complex*,float complex*,
			       const float complex*);

typedef void (*sf_weight)(int,const float*,float*);
typedef void (*sf_cweight)(int,const float complex*,float*);

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
   nprec - size of p (preconditioned model)
   nx - model size
   ny - data size
   x[nx] - model
   dat[ny] - data
   niter - number of iterations
   eps - scaling factor */
void sf_solver_prec (sf_operator oper, sf_solverstep solv, sf_operator prec, 
		     int nprec, int nx, int ny, float* x, const float* dat, 
		     int niter, float eps, ...);

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
void sf_solver_reg (sf_operator oper, sf_solverstep solv, sf_operator reg, 
		    int nreg, int nx, int ny, float* x, const float* dat, 
		    int niter, float eps, ...);

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
void sf_solver (sf_operator oper, sf_solverstep solv, int nx, int ny, 
		float* x, const float* dat, int niter, ...); 

/* csolver
   ------
   Generic linear solver for complex operators.
   Solves
   oper{x} =~ dat

   oper - operator function
   solv - solver step
   nx - model size
   ny - data size
   x[nx] - model
   dat[ny] - data
   niter - number of iterations */
void sf_csolver (sf_coperator oper, sf_csolverstep solv, int nx, int ny, 
		 float complex* x, const float complex* dat, int niter, ...); 


#endif

/* 	$Id: bigsolver.h,v 1.1 2004/06/11 10:46:55 fomels Exp $	 */
