#ifndef _bigsolver_h
#define _bigsolver_h

#include <rsf.h>

typedef void (*operator)(bool,bool,int,int,float*,float*);
typedef void (*solverstep)(bool,int,int,float*,
			   const float*,float*,const float*);

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
void solver_prec (operator oper, solverstep solv, operator prec, 
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
void solver_reg (operator oper, solverstep solv, operator reg, 
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
void solver (operator oper, solverstep solv, int nx, int ny, 
	     float* x, const float* dat, int niter, ...); 

#endif
