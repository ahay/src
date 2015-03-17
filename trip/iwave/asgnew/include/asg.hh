#ifndef __ASG__
#define __ASG__

#include "iwave.h"
#include "except.hh"
//#include "sgcoeffs.h"

/*
int asg_modelinit(PARARRAY *pars,
		  FILE *stream,
		  grid const & g,
		  ireal dt,
		  std::vector<std::string> & active,
		  void ** specs);
*/
int asg_modelinit(PARARRAY pars,
		  FILE *stream,
		  IMODEL & model);

void asg_modeldest(void ** specs);

int asg_timegrid(PARARRAY * pars, 
		 FILE * stream, 
		 grid const & g, 
		 ireal & dt,
		 ireal & rhs);

void asg_timestep(std::vector<RDOM *> dom, 
		  bool fwd, 
		  int iv, 
		  void* fdpars);

int asg_create_sten(void *, 
		    FILE *, 
		    int, 
		    IPNT[RDOM_MAX_NARR], 
		    STENCIL *);

void asg_check(RDOM * dom,
	       void * specs,
	       FILE * stream);

typedef struct s_asg_ts_pars {
  ireal dt;      /* time step - copied from IMODEL.tsinfo */
  RPNT lam;      /* courant params */
  int k;         /* scheme order */
  int ndim;      /* dimension, copied from IMODEL.grid */
  IPNT lbc;      /* flag left boundary conditions */
  IPNT rbc;      /* flag right boundary conditions */
  IPNT nls;      /* number points left pml zone */
  IPNT nrs;      /* number points right pml zone */
  // FD coeffs
  ireal * coeffs[RARR_MAX_NDIM];
  // test params
  ireal cmax;
  ireal cmin;
  // PML
  ireal amp;
  ireal * ep[RARR_MAX_NDIM];     /* precomputed pml arrays */
  ireal * epp[RARR_MAX_NDIM];    /* p = (1-eta*dt^2) */
  ireal * ev[RARR_MAX_NDIM];     /* pp = p/(1+eta*dt^2) */
  ireal * evp[RARR_MAX_NDIM];
} ASG_TS_PARS;  

/* default fraction of max time step */
#define CFL_DEF 0.95

#endif


