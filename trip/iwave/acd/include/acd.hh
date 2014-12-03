#ifndef __IWAVE_ACD_INFO
#define __IWAVE_ACD_INFO

//#include "iwinfo.hh"
#include "iwave.h"
#include "except.hh"
#include "acd_gfdm.h"
#include "acd_gfdm2.h"

int acd_modelinit(PARARRAY *pars,
		  FILE *stream,
		  grid const & g,
		  ireal dt,
		  void ** specs);

void acd_modeldest(void ** specs);

int acd_timegrid(PARARRAY * pars, 
		 FILE * stream, 
		 grid const & g, 
		 ireal & dt);

void acd_timestep(std::vector<RDOM *> dom, 
		  bool fwd, 
		  int iv, 
		  void* fdpars);

int acd_create_sten(void *, 
		    FILE *, 
		    int, 
		    IPNT[RDOM_MAX_NARR], 
		    STENCIL *);

void acd_check(RDOM * dom,
	       void * specs,
	       FILE * stream);

typedef struct {
    ireal dt;      /* time step - copied from IMODEL.tsinfo */
    RPNT lam;      /* courant params */
    int k;         /* scheme order */
    int ndim;      /* dimension, copied from IMODEL.grid */
    IPNT lbc;      /* flag left boundary conditions */
    IPNT rbc;      /* flag right boundary conditions */
    /* coefficient arrays for FD schemes - set in readschemeinfo as they
    // are data-dependent
    // encoding (as RPNT): c[diff index][half-order]
    // thus for 3D, c23[0] is the coefficient of u[.][.][.+2]-u[.][.][.-1]  
    // or u[.][.][.+1]-u[.][.][.-2] in a sixth-order
    // staggered difference formula */
    ireal c0;
    RPNT c1;
    RPNT c2;
    RPNT c3; 
    RPNT c4;
  // test params
  ireal cmax;
  ireal cmin;
} ACD_TS_PARS;  

/*
  Indices of arrays - const density acoustic. 
  D_UC  :  current acoustic potential
  D_UP  :  past acoustic potential
  D_CSQ :  square velocity (multiplier of Laplacian)
*/
#define D_UC  1
#define D_UP  2
#define D_CSQ 0

/* default fraction of max time step */
#define CFL_DEF 0.95

#endif


