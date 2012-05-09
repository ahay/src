#ifndef __IWAVE_SGNPARS__
#define __IWAVE_SGNPARS__

#include "utils.h"

/*----------------------------------------------------------------------------
 * Parameters for time step function.
 *----------------------------------------------------------------------------
 */
typedef struct {
  ireal dt;      // time step - copied from IMODEL.tsinfo
  RPNT lam;      // courant params
  int k;         // scheme order
  int ndim;      // dimension, copied from IMODEL.grid
  IPNT lbc;      // flag left boundary conditions
  IPNT rbc;      // flag right boundary conditions
  ireal * ep0_p;     // precomputed pml arrays
  ireal * ep0_pp;    // p = (1-eta*dt^2)
  ireal * ev0_p;     // pp = p/(1+eta*dt^2)
  ireal * ev0_pp;
  ireal * ep1_p;
  ireal * ep1_pp;   
  ireal * ev1_p;
  ireal * ev1_pp;
  ireal * ep2_p;
  ireal * ep2_pp;   
  ireal * ev2_p;
  ireal * ev2_pp;
  // slightly redundant dim info for aux arrays
  int nep0;
  int nev0;
  int nep1;
  int nev1;
  int nep2;
  int nev2;
} SGN_TS_PARS;  

#endif
