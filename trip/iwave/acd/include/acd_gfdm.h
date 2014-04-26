#ifndef __ACDGFDMFULL__
#define __ACDGFDMFULL__

#include "acd.hh"

/*--- 1st order deriv time step functions --------------------------*/
int acd_tsfm(RDOM * p, RDOM * r, int ia, void * fdpars);

/*--- 1st order deriv adjoint time step functions ------------------*/
int acd_tsam(RDOM * p, RDOM * r, int ia, void * fdpars);

  
#endif
