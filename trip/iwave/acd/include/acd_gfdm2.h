#ifndef __ACDGFDM2__
#define __ACDGFDM2__

#include "acd.hh"

/*--- 2nd order deriv time step functions --------------------------*/
int acd_tsfm2(RDOM * dd, RDOM *d0, RDOM * d, RDOM * r, int ia, void * fdpars);

/*--- 2nd order deriv adjoint time step functions --------------------------*/
int acd_tsam2(RDOM * db, RDOM *b, RDOM * d, RDOM * r, int ia, void * fdpars);

  
#endif
