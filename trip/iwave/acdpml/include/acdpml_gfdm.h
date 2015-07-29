#ifndef __ACDPMLGFDM__
#define __ACDPMLGFDM__

#include "acdpml.hh"

/*--- 1st order deriv time step functions --------------------------*/
int acdpml_tsfm(RDOM * p, RDOM * r, int ia, void * fdpars);

/*--- 1st order deriv adjoint time step functions ------------------*/
int acdpml_tsam(RDOM * p, RDOM * r, int ia, void * fdpars);

  
#endif
