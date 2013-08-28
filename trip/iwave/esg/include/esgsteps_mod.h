#include "esgn.h"

#ifdef __ESG_STEPS_MOD__

#include "rdomain.h"

int esg_2d(RDOM *dom, int iarr, void *pars);

int esg_step_s( RDOM * dom, void * pars );
int esg_step_v( RDOM * dom, void * pars );

#endif
