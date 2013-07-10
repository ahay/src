#ifndef __XWAVE_ESG_STEPS__
#define __XWAVE_ESG_STEPS__

#include "rdomain.h"

/* generalized time step functions for acoustic staggered grid simulation */
/* int esg_ts22(RDOM * unext, RDOM * ucurr, RDOM * coef, int iarr, void *pars); */
/* int esg_ts22_adj(RDOM * unext, RDOM * ucurr, RDOM * coef, int iarr, void *pars); */
/* int esg_ts24(RDOM * unext, RDOM * ucurr, RDOM * coef, int iarr, void *pars); */
/* int esg_ts24_adj(RDOM * unext, RDOM * ucurr, RDOM * coef, int iarr, void *pars); */
/* int esg_ts2k(RDOM * unext, RDOM * ucurr, RDOM * coef, int iarr, void *pars); */
/* int esg_ts2k_adj(RDOM * unext, RDOM * ucurr, RDOM * coef, int iarr, void *pars); */
int esgts_fwd(RDOM * u, int iarr, void *pars);
int esgts_adj(RDOM * u, int iarr, void *pars);
int esgtsm_fwd(RDOM * unext, RDOM * ucurr, RDOM * coef, int iarr, void *pars);
int esgtsm_adj(RDOM * unext, RDOM * ucurr, RDOM * coef, int iarr, void *pars);

#endif
