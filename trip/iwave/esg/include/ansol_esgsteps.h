#ifndef __XWAVE_ANSOL_ESG_STEPS__
#define __XWAVE_ANSOL_ESG_STEPS__

#include "ansol_esgn.h"
#include "rdomain.h"

int ansol_HI_esg_ker2d( RDOM * dom, int iarr, void * pars);
int ansol_HI_esg_ker3d( RDOM * dom, int iarr, void * pars);

#endif
