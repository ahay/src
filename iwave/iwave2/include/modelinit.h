#ifndef __IWAVE_MODELINIT__
#define __IWAVE_MODELINIT__

#include "fd.h"

int im_init(IMODEL * model, 
	    int (*fdinit)(PARARRAY * pars, FILE * stream, FD_MODEL * specs, IMODEL * model), 
	    PARARRAY * par, 
	    FILE * stream);

#endif
