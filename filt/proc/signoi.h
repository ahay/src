#ifndef _signoi_h
#define _signoi_h

#include <rsf.h>

#include "helix.h"

void signoi_init(filter nn_in, filter ss_in, 
		 int niter_in, int nd_in, float eps_in);
void signoi_lop (bool adj, bool add, int n1, int n2, 
		 float *data, float *sign);
void signoi_close(void);

#endif
