#ifndef _cburg_h
#define _cburg_h

#include <rsf.h>

void cburg_init (int n_in, int nc_in, int nf_in);
void cburg_close(void);
void cburg_apply (float complex **x, float complex *a);

#endif
