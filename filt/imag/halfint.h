#ifndef _halfint_h
#define _halfint_h

#include <rsf.h>

void halfint_init (bool adj, bool inv, int n_in, float rho);
void halfint (float* x);
void halfint_close(void);

#endif
