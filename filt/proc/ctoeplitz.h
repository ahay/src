#ifndef _ctorplitz_h
#define _ctorplitz_h

#include <rsf.h>

void ctoeplitz_init (int n);
void ctoeplitz_solve (const float complex *r, float complex *f);
void ctoeplitz_close();

#endif
