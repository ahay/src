#ifndef _lapfill_h
#define _lapfill_h

#include <rsf.h>

void lapfill_init (int m1, int m2, bool grad);
void lapfill_close (void);
void lapfill(int niter, float* mm, bool *known);

#endif
