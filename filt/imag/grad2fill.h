#ifndef _grad2fill_h
#define _grad2fill_h

#include <stdbool.h>

void grad2fill_init (int m1, int m2);
void grad2fill_close (void);
void grad2fill(int niter, float* mm, bool *known);

#endif
