#ifndef _hypotenusei_h
#define _hypotenusei_h

#include <rsf.h>

void hypotenusei_init(int nt1);
void hypotenusei_set(float t0, float dt, float xs);
void hypotenusei_lop(bool adj, bool add, int n1, int n2, float *zz, float *tt);
void hypotenusei_close(void);

#endif
