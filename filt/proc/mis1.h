#ifndef _mis1_h
#define _mis1_h

void mis1_init(int n1, int na, float *aa);
void mis1_close(void);
void mis1(int niter, float *xx, const bool *known);

#endif
