#ifndef _butter_h
#define _butter_h

#include <rsf.h>

typedef struct Butter *butter;

butter butter_init(bool low, float cutoff, int nn);
void butter_close(butter bw);
void butter_apply (const butter bw, int nx, float *x);

#endif
