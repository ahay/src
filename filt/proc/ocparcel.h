#ifndef _ocparcel_h
#define _ocparcel_h

#include <stdio.h>

void ocparcel_init(int dim, int *npatch, int *nwall, int *nwind);
void ocparcel_close(void);
void ocparcel_lop(bool adj, int n, int mw, FILE* wall, float* wind);

#endif
