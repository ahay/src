#ifndef _hipass_h
#define _hipass_h

void hipass_init (float eps);
void hipass_lop(bool adj, bool add, int nx, int ny, float *x, float *y);

#endif
