#ifndef _irls_h
#define _irls_h

void irls_init(int n);
void irls_close(void);
void l1 (int n, const float *res, float *weight);
void cauchy (int n, const float *res, float *weight);

#endif
