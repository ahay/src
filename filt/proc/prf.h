#ifndef _prf_h
#define _prf_h

float prf_burg(int n, float* trace);
float prf_burg2(int n1, int n2, float** trace);
void prf_define (int n, float a, float eps, float* diag, float** offd);

#endif
