#ifndef _burg_h
#define _burg_h

float pef_burg(int n, float* trace);
float pef_burg2(int n1, int n2, float** trace);
void pef_define (int n, float a, float eps, float* diag, float** offd);

#endif
