#ifndef _neighbors_h
#define _neighbors_h

void neighbors_init (int *in, float *rdx, int *n, 
		     int order, float *time);
int  neighbours(int i);
int nearsource(float* xs, int* b, float* d, float* vv);

#endif
