#ifndef _butter_h
#define _butter_h

void butter_init(int nw_in);
void butter_close(void);
void butter_set(bool low, float cutoff, int npoly, float *num, float *den);
void butter (int nx, int ny, const float *num, const float *den, 
	     const float *xx, float *yy);

#endif
