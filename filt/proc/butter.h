#ifndef _butter_h
#define _butter_h

void butter_set(bool low, float cutoff, int na, float *num, float *den);
void butter (bool adj, int na, float *num, float *den, 
	     int nx, int ny, float *xx, float *yy);
#endif
