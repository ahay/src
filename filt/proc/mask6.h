#ifndef _mask6_h
#define _mask6_h

void mask3 (int nw, int nj, int nx, int ny, float **yy, bool **mm);
void mask6 (int nw, int nj1, int nj2, int nx, int ny, float **yy, bool **mm);

void mask32 (int nw, int nj1, int nj2, int nx, int ny, int nz, 
	     float ***yy, bool ***m1, bool ***m2);

#endif
