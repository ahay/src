#ifndef _fdmig_h
#define _fdmig_h

#include <rsf.h>

void fdmig_init(bool hi1, int nx1, int nz1, int nw1, 
		float dx, float dz, float dw1, float vel, float beta);
void fdmig_close(void);
void fdmig (float complex **dat, float **img, sf_file movie);

#endif
