#ifndef _dmo_h
#define _dmo_h

#include <rsf.h>

void dmo_init (float vel1, bool inv1, float t01, float dt1, float dx1,
	       int nt1, int nx1, float h1, int mint1, int n1, int type1);

void dmo_close(void);

void dmo_lop (bool adj, bool add, int n1, int n2, float *dat1, float *dat2);

#endif
