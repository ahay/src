#ifndef _twoplane2_h
#define _twoplane2_h

void twoplane2_init (int nw_in, int nj1, int nj2, int nx_in, int ny_in,
		     float **pp_in, float **qq_in);
void twoplane2_close(void);
void twoplane2_lop (bool adj, bool add, int n1, int n2, float *xx, float *yy);

#endif
