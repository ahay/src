#ifndef _parcel_h
#define _parcel_h

void parcel_init(int dim, int *npatch, int *nwall, int *nwind);
void parcel_lop(bool adj, bool add, int n, int mw, 
		float* wall, float* wind);

#endif
