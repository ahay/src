#ifndef _coh_h
#define _coh_h

typedef struct Coh *coh;

coh coh_init(int nw, int nj, int nx, int ny, int nz, float ***pp);
float coh1 (const coh ap, int jx, int jy, int jz, int nx, int ny, int nz,
	    float*** xx);
float coh2 (const coh ap, int jx, int jy, int jz, int nx, int ny, int nz,
	    float*** xx);

#endif

/* 	$Id: coh.h,v 1.1 2004/05/25 00:46:12 fomels Exp $	 */
