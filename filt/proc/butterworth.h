#ifndef _butterworth_h
#define _butterworth_h

void bfhighpass (int npoles, float f3db, int n, float p[], float q[]);
void bflowpass (int npoles, float f3db, int n, float p[], float q[]);
void bfdesign (float fpass, float apass, float fstop, float astop,
	int *npoles, float *f3db);

#endif

/* 	$Id: butterworth.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
