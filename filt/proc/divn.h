#ifndef _divn_h
#define _divn_h

void divn_init(int ndim, int nd, int *ndat, int *nbox, int niter);
void divn_close (void);
void divn (float* num, float* den,  float* rat);

#endif

/* 	$Id: divn.h,v 1.1 2004/04/19 21:55:10 fomels Exp $	 */

