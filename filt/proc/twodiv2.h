#ifndef _twodiv2_h
#define _twodiv2_h

void twodiv2_init(int nw, int n1, int n2, float f1, float f2, int niter1, 
		  bool gauss, float *den);
void twodiv2_close (void);
void twodiv2 (float* num, float* rat);

#endif

/* 	$Id$	 */

