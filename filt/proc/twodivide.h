#ifndef _twodivide_h
#define _twodivide_h

void twodivide_init(int n1, int n2, float f1, float f2, int niter1, 
		    bool gauss);
void twodivide_close (void);
void twodivide (const float* num, float* den,  float* rat);

#endif

/* 	$Id: twodivide.h,v 1.1 2004/02/27 21:07:59 fomels Exp $	 */

