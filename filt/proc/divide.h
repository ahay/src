#ifndef _divide_h
#define _divide_h

void divide_init(int n1, int n2, float f1, float f2, int niter1, bool gauss);
void divide_close (void);
void divide (const float* num, float* den,  float* rat);

#endif

/* 	$Id: divide.h,v 1.3 2004/02/27 20:59:57 fomels Exp $	 */

